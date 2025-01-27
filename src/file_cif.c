/*
Copyright (C) 2003 by Sean David Fleming

sean@ivec.org

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

The GNU GPL can also be found at http://www.gnu.org
*/

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "gdis.h"
#include "coords.h"
#include "file.h"
#include "parse.h"
#include "matrix.h"
#include "model.h"
#include "interface.h"

#define DEBUG_MORE 0
#define MAX_KEYS 15

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/***************/
/* CIF writing */
/***************/
gint write_cif(gchar *filename, struct model_pak *data)
{
gint flag=0;
gdouble tmat[9], x[3], depth=1.0;
GSList *list;
FILE *fp;
time_t t1;
struct core_pak *core;

/* init & checks */
g_return_val_if_fail(data != NULL, 1);
fp = fopen(filename, "wt");
g_return_val_if_fail(fp != NULL, 2);
  
/* is it a surface with negative z? */
if (data->periodic == 2)
  {
  if (g_slist_length(data->cores))
    {
    core = (struct core_pak *) g_slist_nth_data(data->cores, 0);
    if (core->x[2] < 0.0)
      flag++;
    }
  }

/* rot to make z positive */
matrix_rotation(&tmat[0], G_PI, PITCH);

t1 = time(NULL);
fprintf(fp, "data_%s\n", data->basename);
fprintf(fp, "_audit_creation_date         '%s'\n", g_strstrip(ctime(&t1)));
fprintf(fp, "_audit_creation_method       'generated by GDIS v%4.2f'\n", VERSION);
fprintf(fp, "\n\n");
if (data->periodic)
  {
  fprintf(fp, "_cell_length_a            %8.4f\n", data->pbc[0]);
  fprintf(fp, "_cell_length_b            %8.4f\n", data->pbc[1]);

  if (data->periodic == 2)
    {
/* get depth info */
    depth = (data->surface.region[0]+data->surface.region[1])*data->surface.depth;
/* no depth info - make it large enough to fit everything */
    if (depth < POSITION_TOLERANCE)
      depth = 2.0*data->rmax;
    fprintf(fp, "_cell_length_c            %8.6f\n", depth);
    }
  else
    {
    fprintf(fp, "_cell_length_c            %8.6f\n", data->pbc[2]);
    }

    fprintf(fp, "_cell_angle_alpha       %8.4f\n", fabs(R2D*data->pbc[3]));
    fprintf(fp, "_cell_angle_beta        %8.4f\n", fabs(R2D*data->pbc[4]));
    fprintf(fp, "_cell_angle_gamma       %8.4f\n", fabs(R2D*data->pbc[5]));
/* VASP doesn't like negative angles */
    fprintf(fp, "\n\n");
    fprintf(fp, "_symmetry_space_group_name_H-M   '%s", data->sginfo.spacename);

    if (data->sginfo.originchoice > 1)
      fprintf(fp, " Z");
    fprintf(fp, "'\n");

    fprintf(fp, "_symmetry_Int_Tables_number       %d\n", data->sginfo.spacenum);
    fprintf(fp, "\n\n");
    }

/* charges - format section */
fprintf(fp, "loop_\n");
fprintf(fp, "_atom_type_symbol\n");
fprintf(fp, "_atom_type_oxidation_number\n");

/* charges - data section */
data->unique_atom_list=find_unique(ELEMENT,data);
for (list=data->unique_atom_list ; list ; list=g_slist_next(list))
  {
  fprintf(fp,"%-4s %8.4f\n",
     elements[GPOINTER_TO_INT(list->data)].symbol,
     elements[GPOINTER_TO_INT(list->data)].charge);
  }
fprintf(fp, "\n\n");

/* coords - format section */
fprintf(fp, "loop_\n");
fprintf(fp, "_atom_site_label\n");
fprintf(fp, "_atom_site_type_symbol\n");
if (data->periodic)
  {
  fprintf(fp, "_atom_site_fract_x\n");
  fprintf(fp, "_atom_site_fract_y\n");
  fprintf(fp, "_atom_site_fract_z\n");
  }
else
  {
  fprintf(fp, "_atom_site_Cartn_x\n");
  fprintf(fp, "_atom_site_Cartn_y\n");
  fprintf(fp, "_atom_site_Cartn_z\n");
  }

fprintf(fp, "_atom_site_occupancy\n");

/* coords - data section */
for (list=data->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;

/* ignore deleted/hidden atoms */
  if (core->status & (DELETED | HIDDEN))
    continue;

/* only the asymmetric cell if periodic */
  if (data->periodic && !core->primary)
    continue;

/* NB: want fractional if 3D periodic, otherwise cartesian */
  ARR3SET(x, core->x);
/* transformation needed? */
    if (flag)
      vecmat(tmat, x);
/* make the z coordinate "fractional" for a surface */
  if (data->periodic == 2)
    x[2] /= depth;

  fprintf(fp, "%6s %2s %11.6f %11.6f %11.6f", 
          core->atom_label, elements[core->atom_code].symbol, x[0], x[1], x[2]);
  if (core->has_sof)
    fprintf(fp, " %8.6f\n", core->sof); 
  else
    fputs(" 1.0\n", fp);
  }
fprintf(fp, "\n\n");

fclose(fp);
return(0);
}

/***************/
/* CIF loading */
/***************/
#define DEBUG_LOAD_CIF 0
gint read_cif(gchar *filename, struct model_pak *data)
{
gint i, j, n, natom, new, order, pos, keyword, loop_count=0;
gint min_tokens, num_tokens, len, flag;
gint s_pos, l_pos, x_pos, y_pos, z_pos, o_pos;
gint c_pos;
gchar **buff, *tmp, *name=NULL, *line;
gdouble sign, mat[9], off[3];
gboolean comment, mstring = FALSE;
GSList *list=NULL;
GSList *label=NULL, *oxidation=NULL, *clist=NULL;
struct core_pak *core;
FILE *fp;
gchar *symbol;

/* checks */
g_return_val_if_fail(data != NULL, 1);
g_return_val_if_fail(filename != NULL, 2);

fp = fopen(filename, "rt");
if (!fp)
  return(3);

new = 0;
natom = 0;
data->id = -1;
for(;;)
  {
  comment = FALSE;
/* end if EOF */
  line = file_read_line(fp);
  if (!line)
    break;

/* check for multiline strings beginning with semicolon */
  if (strncmp(g_strstrip(line), ";", 1) == 0)
    mstring ^= TRUE;

/* check for comment beginning wiht #*/
  if (strncmp(g_strstrip(line), "#", 1) == 0)
    comment ^= TRUE;

/* search for data */
  if (!mstring && !comment)
    list = get_keywords(line);

  cif_parse:;
  if (list != NULL)
    {
    keyword = GPOINTER_TO_INT(list->data);
    switch(keyword)
      {
/* model labels */
/* FIXME - needs to search for the 1st occurrence of ' or " & then get the string */

/*    case CIF_CHEMICAL_NAME:
        tmp = g_strdup(get_token_pos(line,1));
        for (i=0 ; i<strlen(tmp) ; i++)
          if (*(tmp+i) == '\'')
            *(tmp+i) = ' ';
        name = g_strdup(g_strstrip(tmp));
        #g_free(data->basename);
        #data->basename = g_strdup(g_strstrip(tmp));
        #g_free(tmp);
        break;
      case CIF_MINERAL_NAME:
        tmp = g_strdup(get_token_pos(line,1));
        for (i=0 ; i<strlen(tmp) ; i++)
          if (*(tmp+i) == '\'')
            *(tmp+i) = ' ';
        name = g_strdup(g_strstrip(tmp));
      #  g_free(data->basename);
      #  data->basename = g_strdup(g_strstrip(tmp));
        g_free(tmp);
        break;
*/
/* stopgap model name */
/* candidate new model trigger */
      case CIF_DATA_START:
        buff = tokenize(line, &num_tokens);
        if (num_tokens == 1)
          {
          if (new > 0)
            {
/* if no cell parameters found yet or parameters read incorrectly
   don't create a new model */
            if (!data->periodic || data->periodic > 3)
              {
              break;
	      }
/* NEW - some dodgy models have no atom data - so allow new model */
/* creation if we have (the bare minimum) some new cell parameters */
            else
              {
              new++;
              data = model_new();
# if DEBUG_LOAD_CIF
printf("Found %d atoms. [reset]\n", natom);
#endif
              }
            }
/* alloc new pointer */
          if (data == NULL)
            goto cif_done;
          name = g_strdup(g_strstrip(&line[5]));
          }
          g_strfreev(buff);
#if DEBUG_LOAD_CIF
printf("Start of new model: %s.\n", data->basename);
#endif
        break;

/* model specific data - the *data ptr MUST be allocated */
      case CIF_CELL_A:
	if (data->periodic < 3) /* for safety */
          {
          buff = get_tokens(line, 3);
          data->pbc[0] = str_to_float(*(buff+1));
          data->periodic++;
          g_strfreev(buff);
	  }
        break;
      case CIF_CELL_B:
	if (data->periodic < 3) /* for safety */
          {
          buff = get_tokens(line, 3);
          data->pbc[1] = str_to_float(*(buff+1));
          data->periodic++;
          g_strfreev(buff);
          }
        break;
      case CIF_CELL_C:
	if (data->periodic < 3) /* for safety */
          {
          buff = get_tokens(line, 3);
          data->pbc[2] = str_to_float(*(buff+1));
          data->periodic++;
          g_strfreev(buff);
          }
        break;
      case CIF_CELL_ALPHA:
        buff = get_tokens(line, 3);
        data->pbc[3] = D2R*str_to_float(*(buff+1));
        g_strfreev(buff);
        break;
      case CIF_CELL_BETA:
        buff = get_tokens(line, 3);
        data->pbc[4] = D2R*str_to_float(*(buff+1));
        g_strfreev(buff);
        break;
      case CIF_CELL_GAMMA:
        buff = get_tokens(line, 3);
        data->pbc[5] = D2R*str_to_float(*(buff+1));
        g_strfreev(buff);
        break;
      case CIF_SPACE_NAME:
/* remove the enclosing ' and " characters */
        tmp = g_strdup(get_token_pos(line,1));
        for (i=0 ; i<strlen(tmp) ; i++) 
          {
          if (*(tmp+i) == '\'')
            *(tmp+i) = ' ';
          if (*(tmp+i) == '\"')
            *(tmp+i) = ' ';
          }

/* store the name, stripping spaces */
        data->sginfo.spacename = g_strdup(g_strstrip(tmp));
/* indicate that name should be used in lookup */
        data->sginfo.spacenum = -1;
#if DEBUG_LOAD_CIF
printf("space group: [%s]\n", data->sginfo.spacename);
#endif
        g_free(tmp);
        break;
      case CIF_SPACE_NUM:
        buff = get_tokens(line, 3);
        data->sginfo.spacenum = (gint)g_ascii_strtod(*(buff+1), NULL);
#if DEBUG_LOAD_CIF
printf("space group number: [%d]\n", data->sginfo.spacenum);
#endif
        g_strfreev(buff);

/* If no space group given, or unrecognizable, use number */
       if( !data->sginfo.spacename || g_ascii_strcasecmp("",data->sginfo.spacename) == 0)
         {
         if( data->sginfo.spacenum > 1)
           {
           if( data->sginfo.spacename)
             {
             g_free(data->sginfo.spacename);
             data->sginfo.spacename = g_strdup_printf("%d", data->sginfo.spacenum);
             }
           }
         else
           {
           data->sginfo.spacename =g_strdup("P 1");
           }
         }
       break;

      case CIF_EQUIV_SITE:
        loop_count++;
        break;

/* NEW - reinserted the lost symmetry matrix code */
      case CIF_EQUIV_POS:
/* allocate for order number of pointers (to matrices) */
        data->sginfo.matrix = (gdouble **) g_malloc(sizeof(gdouble *));
        data->sginfo.offset = (gdouble **) g_malloc(sizeof(gdouble *));
        data->sginfo.order = order = 0;
        for(;;)
          {
/* terminate on EOF, */
          g_free(line);
          line = file_read_line(fp);
          if (!line)
            break;
/* blank line, */
          g_strstrip(line);
          if (!strlen(line))
            break;
/* check for multiline strings beginning with semicolon */
          if (strncmp(g_strstrip(line), ";", 1) == 0)
            mstring ^= TRUE;
/* or a new command */
          if (!mstring && !comment)
            list = get_keywords(line);
          if (list)
            goto cif_parse;

/* TODO - make this parsing a subroutine */
          for(i=0 ; i<strlen(line) ; i++)
            if (*(line+i) == '\'' || *(line+i) == ',')
              *(line+i) = ' ';
          g_strstrip(line);
          buff = tokenize(line, &num_tokens);

          n = loop_count;
          while (n < num_tokens-2)
            {
/* FIXME - yet another mess that a linked list would greatly simplify */
/* number of ops */
            data->sginfo.matrix = (gdouble **) g_renew
                                (gdouble *, data->sginfo.matrix , (order+1));
            data->sginfo.offset = (gdouble **) g_renew
                                (gdouble *, data->sginfo.offset , (order+1));
/* actual op */
            *(data->sginfo.matrix+order) = (gdouble *) g_malloc(9*sizeof(gdouble));
            *(data->sginfo.offset+order) = (gdouble *) g_malloc(3*sizeof(gdouble));

#if DEBUG_LOAD_CIF
printf("[%s] [%s] [%s]\n", *(buff+n+0), *(buff+n+1), *(buff+n+2));
#endif
            VEC3SET(&mat[0], 0.0, 0.0, 0.0);
            VEC3SET(&mat[3], 0.0, 0.0, 0.0);
            VEC3SET(&mat[6], 0.0, 0.0, 0.0);
            VEC3SET(&off[0], 0.0, 0.0, 0.0);

            for (i=0 ; i<3 ; i++)
              {
              pos = 0;
              sign = 1.0;
              for (j=0 ; j<strlen(*(buff+i+n)) ; j++)
                {
                switch (g_ascii_tolower(*(*(buff+i+n)+j)))
                  {
                  case '-':
                    sign = -1.0;
                    break;
                  case '+':
                    sign = +1.0;
                    break;
                  case 'x':
                    mat[i*3] = sign*1.0;
                    pos++;
                    break;
                  case 'y':
                    mat[i*3 + 1] = sign*1.0;
                    pos++;
                    break;
                  case 'z':
                    mat[i*3 + 2] = sign*1.0;
                    pos++;
                    break;
/* FIXME - a bit crude */
                  case '/':
                    g_free(line);
                    line = g_strndup(*(buff+i+n)+j-1, 3);

                    off[i] = sign * str_to_float(line);
                    break;
/* TODO: better way to parse? */
                  case '0':
                    off[i] = sign * str_to_float(*buff);
                    break;
                  }
                }
              }
            ARR3SET((*(data->sginfo.matrix+order)+0), &mat[0]);
            ARR3SET((*(data->sginfo.matrix+order)+3), &mat[3]);
            ARR3SET((*(data->sginfo.matrix+order)+6), &mat[6]);
            ARR3SET((*(data->sginfo.offset+order)+0), &off[0]);

#if DEBUG_LOAD_CIF
P3MAT("output: ", *(data->sginfo.matrix+order));
P3VEC("output: ", *(data->sginfo.offset+order));
printf("\n");
#endif
            order++;
            // data->sginfo.order++;
            n += 3;
            }
          g_strfreev(buff);
          }

#if DEBUG_LOAD_CIF
printf("Found %d symmetry matrices.\n", order);
#endif
        data->sginfo.order = order;
        break;

      case CIF_LOOP_START:
        loop_count = 0;
        break;

/* parsing for atom charges. */
      case CIF_ATOM_TYPE:
        s_pos = c_pos = -1;

        while (g_strrstr(line, "_atom_type") != NULL)
          {
          if (g_strrstr(line, "_atom_type_symbol"))
            s_pos = loop_count;
          if (g_strrstr(line, "_atom_type_oxidation_number"))
            c_pos = loop_count;

/* get next line and keyword list */
          loop_count++;

          g_free(line);
          line = file_read_line(fp);
          if (!line)
            goto cif_done;
          }
/* check for multiline strings beginning with semicolon */
        if (strncmp(g_strstrip(line), ";", 1) == 0)
          mstring ^= TRUE;
/* while no new keywords found */
        if (!mstring)
          list = get_keywords(line);

        while (list == NULL)
          {
          buff = tokenize(line, &num_tokens);
          if (num_tokens >= c_pos && num_tokens >= s_pos)
            {
            if (c_pos > -1)
              oxidation = g_slist_prepend(oxidation, g_strdup(*(buff+c_pos)));
            if (c_pos > -1 && s_pos > -1)
              label = g_slist_prepend(label, g_strdup(*(buff+s_pos)));
            }
#if DEBUG_LOAD_CIF
          else
            printf("Not enough tokens found. ('%s')\n", g_strstrip(line));
#endif
/* get next line */
          g_strfreev(buff);

          g_free(line);
          line = file_read_line(fp);
          if (!line)
            goto cif_done;

/* check for multiline strings beginning with semicolon */
          if (strncmp(g_strstrip(line), ";", 1) == 0)
            mstring ^= TRUE;
          if (!mstring)
            list = get_keywords(line);
/* CURRENT - we really want this keyword list parsed again, just in case */
/* a new model trigger (audit_creation_date or data_xxxx) was found */
          if (list)
            goto cif_parse;
          }
        break;

/* parsing for column info ie x,y,z positions frac/cart etc. */
      case CIF_ATOM_SITE:
        s_pos = l_pos = x_pos = y_pos = z_pos = o_pos = -1;
	line = g_strstrip(line);

        while (g_strrstr(line, "_atom_site") != NULL && strncmp(line, "_atom_", 1) == 0)
          {
          if (g_strrstr(line, "_atom_site_type_symbol"))
	    s_pos = loop_count;
          if (g_strrstr(line, "_atom_site_label"))
            l_pos = loop_count;
          if (g_strrstr(line, "_Cartn_x"))
            {
            x_pos = loop_count;
            data->fractional = FALSE;
            }
          if (g_strrstr(line, "_Cartn_y"))
            {
            y_pos = loop_count;
            data->fractional = FALSE;
            }
          if (g_strrstr(line, "_Cartn_z"))
            {
            z_pos = loop_count;
            data->fractional = FALSE;
            }
          if (g_strrstr(line, "_fract_x"))
            {
            x_pos = loop_count;
            data->fractional = TRUE;
            }
          if (g_strrstr(line, "_fract_y"))
            {
            y_pos = loop_count;
            data->fractional = TRUE;
            }
          if (g_strrstr(line, "_fract_z"))
            {
            z_pos = loop_count;
            data->fractional = TRUE;
            }
          if (g_strrstr(line, "_occupancy"))
            o_pos = loop_count;

/* get next line and keyword list */
          loop_count++;

          g_free(line);
          line = file_read_line(fp);

          if (strncmp(g_strstrip(line), ";", 1) == 0)
            mstring ^= TRUE;

          if (!line)
            goto cif_done;
          }

/* either symbol or label can be present & used for atom id purposes */
        if (s_pos < 0)
          s_pos = l_pos;
        if (l_pos < 0)
          l_pos = s_pos;
/* check for minimum data */
        if (s_pos < 0 || x_pos < 0 || y_pos < 0 || z_pos < 0)
          {
#if DEBUG_LOAD_CIF
printf("read_cif() warning: incomplete cif file? [%d:%d:%d:%d]\n",
                                        s_pos, x_pos, y_pos, z_pos);
#endif
          break;
          }

/* the expected number of tokens */
        min_tokens = loop_count;

#if DEBUG_LOAD_CIF
printf(" min tokens: %d\n", min_tokens);
printf("data format: [%d] (%d) - [%d] [%d] [%d] (%d)",
                  s_pos, l_pos, x_pos, y_pos, z_pos, o_pos);
if (data->fractional)
    printf(" (frac)\n");
else
    printf(" (cart)\n");
#endif

/* while no new keywords found */
        natom = 0;
/* check for multiline strings beginning with semicolon */
        if (strncmp(g_strstrip(line), ";", 1) == 0)
          mstring ^= TRUE;
        if (!mstring)
          list = get_keywords(line);

        while (list == NULL)
          {
          buff = tokenize(line, &num_tokens);

#if DEBUG_LOAD_CIF
	  if (num_tokens < min_tokens)
            printf("Need to read next line. ('%s')\n", g_strstrip(line));
#endif
/* NB: cif is stupid - it allows data to continue on the next line, */
/* until it gets its number of tokens */
/* hopefully, this'll allow us to get something, even on short lines */
/* TODO: Cope with multiline entries */
          if (num_tokens > z_pos)
            {
/* NEW - ignore labelled (*) symmetry equiv. atoms */
            flag = TRUE;
            if (l_pos >= 0 && l_pos < num_tokens)
              {
              len = strlen(*(buff+l_pos));
              if ((*(buff+l_pos))[len-1] == '*')
                flag = FALSE;

              if (s_pos < num_tokens)
	        {
                if (elem_symbol_test(*(buff+s_pos)) && flag)
                  core = core_new(*(buff+s_pos), *(buff+l_pos), data);
		else
                  core = core_new(*(buff+s_pos), NULL, data);
		}
              else if (flag)
                core = core_new(*(buff+l_pos), NULL, data);
	      else
                core = core_new("X", NULL, data);

              data->cores = g_slist_prepend(data->cores, core);

              core->x[0] = str_to_float(*(buff+x_pos));
              core->x[1] = str_to_float(*(buff+y_pos));
              core->x[2] = str_to_float(*(buff+z_pos));

              for (clist = label, i=0; clist; clist = g_slist_next(clist), i++)
                {
                symbol = g_strdup(clist->data);
                if (g_ascii_strcasecmp(symbol, *(buff+s_pos)) == 0)
                  {
                  core->charge = g_ascii_strtod(g_slist_nth_data(oxidation,i), NULL);
                  core->lookup_charge = FALSE;
                  }
                }

/* only get occupancy if we're sure we have enough tokens */
              if (o_pos > -1 && num_tokens >= min_tokens)
                {
                core->sof = str_to_float(*(buff+o_pos));
                core->has_sof = TRUE;
                }
              natom++;
              }
            }
#if DEBUG_LOAD_CIF
          else
            printf("Not enough tokens found. ('%s')\n", g_strstrip(line));
#endif
/* get next line */
          g_strfreev(buff);

          g_free(line);
          line = file_read_line(fp);
          if (!line)
            goto cif_done;

/* check for multiline strings beginning with semicolon */
          if (strncmp(g_strstrip(line), ";", 1) == 0)
            mstring ^= TRUE;
          if (!mstring)
            list = get_keywords(line);

/* CURRENT - we really want this keyword list parsed again, just in case */
/* a new model trigger (ie audit_creation_date) was found */
          if (list)
            goto cif_parse;
          }
        break;
      }
    }
  }
  cif_done:;

  g_free(line);

/* yet another hack to support the dodgy CIF standard */
  if (!new && natom)
    {
/* no new model was triggered - but we found atoms, so we'll */
/* assume only one model was in the file & hope for the best */
    new++;
    }

/* if we found a name - assign it now */
  if (name)
    {
    g_free(data->basename);
    data->basename = name;
    name = NULL;
    }
        
  if (strlen(filename) < FILELEN)
    strcpy(data->filename, filename);
  else
    {
    gui_text_show(ERROR, "File name is too long.\n");
    return(-1);
    }
#if DEBUG_LOAD_CIF
printf("Found %d atom%s\n", natom, (natom != 1 ? "s.":"."));
printf("Found %d symmetry matrices.\n", data->sginfo.order);
printf("Found %d model(s)\n", new);
printf("Periodicity is %d\n", data->periodic);
#endif

/* Move identity matrix to first in list */
  j = data->sginfo.order;
  for (i=1; i < j; i++)
    {
    if (matrix_is_identity(*(data->sginfo.matrix+i)) &&
        matrix_is_empty(*(data->sginfo.offset+i), 3) )
      {
       ARR3SET((*(data->sginfo.matrix+i)+0), (*(data->sginfo.matrix)+0));
       ARR3SET((*(data->sginfo.matrix+i)+3), (*(data->sginfo.matrix)+3));
       ARR3SET((*(data->sginfo.matrix+i)+6), (*(data->sginfo.matrix)+6));
       ARR3SET(*(data->sginfo.offset+i), *(data->sginfo.offset));
       matrix_identity((*data->sginfo.matrix));
       VEC3SET(*(data->sginfo.offset), 0.0, 0.0, 0.0);
       }
    }

/* set up for display */
  for (list=sysenv.mal ; list ; list=g_slist_next(list))
    {
    data = list->data;
    if (data->id == -1)
      {
      data->id = CIF;
      data->cores = g_slist_reverse(data->cores);
      model_prep(data);
      }
    }

#if DEBUG_LOAD_CIF
printf("Setting up %d model(s).\n", new);
#endif

/* clean up & exit */
  if (list)
    g_slist_free(list);

return(0);
}

