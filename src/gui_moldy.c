/*
Copyright (C) 2022 by Craig Andrew James Fisher

c_fisher@jfcc.or.jp

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

/**************************************************************************
* GUI_MOLDY.c	Functions for handling Moldy GUI			  *
*									  *
* pot_titles[]		array of column titles for potential window	  *
* moldy_cleanup()       frees widget variables on closing		  *
* change_basename()     allows name of model to be changed		  *
* choose_potential()	select potential name for display		  *
* update_pot_dialog()	redo potentials box				  *
* strain_mask_toggle()  strain mask for constant strain calcs		  *
* auto_toggle() 	toggle auto calculation of cutoff parms		  *
* surf_toggle() 	toggle surface dipole on/off			  *
* strict_toggle() 	toggle strict cutoff on/off			  *
* ewald_toggle() 	Ewald calculation on/off			  *
* real_toggle() 	real space component only on/off		  *
* scale_opt_toggle() 	options for scaling				  * 
* dump_level_toggle()	dump file info toggles				  *
* gtk_read_potlib() 	hook for reading potential library files	  *
* gtk_write_potlib() 	hook for writing potential library files	  *
* gtk_delete_pots() 	hook for removing potential parameters     	  *
* exec_moldy()	 	generate moldy command and submit to task list	  *
* exec_moldy_task() 	generate input/output filenames for moldy job	  *
* proc_moldy_task() 	process output of moldy run			  *
* run_moldy()		controls moldy job tasks			  *
* thermo_keyword() 	interprets thermostat settings			  *
* pressure_keyword() 	interprets barostat settings			  *
* ewald_keyword() 	interprets ewald sum settings			  *
* unit_keyword()	interprets unit entries                           *
* moldy_potential() 	toggles moldy potential types			  *
* moldy_keyword() 	toggles moldy parameter settings		  *
* moldy_unit_dialog()	opens dialog for inputing units explicitly        *
* reset_units()		restore units to default		          *
* Toggles for energy units:						  *
*    moldy_kjmol(), moldy_kcal(), moldy_ev(), moldy_e2a()		  *
* gui_moldy_widget() 	sets up moldy widget inside dialog		  *
* moldy_dialog()        sets up moldy dialog                              *
**************************************************************************/
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

#include "gdis.h"
#include "coords.h"
#include "dialog.h"
#include "edit.h"
#include "file.h"
#include "task.h"
#include "matrix.h"
#include "spatial.h"
#include "parse.h"
#include "interface.h"
#include "gui_shorts.h"
#include "moldy.h"

extern struct sysenv_pak sysenv;
GtkWidget *pot_list;
gint num_pot_rows = 0;

/* Column titles for potential table */
gchar *pot_titles[] = {"Atom 1", "Chg 1", "Atom 2", "Chg 2",
                       "A", "B", "C", "D", "E", "F", "G", "H"};

/* Name and no of parms for known potentials */
const pots_pak potspec[]  = {{"Lennard-Jones", 2},
                             {"Buckingham", 3},
                             {"Mxdorto", 6},
                             {"MCY", 4},
                             {"HIW", 3},
                             {"Generic", 6},
                             {"Morse", 4},
                             {0, 0}};

#define WIDTH 60  /* Width of entry box */

/*********************************/
/* cleanup for moldy dialog      */
/*********************************/
void moldy_cleanup(struct model_pak *model)
{
g_assert(model != NULL);

model->moldy.pot_label = NULL;
/* model->moldy.cutoff_frame = NULL; */
model->moldy.lib_file_entry = NULL;
}

/**********************************************************/
/* determine potential type label for display             */
/**********************************************************/
gchar* choose_potential(gint n)
{
gchar *label;

if( n < 0 )
    label = g_strdup("none");
else
    label = g_strdup(potspec[n].name);

return label;
}

/**********************************************************/
/* write potential parameters to the current clist widget */
/**********************************************************/
#define DEBUG_UPDATE_POTENTIALS_DIALOG 0
void update_pot_dialog(struct model_pak *data)
{
gint i, j, alt, old_rows=0;
gchar *info[MAX_POT_PARMS+4];
GSList *plist=NULL;
GtkStyle *style, *pot_style;
GdkColor alt_colour = {0, 0x0000, 0x6000, 0x8000};
struct pot_pak *pdata;
GString *label;

/* check */
g_return_if_fail(data != NULL);

#if DEBUG_UPDATE_POTENTIALS_DIALOG
printf("Updating potentials for model %d\n", data->id);
#endif

for( i=0; i<MAX_POT_PARMS+4; i++)
   info[i] = g_strdup("");

/* ensure the list is valid */
if (!pot_list)
  return;
if (!GTK_IS_CLIST(pot_list))
  return;
pot_style = gtk_widget_get_style(pot_list);
if (!pot_style)
  return;

/* pot list update */
gtk_clist_freeze((GtkCList *) pot_list);
gtk_clist_clear((GtkCList *) pot_list);

/* alternate row text colour */
alt=0;
style = gtk_style_copy(pot_style);
style->fg[GTK_STATE_NORMAL] = alt_colour;

/* go through the potential list */
plist = data->moldy.pots;
num_pot_rows = 0;

while( plist != NULL )
     {
     pdata = plist->data;
     old_rows = num_pot_rows;

     for( j=0; j<4; j++)
       g_free(info[j]);
     info[0] = g_strdup(pdata->site1_label);
     info[1] = g_strdup_printf("%g", pdata->site1_charge);
     info[2] = g_strdup(pdata->site2_label);
     info[3] = g_strdup_printf("%g", pdata->site2_charge);

     for(j=0; j < MAX_POT_PARMS; j++)
       info[4+j] = g_strdup_printf("%7.4f", pdata->parm[j]);
     if( g_ascii_strcasecmp(info[0],"") != 0 || g_ascii_strcasecmp(info[2],"") != 0)
       {
       gtk_clist_append((GtkCList *) pot_list, &info[0]);
/* FIXME - be careful not to inadvertently free this pointer */
       gtk_clist_set_row_data((GtkCList *) pot_list, num_pot_rows, &info[0]);
       if (alt)
         gtk_clist_set_row_style(GTK_CLIST(pot_list), num_pot_rows, style);

       num_pot_rows++;
       }
     plist = g_slist_next(plist);

     if (num_pot_rows > old_rows)
       alt ^= 1;
     }
/* unfreeze */
gtk_clist_thaw((GtkCList *) pot_list);

/* is there a dialog entry to be updated? */
if (GTK_IS_ENTRY(data->moldy.lib_file_entry))
  gtk_entry_set_text(GTK_ENTRY(data->moldy.lib_file_entry), data->moldy.lib_file);
/* update potential type */
label = g_string_new(NULL);
label->str = choose_potential(data->moldy.pot_type);
gtk_entry_set_text(GTK_ENTRY(data->moldy.pot_label), label->str);
gtk_entry_set_editable(GTK_ENTRY(data->moldy.pot_label), FALSE);
}

/* Calculate strain mask from clicked buttons */
gint strain_mask_toggle(GtkWidget *w, gpointer *obj)
{
gint index;
struct model_pak *data;
index = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(obj), "key"));

data = (struct model_pak *) g_object_get_data(G_OBJECT(obj), "ptr");
g_return_val_if_fail(data != NULL, FALSE);

data->moldy.strain_mask ^= index;

return(FALSE);
}

/* Auto cutoff toggle button */
void auto_toggle(struct model_pak *model)
{
model->moldy.auto_cutoff ^= 1;
/* if (model->moldy.auto_cutoff)
  gtk_widget_set_sensitive(GTK_WIDGET(model->moldy.cutoff_frame), FALSE);
else
  gtk_widget_set_sensitive(GTK_WIDGET(model->moldy.cutoff_frame), TRUE); */
}

/* Strict cutoff toggle button */
void strict_toggle(struct model_pak *data)
{
data->moldy.strict ^= 1;
}

/* Surface dipole toggle button */
void surf_toggle(struct model_pak *data)
{
data->moldy.surf ^= 1;
}

/* Include Ewald Sum toggle button */
void ewald_toggle(struct model_pak *data)
{
data->moldy.ewald ^= 1;
}

/* Include real space toggle button */
void real_toggle(struct model_pak *data)
{
data->moldy.real_only ^= 1;
}

/* Scaling/thermostat options */
gint scale_opt_toggle(GtkWidget *w, gpointer *obj)
{
gint index;
struct model_pak *data;
index = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(obj), "key"));

data = (struct model_pak *) g_object_get_data(G_OBJECT(obj), "ptr");
g_return_val_if_fail(data != NULL, FALSE);

data->moldy.scale_options ^= index;

return(FALSE);
}

gint dump_level_toggle(GtkWidget *w, gpointer *obj)
{
gint index;
struct model_pak *data;
index = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(obj), "key"));

data = (struct model_pak *) g_object_get_data(G_OBJECT(obj), "ptr");
g_return_val_if_fail(data != NULL, FALSE);

data->moldy.dump_level ^= index;

return(FALSE);
}

void gtk_read_potlib(GtkWidget *w, gpointer *obj)
{
struct model_pak *data;

data = (struct model_pak *) g_object_get_data(G_OBJECT(obj), "ptr");
g_return_if_fail(data != NULL);

moldy_load_dialog("Load Moldy potentials file", FILE_LOAD, (gpointer) read_moldy_potentials, data, MOLDY_POT);
}

void gtk_write_potlib(GtkWidget *w, gpointer *obj)
{
struct model_pak *data;

data = (struct model_pak *) g_object_get_data(G_OBJECT(obj), "ptr");
g_return_if_fail(data != NULL);

if( data->moldy.pots )
  moldy_load_dialog("Save Moldy potentials file", FILE_SAVE, (gpointer) write_moldy_potentials, data, MOLDY_POT);
}

void gtk_delete_pots(GtkWidget *w, gpointer *obj)
{
struct model_pak *data;

data = (struct model_pak *) g_object_get_data(G_OBJECT(obj), "ptr");
g_return_if_fail(data != NULL);

delete_moldy_potentials(data);
}

/*void gtk_load_moldy(GtkWidget *w, gpointer *obj )
{
struct model_pak *data;

data = (struct model_pak *) g_object_get_data(G_OBJECT(obj), "ptr");
g_return_if_fail(data != NULL);

moldy_load_dialog("Load Moldy restart file", FILE_LOAD, (gpointer) read_moldy_restart, data, MOLDY_RES);
} */

/********************/
/* run a moldy file */
/********************/
#define DEBUG_EXEC_MOLDY 0
gint exec_moldy(const gchar *input, const gchar *output)
{
gchar *cmd;

/* checks */
if (!sysenv.moldy_path)
  return(-1);

gui_text_show(WARNING, "Moldy calc started.\n");

/* stop moldy renaming (if it exists already) */
unlink(output);

cmd = g_strdup_printf("%s %s %s", sysenv.moldy_path, input, output);

#if DEBUG_EXEC_MOLDY
printf("executing: [%s]\n",cmd);
#endif

task_sync(cmd);

/* done */
g_free(cmd);
return(0);
}
/******************************/
/* execute a moldy run (task) */
/******************************/
#define DEBUG_EXEC_MOLDY_TASK 1
void exec_moldy_task(gpointer *ptr)
{
gchar *control_file, *sysspec_file;
GString *outfile;
gchar **buff;
struct model_pak *data;

g_return_if_fail(ptr != NULL);
data = (struct model_pak *) ptr;

outfile = g_string_new(NULL);

control_file = construct_filename(data->moldy.control_file);
sysspec_file = construct_filename(data->moldy.sysspec_file);

write_moldy_control(control_file, data);

#if DEBUG_EXEC_MOLDY_TASK
printf("control file : %s\n", control_file);
printf("sys-spec file: %s\n", sysspec_file);
#endif

/* Only write sys-spec info if no restart file */
if( !(data->moldy.restart_file) || strlen(data->moldy.restart_file) == 0 ||
    !g_ascii_strcasecmp(data->moldy.restart_file, "none"))
  write_moldy_sys(sysspec_file, data);

if( g_ascii_strcasecmp(data->moldy.out_file, "none") == 0 ||
      g_ascii_strcasecmp(data->moldy.out_file, "") == 0 )
   {
/* copy first token as the base name */
   buff = g_strsplit(data->moldy.sysspec_file, ".min", 0);
   g_string_assign(outfile, *buff);
   g_string_append(outfile, ".mdo");
   g_strfreev(buff);
   }
else 
   g_string_assign(outfile, data->moldy.out_file);

#if DEBUG_EXEC_MOLDY_TASK
printf("output file  : %s\n", outfile->str);
#endif

/* delete the old file to be sure output is only data from current run */
unlink(outfile->str);

/* are we supposed to execute Moldy? */
if (data->moldy.no_exec)
  {
#if DEBUG_EXEC_MOLDY_TASK
printf("Skipping Moldy execution on user request.\n");
#endif
  return;
  }

/* run */
exec_moldy(control_file, outfile->str);
}

/******************************/
/* process a moldy run (task) */
/******************************/
#define DEBUG_PROC_MOLDY 0
void proc_moldy_task(gpointer *ptr)
{
struct model_pak *data;

/* TODO - model locking (model_ptr RO/RW etc) to prevent screw-ups */
g_return_if_fail(ptr != NULL);
data = (struct model_pak *) ptr;

/* don't attempt to process if moldy execution was turned off */
if (data->moldy.no_exec)
  return;

gui_text_show(WARNING, "Moldy calc completed.\n");
/* TODO - process output from Moldy */

return;
}

/************************************************/
/* controls the use of MOLDY MD runs            */
/************************************************/
#define DEBUG_RUN_MOLDY 0
void run_moldy(GtkWidget *w)
{
gchar **buff;
struct model_pak *data;

data = sysenv.active_model;

/* checks */
g_return_if_fail(data != NULL);
if (!sysenv.moldy_path && !data->moldy.no_exec)
  {
  gui_text_show(ERROR, "MOLDY executable was not found.\n");
  return;
  }

if( g_ascii_strcasecmp(data->moldy.out_file, "none") == 0 || 
       g_ascii_strcasecmp(data->moldy.out_file, "") == 0 )
   {
/* copy first token as the base name */
   g_free(data->moldy.out_file);
   buff = g_strsplit(data->moldy.sysspec_file, ".min", 0);
   data->moldy.out_file = g_strconcat(g_strdup(*buff), ".mdo", NULL);
   g_strfreev(buff);
   }

#if DEBUG_RUN_MOLDY
printf("output file: %s\n", data->moldy.out_file);
#endif

/* run as a background process */
task_new("Moldy", &exec_moldy_task, data, &proc_moldy_task, data, data);

}

gint thermo_keyword(GtkWidget *w, gpointer *obj)
{
gint keyword;
struct model_pak *data;

keyword = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(obj), "key"));

data = (struct model_pak *) g_object_get_data(G_OBJECT(obj), "ptr");
g_return_val_if_fail(data != NULL, FALSE);

data->moldy.const_temp = keyword;

return(FALSE);
}
gint pressure_keyword(GtkWidget *w, gpointer *obj)
{
gint keyword;
struct model_pak *data;

keyword = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(obj), "key"));

data = (struct model_pak *) g_object_get_data(G_OBJECT(obj), "ptr");
g_return_val_if_fail(data != NULL, FALSE);

data->moldy.const_press = keyword;

return(FALSE);
}

gint ewald_keyword(GtkWidget *w, gpointer *obj)
{
gint keyword;
struct model_pak *data;

keyword = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(obj), "key"));

data = (struct model_pak *) g_object_get_data(G_OBJECT(obj), "ptr");
g_return_val_if_fail(data != NULL, FALSE);

switch(keyword)
  {
  case EWALD_INC:
    data->moldy.ewald ^= 1;
    break;
  case REAL_INC:
    data->moldy.real_only ^= 1;
    break;
  }

return(FALSE);
}

gint unit_keyword(GtkWidget *w, gpointer *obj)
{
gint keyword;
struct model_pak *data;

keyword = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(obj), "key"));

data = (struct model_pak *) g_object_get_data(G_OBJECT(obj), "ptr");
g_return_val_if_fail(data != NULL, FALSE);

switch(keyword)
  {
  case TIME_UNIT:
    data->moldy.time_unit = str_to_float(gtk_entry_get_text
             (GTK_ENTRY(obj)));
    break;
  case LENGTH_UNIT:
    data->moldy.length_unit = str_to_float(gtk_entry_get_text
             (GTK_ENTRY(obj)));
    break;
  case MASS_UNIT:
    data->moldy.mass_unit = str_to_float(gtk_entry_get_text
             (GTK_ENTRY(obj)));
    break;
  case CHARGE_UNIT:
    data->moldy.mass_unit = str_to_float(gtk_entry_get_text
             (GTK_ENTRY(obj)));
    break;
  }
return(FALSE);
}

/*******************************/
/* MOLDY potential type toggle */
/*******************************/
gint moldy_potential(GtkWidget *w, gpointer *obj)
{
gint pot_type;
struct model_pak *data;

pot_type = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(obj), "key"));

data = (struct model_pak *) g_object_get_data(G_OBJECT(obj), "ptr");
g_return_val_if_fail(data != NULL, FALSE);

data->moldy.pot_type = pot_type;

return(FALSE);
}
/************************/
/* MOLDY keyword toggle */
/************************/
gint moldy_keyword(GtkWidget *w, gpointer *obj)
{
gint keyword;
struct model_pak *data;

keyword = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(obj), "key"));

data = (struct model_pak *) g_object_get_data(G_OBJECT(obj), "ptr");
g_return_val_if_fail(data != NULL, FALSE);

switch(keyword)
  {
  /* Simulation parms */
  case DELTAT:
    data->moldy.deltat = str_to_float(gtk_entry_get_text
             (GTK_ENTRY(obj)))/1000;
    break;
  case STEPS:
    data->moldy.nsteps = (gint) GTK_ADJUSTMENT(obj)->value;
    break;
  case TEMPERATURE:
    data->moldy.temperature = (gdouble) GTK_ADJUSTMENT(obj)->value;
    break;
  case TITLE:
    g_free(data->moldy.title);
    data->moldy.title = g_strstrip(g_strdup(gtk_entry_get_text(GTK_ENTRY(obj))));
    break;
  case NMOLS:
    data->moldy.num_mols = (gint) GTK_ADJUSTMENT(obj)->value;
    break;
  case CONTROL_FILE:
    g_free(data->moldy.control_file);
    data->moldy.control_file = g_strstrip(g_strdup(gtk_entry_get_text(GTK_ENTRY(obj))));
    break;
  case SYSSPEC_FILE:
    g_free(data->moldy.sysspec_file);
    data->moldy.sysspec_file = g_strstrip(g_strdup(gtk_entry_get_text(GTK_ENTRY(obj))));
    break;
  case SAVE_FILE:
    g_free(data->moldy.save_file);
    data->moldy.save_file = g_strstrip(g_strdup(gtk_entry_get_text(GTK_ENTRY(obj))));
    break;
  case RESTART_FILE:
    g_free(data->moldy.restart_file);
    data->moldy.restart_file = g_strstrip(g_strdup(gtk_entry_get_text(GTK_ENTRY(obj))));
    break;
  case DUMP_FILE:
    g_free(data->moldy.dump_file);
    data->moldy.dump_file = g_strstrip(g_strdup(gtk_entry_get_text(GTK_ENTRY(obj))));
    break;
  case OUTPUT:
    g_free(data->moldy.out_file);
    data->moldy.out_file = g_strstrip(g_strdup(gtk_entry_get_text(GTK_ENTRY(obj))));
    break;
  case LIB_FILE:
    g_free(data->moldy.lib_file);
    data->moldy.lib_file = g_strstrip(g_strdup(gtk_entry_get_text(GTK_ENTRY(obj))));
    break;
  case DENSITY:
    data->moldy.density = str_to_float(gtk_entry_get_text(GTK_ENTRY(obj)));
    break;
  case PRESSURE:
    data->moldy.pressure = str_to_float(gtk_entry_get_text
             (GTK_ENTRY(obj)));
    break;
  case TTMASS:
    data->moldy.ttmass = str_to_float(gtk_entry_get_text
             (GTK_ENTRY(obj)));
    break;
  case RTMASS:
    data->moldy.rtmass = str_to_float(gtk_entry_get_text
             (GTK_ENTRY(obj)));
    break;
  case MASSPARM:
    data->moldy.mass_parm = str_to_float(gtk_entry_get_text
             (GTK_ENTRY(obj)));
    break;
  case SUBCELL: 
    data->moldy.subcell = str_to_float(gtk_entry_get_text
             (GTK_ENTRY(obj)));
    break;
  case ACCURACY: 
    data->moldy.ewald_accuracy = str_to_float(gtk_entry_get_text
             (GTK_ENTRY(obj)));
    break;
  case DUMPSTART:
    data->moldy.dump_begin = (gint) GTK_ADJUSTMENT(obj)->value;
    break;
  case DUMPINT:
    data->moldy.dump_int = (gint) GTK_ADJUSTMENT(obj)->value;
    break;
  case NDUMPS:
    data->moldy.maxdumps = (gint) GTK_ADJUSTMENT(obj)->value;
    break;
  case SCALEINT:
    data->moldy.scale_int = (gint) GTK_ADJUSTMENT(obj)->value;
    break;
  case SCALEEND:
    data->moldy.scale_end = (gint) GTK_ADJUSTMENT(obj)->value;
    break;
  case PRINTINT:
    data->moldy.print_int = (gint) GTK_ADJUSTMENT(obj)->value;
    break;
  case AVSTART:
    data->moldy.av_start = (gint) GTK_ADJUSTMENT(obj)->value;
    break;
  case AVINT:
    data->moldy.av_int = (gint) GTK_ADJUSTMENT(obj)->value;
    break;
  case ROLLINT:
    data->moldy.roll_int = (gint) GTK_ADJUSTMENT(obj)->value;
    break;
  case BACKINT:
    data->moldy.backup_int = (gint) GTK_ADJUSTMENT(obj)->value;
    break;
/*  case BACKFILE:
    g_free(data->moldy.backup_file);
    data->moldy.backup_file = g_strstrip(g_strdup(gtk_entry_get_text(GTK_ENTRY(obj))));
    break; */
  case CUTOFF:
    data->moldy.cutoff = str_to_float(gtk_entry_get_text
             (GTK_ENTRY(obj)));
    break;
  case ALPHA:
    data->moldy.alpha = str_to_float(gtk_entry_get_text
             (GTK_ENTRY(obj)));
    break;
  case KCUTOFF:
    data->moldy.kcutoff = str_to_float(gtk_entry_get_text
             (GTK_ENTRY(obj)));
    break;
  case RDFSTART:
    data->moldy.rdf_begin = (gint) GTK_ADJUSTMENT(obj)->value;
    break;
  case RDFINT:
    data->moldy.rdf_int = (gint) GTK_ADJUSTMENT(obj)->value;
    break;
  case RDFOUT:
    data->moldy.rdf_out = (gint) GTK_ADJUSTMENT(obj)->value;
    break;
  case RDFLIM:
    data->moldy.rdf_limit = str_to_float(gtk_entry_get_text
             (GTK_ENTRY(obj)));
    break;
  case NBINS:
    data->moldy.nbins = (gint) GTK_ADJUSTMENT(obj)->value;
    break;
  }

return(FALSE);
}

void moldy_unit_dialog(struct model_pak *data)
{
gpointer dialog;
GtkWidget *hbox, *vbox;
GtkWidget *frame, *label, *entry;
GtkWidget *window;

if (!data)
  return;

/* request a moldy unit dialog */
dialog = dialog_request(MOLDY_ENERGY, "Moldy units", NULL, NULL, data);
if (!dialog)
  return;
window = dialog_window(dialog);
gtk_window_set_default_size(GTK_WINDOW(window), 300, 140);

/* Frame */
frame = gtk_frame_new (NULL);
gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox), frame, FALSE, FALSE, 0); 
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);

/* create an hbox in the frame */
/* Time unit */
hbox = gtk_hbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(frame), hbox);
gtk_container_set_border_width(GTK_CONTAINER(hbox), PANEL_SPACING);

/* left vbox */
vbox = gtk_vbox_new(TRUE, 0);
gtk_box_pack_start (GTK_BOX (hbox), vbox, TRUE, TRUE, 0);

label = gtk_label_new("Time (s)");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);

label = gtk_label_new("Length (m)");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);

label = gtk_label_new("Mass (g)");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);

label = gtk_label_new("Charge (C)");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);

/* right vbox */
vbox = gtk_vbox_new(TRUE, 0);
gtk_box_pack_start (GTK_BOX (hbox), vbox, TRUE, TRUE, 0);

entry = gtk_entry_new_with_max_length(WIDTH);
gtk_entry_set_text(GTK_ENTRY(entry),
                   g_strdup_printf("%-12.12g", data->moldy.time_unit));
gtk_box_pack_start (GTK_BOX (vbox), entry, TRUE, TRUE, 0);
/* event */
g_signal_connect(GTK_OBJECT (entry), "changed",
                          GTK_SIGNAL_FUNC (unit_keyword),
                         (gpointer) entry);
g_object_set_data(G_OBJECT(entry), "key", (gpointer) TIME_UNIT);
g_object_set_data(G_OBJECT(entry), "ptr", (gpointer) data);

entry = gtk_entry_new_with_max_length(WIDTH);
gtk_entry_set_text(GTK_ENTRY(entry),
                   g_strdup_printf("%-12.12g", data->moldy.length_unit));
gtk_box_pack_start (GTK_BOX (vbox), entry, TRUE, TRUE, 0);
/* event */
g_signal_connect(GTK_OBJECT (entry), "changed",
                          GTK_SIGNAL_FUNC (unit_keyword),
                         (gpointer) entry);
g_object_set_data(G_OBJECT(entry), "key", (gpointer) LENGTH_UNIT);
g_object_set_data(G_OBJECT(entry), "ptr", (gpointer) data);

entry = gtk_entry_new_with_max_length(WIDTH);
gtk_entry_set_text(GTK_ENTRY(entry),
                   g_strdup_printf("%-12.12g", data->moldy.mass_unit));
gtk_box_pack_start (GTK_BOX (vbox), entry, TRUE, TRUE, 0);
/* event */
g_signal_connect(GTK_OBJECT (entry), "changed",
                          GTK_SIGNAL_FUNC (unit_keyword),
                         (gpointer) entry);
g_object_set_data(G_OBJECT(entry), "key", (gpointer) MASS_UNIT);
g_object_set_data(G_OBJECT(entry), "ptr", (gpointer) data);

entry = gtk_entry_new_with_max_length(WIDTH);
gtk_entry_set_text(GTK_ENTRY(entry),
                   g_strdup_printf("%-12.12g", data->moldy.charge_unit));
gtk_box_pack_start (GTK_BOX (vbox), entry, TRUE, TRUE, 0);
/* event */
g_signal_connect(GTK_OBJECT (entry), "changed",
                          GTK_SIGNAL_FUNC (unit_keyword),
                         (gpointer) entry);
g_object_set_data(G_OBJECT(entry), "key", (gpointer) CHARGE_UNIT);
g_object_set_data(G_OBJECT(entry), "ptr", (gpointer) data);

/* terminating buttons */
gui_stock_button(GTK_STOCK_CLOSE, dialog_destroy, dialog,
                   GTK_DIALOG(window)->action_area);

/* done */
gtk_widget_show_all(window);
}

/*****************************/
/* MOLDY energy unit toggles */
/*****************************/
void reset_units(struct model_pak *data)
{
dialog_destroy_type(MOLDY_ENERGY);
data->moldy.length_unit = ANGST;
data->moldy.mass_unit = AMU;
data->moldy.charge_unit = ELCHARGE;
}
void moldy_kjmol(struct model_pak *data)
{
reset_units(data);
data->moldy.energy_unit = KJMOL;
data->moldy.time_unit = TIME_KJ;
}
void moldy_kcal(struct model_pak *data)
{
reset_units(data);
data->moldy.energy_unit = KCAL;
data->moldy.time_unit = TIME_KC;
}
void moldy_ev(struct model_pak *data)
{
reset_units(data);
data->moldy.energy_unit = EV;
data->moldy.time_unit = TIME_EV;
}
void moldy_e2a(struct model_pak *data)
{
reset_units(data);
data->moldy.energy_unit = E2A;
data->moldy.time_unit = TIME_E2A;
}

/********************************/
/* MOLDY job setup/results page */
/********************************/
void gui_moldy_widget(GtkWidget *box, gpointer dialog)
{
gint i, n_cols;
GtkWidget *hbox, *hbox1, *vbox, *vbox1, *vbox2, *page;
GtkWidget *frame, *button, *label, *entry, *notebook;
GtkWidget *table, *scr_win;
GtkAdjustment *adj;
GSList *radio_group[3];
GString *line;
struct model_pak *model;
/* Sensitive frames */
GtkWidget *cutoff_box, *dump_frame;

model = dialog_model(dialog);
g_assert(model != NULL);

if (model->periodic == 2)
  return;

/* string manipulation scratchpad */
line = g_string_new(NULL);

/* create notebook */
notebook = gtk_notebook_new();
gtk_notebook_set_tab_pos(GTK_NOTEBOOK(notebook), GTK_POS_TOP);
gtk_container_add(GTK_CONTAINER(box), notebook);
gtk_notebook_set_show_border(GTK_NOTEBOOK(notebook), FALSE);

/* system page */
page = gtk_vbox_new(FALSE,0);
label = gtk_label_new(" System ");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);

/* Simulation conditions */
frame = gtk_frame_new(" Conditions ");
gtk_box_pack_start(GTK_BOX(page), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);

hbox = gtk_hbox_new(FALSE, 0);
gtk_container_add (GTK_CONTAINER(frame),hbox);
gtk_container_set_border_width(GTK_CONTAINER(hbox), PANEL_SPACING);

/* left box */
vbox = gtk_vbox_new(TRUE, 2);
gtk_box_pack_start (GTK_BOX (hbox), vbox, TRUE, TRUE, 0);

label = gtk_label_new("Temperature (K)");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);
label = gtk_label_new("Pressure (MPa)");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);

/* right box */
vbox = gtk_vbox_new(TRUE, 2);
gtk_box_pack_start (GTK_BOX (hbox), vbox, TRUE, TRUE, 0);

adj = (GtkAdjustment *) gtk_adjustment_new
                        (model->moldy.temperature, 0, 10000, 1, 1, 0);
button = gtk_spin_button_new(adj, 0, 0);
gtk_spin_button_set_wrap(GTK_SPIN_BUTTON(button), FALSE);
gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT (adj), "value_changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) adj);
g_object_set_data(G_OBJECT(adj), "key", (gpointer) TEMPERATURE);
g_object_set_data(G_OBJECT(adj), "ptr", (gpointer) model);
gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(button), TRUE);

entry = gtk_entry_new_with_max_length(WIDTH);
gtk_entry_set_text(GTK_ENTRY(entry),
                   g_strdup_printf("%g", model->moldy.pressure));
gtk_box_pack_start (GTK_BOX (vbox), entry, TRUE, TRUE, 0);

/* event */
g_signal_connect(GTK_OBJECT (entry), "changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) entry);
g_object_set_data(G_OBJECT(entry), "key", (gpointer) PRESSURE);
g_object_set_data(G_OBJECT(entry), "ptr", (gpointer) model);

if (model->periodic == 0 )
  {
  frame = gtk_frame_new(" Skew start ");
  gtk_box_pack_start(GTK_BOX(page), frame, FALSE, FALSE, 0);
  gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
  hbox = gtk_hbox_new(FALSE, 0);
  gtk_container_add(GTK_CONTAINER(frame),hbox);
  gtk_container_set_border_width(GTK_CONTAINER(hbox), PANEL_SPACING);

  /* left vbox */
  vbox = gtk_vbox_new(TRUE, 2);
  gtk_box_pack_start (GTK_BOX (hbox), vbox, TRUE, TRUE, 0);

  label = gtk_label_new("No of molecules");
  gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);
  label = gtk_label_new("Density (g/cc)");
  gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);

  /* right vbox */
  vbox = gtk_vbox_new(TRUE, 2);
  gtk_box_pack_start (GTK_BOX (hbox), vbox, TRUE, TRUE, 0);

  adj = (GtkAdjustment *) gtk_adjustment_new
                          (model->moldy.num_mols, 1, 1000000, 1, 1, 0);
  button = gtk_spin_button_new(adj, 0, 0);
  gtk_spin_button_set_wrap(GTK_SPIN_BUTTON(button), FALSE);
  gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
  gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(button), TRUE);
  g_signal_connect(GTK_OBJECT (adj), "value_changed",
                            GTK_SIGNAL_FUNC (moldy_keyword),
                           (gpointer) adj);
  g_object_set_data(G_OBJECT(adj), "key", (gpointer) NMOLS);
  g_object_set_data(G_OBJECT(adj), "ptr", (gpointer) model);

  entry = gtk_entry_new_with_max_length(WIDTH);
  gtk_entry_set_text(GTK_ENTRY(entry),
                     g_strdup_printf("%g", model->moldy.density));
  gtk_box_pack_start (GTK_BOX (vbox), entry, TRUE, TRUE, 0);
  /* event */
  g_signal_connect(GTK_OBJECT (entry), "changed",
                            GTK_SIGNAL_FUNC (moldy_keyword),
                           (gpointer) entry);
  g_object_set_data(G_OBJECT(entry), "key", (gpointer) DENSITY);
  g_object_set_data(G_OBJECT(entry), "ptr", (gpointer) model);
  }

/* set up the time specs frame */
frame = gtk_frame_new(" Time ");
gtk_box_pack_start(GTK_BOX(page), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
hbox = gtk_hbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(frame), hbox);
gtk_container_set_border_width(GTK_CONTAINER(hbox), PANEL_SPACING);

/* left vbox */
vbox = gtk_vbox_new(TRUE, 0);
gtk_box_pack_start (GTK_BOX (hbox), vbox, TRUE, TRUE, 0);

label = gtk_label_new("Step size (fs)");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);
label = gtk_label_new("Time steps");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);

/* right vbox */
vbox = gtk_vbox_new(TRUE, 0);
gtk_box_pack_start (GTK_BOX (hbox), vbox, TRUE, TRUE, 0);

entry = gtk_entry_new_with_max_length(WIDTH);
gtk_entry_set_text(GTK_ENTRY(entry),
                   g_strdup_printf("%g", model->moldy.deltat*1000));
gtk_box_pack_start (GTK_BOX (vbox), entry, TRUE, TRUE, 0);
/* event */
g_signal_connect(GTK_OBJECT (entry), "changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) entry);
g_object_set_data(G_OBJECT(entry), "key", (gpointer) DELTAT);
g_object_set_data(G_OBJECT(entry), "ptr", (gpointer) model);

adj = (GtkAdjustment *) gtk_adjustment_new
                        (model->moldy.nsteps, 0, 1000000, 1, 1, 0);
button = gtk_spin_button_new(adj, 0, 0);
gtk_spin_button_set_wrap(GTK_SPIN_BUTTON(button), FALSE);
gtk_box_pack_start(GTK_BOX(vbox), button, TRUE, TRUE, 0);
gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT (adj), "value_changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) adj);
g_object_set_data(G_OBJECT(adj), "key", (gpointer) STEPS);
g_object_set_data(G_OBJECT(adj), "ptr", (gpointer) model);

frame = gtk_frame_new(" Energy units ");
gtk_box_pack_start(GTK_BOX(page), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
hbox = gtk_hbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(frame), hbox);

/* do the radio buttons */
new_radio_group(0, hbox, TF);
button = add_radio_button("kJ/mol", (gpointer) moldy_kjmol, model);
if (model->moldy.energy_unit == KJMOL)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

button = add_radio_button("kcal/mol", (gpointer) moldy_kcal, model);
if (model->moldy.energy_unit == KCAL)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

button = add_radio_button("eV", (gpointer) moldy_ev, model);
if (model->moldy.energy_unit == EV)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

button = add_radio_button("e^2/A", (gpointer) moldy_e2a, model);
if (model->moldy.energy_unit == E2A)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

button = add_radio_button("Other", (gpointer) moldy_unit_dialog, model);
if (model->moldy.energy_unit == OTHERS)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* Moldy files */
page = gtk_vbox_new(FALSE,0);
label = gtk_label_new(" Files ");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);

frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(page), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);

hbox = gtk_hbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(frame),hbox);
gtk_container_set_border_width(GTK_CONTAINER(hbox), PANEL_SPACING);

/* set up the filename entry */
/* left vbox */
vbox = gtk_vbox_new(TRUE, 2);
gtk_box_pack_start (GTK_BOX (hbox), vbox, TRUE, TRUE, 0);

/* temporary MOLDY input filename */
label = gtk_label_new(" Control file ");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);
label = gtk_label_new(" Sys-spec file ");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);
label = gtk_label_new(" Save file ");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);
label = gtk_label_new(" Restart file ");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);
label = gtk_label_new(" Dump file ");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);
label = gtk_label_new(" Output file ");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);
/* label = gtk_label_new("Backup file");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0); */

/* right vbox */
vbox = gtk_vbox_new(TRUE, 2);
gtk_box_pack_start (GTK_BOX (hbox), vbox, TRUE, TRUE, 0);

/* control file */
entry = gtk_entry_new_with_max_length(WIDTH);
gtk_entry_set_text(GTK_ENTRY(entry), model->moldy.control_file);
gtk_box_pack_start (GTK_BOX (vbox), entry, TRUE, TRUE, 0);
/* event */
g_signal_connect(GTK_OBJECT (entry), "changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) entry);
g_object_set_data(G_OBJECT(entry), "key", (gpointer) CONTROL_FILE);
g_object_set_data(G_OBJECT(entry), "ptr", (gpointer) model);

/* sys-spec file */
entry = gtk_entry_new_with_max_length(WIDTH);
gtk_entry_set_text(GTK_ENTRY(entry), model->moldy.sysspec_file);
gtk_box_pack_start (GTK_BOX (vbox), entry, TRUE, TRUE, 0);
/* event */
g_signal_connect(GTK_OBJECT (entry), "changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) entry);
g_object_set_data(G_OBJECT(entry), "key", (gpointer) SYSSPEC_FILE);
g_object_set_data(G_OBJECT(entry), "ptr", (gpointer) model);

/* save file */
entry = gtk_entry_new_with_max_length(WIDTH);
gtk_entry_set_text(GTK_ENTRY(entry), model->moldy.save_file);
gtk_box_pack_start (GTK_BOX (vbox), entry, TRUE, TRUE, 0);
/* event */
g_signal_connect(GTK_OBJECT (entry), "changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) entry);
g_object_set_data(G_OBJECT(entry), "key", (gpointer) SAVE_FILE);
g_object_set_data(G_OBJECT(entry), "ptr", (gpointer) model);

/* restart file */
entry = gtk_entry_new_with_max_length(WIDTH);
gtk_entry_set_text(GTK_ENTRY(entry), model->moldy.restart_file);
gtk_box_pack_start (GTK_BOX (vbox), entry, TRUE, TRUE, 0);
/* event */
g_signal_connect(GTK_OBJECT (entry), "changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) entry);
g_object_set_data(G_OBJECT(entry), "key", (gpointer) RESTART_FILE);
g_object_set_data(G_OBJECT(entry), "ptr", (gpointer) model);

/* dump file */
entry = gtk_entry_new_with_max_length(WIDTH);
gtk_entry_set_text(GTK_ENTRY(entry), model->moldy.dump_file);
gtk_box_pack_start (GTK_BOX (vbox), entry, TRUE, TRUE, 0);
/* event */
g_signal_connect(GTK_OBJECT (entry), "changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) entry);
g_object_set_data(G_OBJECT(entry), "key", (gpointer) DUMP_FILE);
g_object_set_data(G_OBJECT(entry), "ptr", (gpointer) model);

/* output file */
entry = gtk_entry_new_with_max_length(WIDTH);
gtk_entry_set_text(GTK_ENTRY(entry), model->moldy.out_file);
gtk_box_pack_start (GTK_BOX (vbox), entry, TRUE, TRUE, 0);
/* event */
g_signal_connect(GTK_OBJECT (entry), "changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) entry);
g_object_set_data(G_OBJECT(entry), "key", (gpointer) OUTPUT);
g_object_set_data(G_OBJECT(entry), "ptr", (gpointer) model);

/* input file options */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(page), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(frame), hbox);
gtk_container_set_border_width(GTK_CONTAINER(hbox), PANEL_SPACING);

/* write input files button */
gui_direct_check("Create input files then stop", &model->moldy.no_exec, NULL, NULL, hbox);

/* Load restart file button */
/*
button = gtk_button_new_with_label(" Load restart file ");
gtk_box_pack_start(GTK_BOX (hbox), button, TRUE, FALSE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                    GTK_SIGNAL_FUNC(gtk_load_moldy), (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) 1);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);
*/

/* starting new page */
page = gtk_vbox_new(FALSE,0);
label = gtk_label_new(" Ensemble ");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);

/* Ensemble settings */
hbox = gtk_hbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(page),hbox);
vbox1 = gtk_vbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(hbox),vbox1);
vbox2 = gtk_vbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(hbox),vbox2);

/* control parameters */
frame = gtk_frame_new(" Temperature control ");
gtk_box_pack_start(GTK_BOX(vbox1), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);

/* do the first (DEFAULT MODE) radio button */
button = gtk_radio_button_new_with_label (NULL, "Scaling");
gtk_box_pack_start (GTK_BOX (vbox), button, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (thermo_keyword),
                         (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) SCALING);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);
if (model->moldy.const_temp == SCALING || model->moldy.const_temp == -1)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* make a radio group */
radio_group[1] = gtk_radio_button_group (GTK_RADIO_BUTTON (button));
/* do the next button */
button = gtk_radio_button_new_with_label (radio_group[1], "Nose-Hoover");
gtk_box_pack_start (GTK_BOX (vbox), button, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (thermo_keyword),
                         (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) NOSEHOOVER);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);
if (model->moldy.const_temp == NOSEHOOVER)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* subsequent button... */
button = gtk_radio_button_new_with_label(gtk_radio_button_group
       (GTK_RADIO_BUTTON (button)), "Gaussian");
gtk_box_pack_start (GTK_BOX (vbox), button, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (thermo_keyword),
                         (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) GAUSSIAN);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);
if (model->moldy.const_temp == GAUSSIAN)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* framework option */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(page), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(frame), hbox);
gtk_container_set_border_width(GTK_CONTAINER(hbox), PANEL_SPACING);

/* set framework option */
gui_direct_check("Set as framework", &model->moldy.set_framework, NULL, NULL, hbox);

/* pressure management */
frame = gtk_frame_new(" Pressure control ");
gtk_box_pack_start(GTK_BOX(vbox2), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);

/* do the first (DEFAULT MODE) radio button */
button = gtk_radio_button_new_with_label (NULL, "Constant volume");
gtk_box_pack_start (GTK_BOX (vbox), button, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (pressure_keyword),
                         (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) NOPRESS);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);
if (model->moldy.const_press == NOPRESS || model->moldy.const_press == -1)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* make a radio group */
radio_group[2] = gtk_radio_button_group (GTK_RADIO_BUTTON (button));
/* do the next button */
button = gtk_radio_button_new_with_label (radio_group[2], "Andersen constant pressure");
gtk_box_pack_start (GTK_BOX (vbox), button, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (pressure_keyword),
                         (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) ANDERSEN);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);
if (model->moldy.const_press == ANDERSEN)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* do the next button */
button = gtk_radio_button_new_with_label(gtk_radio_button_group
       (GTK_RADIO_BUTTON (button)), "P-R constant pressure");
gtk_box_pack_start (GTK_BOX (vbox), button, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (pressure_keyword),
                         (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) PRPRESS);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);
if (model->moldy.const_press == PRPRESS)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* do the next button */
button = gtk_radio_button_new_with_label (gtk_radio_button_group
       (GTK_RADIO_BUTTON (button)), "P-R constant stress");
gtk_box_pack_start (GTK_BOX (vbox), button, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (pressure_keyword),
                         (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) PRSTRESS);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);
if (model->moldy.const_press == PRSTRESS)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* do the next button */
button = gtk_radio_button_new_with_label(gtk_radio_button_group
       (GTK_RADIO_BUTTON (button)), "C-W constant stress");
gtk_box_pack_start (GTK_BOX (vbox), button, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (pressure_keyword),
                         (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) WCSTRESS);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);
if (model->moldy.const_press == WCSTRESS)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);


/* starting new *page */
page = gtk_vbox_new(FALSE,0);
label = gtk_label_new(" Options ");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);

/* Next section */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(page), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
hbox1 = gtk_hbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(frame), hbox1);

/* left vbox */
vbox1 = gtk_vbox_new(TRUE, 2);
gtk_box_pack_start (GTK_BOX (hbox1), vbox1, TRUE, TRUE, 0);

label = gtk_label_new("Scaling interval");
gtk_box_pack_start (GTK_BOX (vbox1), label, TRUE, TRUE, 0);
label = gtk_label_new("Scaling end");
gtk_box_pack_start (GTK_BOX (vbox1), label, TRUE, TRUE, 0);
label = gtk_label_new("Translational inertia");
gtk_box_pack_start (GTK_BOX (vbox1), label, TRUE, TRUE, 0);
label = gtk_label_new("Rotational inertia");
gtk_box_pack_start (GTK_BOX (vbox1), label, TRUE, TRUE, 0);
label = gtk_label_new("Mass parameter");
gtk_box_pack_start (GTK_BOX (vbox1), label, TRUE, TRUE, 0);

/* right vbox */
vbox2 = gtk_vbox_new(TRUE, 2);
gtk_box_pack_start (GTK_BOX (hbox1), vbox2, TRUE, TRUE, 0);

adj = (GtkAdjustment *) gtk_adjustment_new
                        (model->moldy.scale_int, 0, 1e6, 1, 1, 0);
button = gtk_spin_button_new(adj, 0, 0);
gtk_spin_button_set_wrap(GTK_SPIN_BUTTON(button), FALSE);
gtk_box_pack_start(GTK_BOX(vbox2), button, FALSE, FALSE, 0);
gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT (adj), "value_changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) adj);
g_object_set_data(G_OBJECT(adj), "key", (gpointer) SCALEINT);
g_object_set_data(G_OBJECT(adj), "ptr", (gpointer) model);

adj = (GtkAdjustment *) gtk_adjustment_new
                        (model->moldy.scale_end, 0, 1e6, 1, 1, 0);
button = gtk_spin_button_new(adj, 0, 0);
gtk_spin_button_set_wrap(GTK_SPIN_BUTTON(button), FALSE);
gtk_box_pack_start(GTK_BOX(vbox2), button, FALSE, FALSE, 0);
gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT (adj), "value_changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) adj);
g_object_set_data(G_OBJECT(adj), "key", (gpointer) SCALEEND);
g_object_set_data(G_OBJECT(adj), "ptr", (gpointer) model);

entry = gtk_entry_new_with_max_length(WIDTH);
gtk_entry_set_text(GTK_ENTRY(entry),
                   g_strdup_printf("%g", model->moldy.ttmass));
gtk_box_pack_start (GTK_BOX (vbox2), entry, TRUE, TRUE, 0);
/* event */
g_signal_connect(GTK_OBJECT (entry), "changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) entry);
g_object_set_data(G_OBJECT(entry), "key", (gpointer) TTMASS);
g_object_set_data(G_OBJECT(entry), "ptr", (gpointer) model);

entry = gtk_entry_new_with_max_length(WIDTH);
gtk_entry_set_text(GTK_ENTRY(entry),
                   g_strdup_printf("%g", model->moldy.rtmass));
gtk_box_pack_start (GTK_BOX (vbox2), entry, TRUE, TRUE, 0);
/* event */
g_signal_connect(GTK_OBJECT (entry), "changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) entry);
g_object_set_data(G_OBJECT(entry), "key", (gpointer) RTMASS);
g_object_set_data(G_OBJECT(entry), "ptr", (gpointer) model);

entry = gtk_entry_new_with_max_length(WIDTH);
gtk_entry_set_text(GTK_ENTRY(entry),
                   g_strdup_printf("%g", model->moldy.mass_parm));
gtk_box_pack_start (GTK_BOX (vbox2), entry, TRUE, TRUE, 0);
/* event */
g_signal_connect(GTK_OBJECT (entry), "changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) entry);
g_object_set_data(G_OBJECT(entry), "key", (gpointer) MASSPARM);
g_object_set_data(G_OBJECT(entry), "ptr", (gpointer) model);

/* split panel */
hbox = gtk_hbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(page),hbox);
vbox1 = gtk_vbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(hbox),vbox1);
vbox2 = gtk_vbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(hbox),vbox2);

/* Create strain mask table */
frame = gtk_frame_new(" Strain mask ");
gtk_box_pack_start(GTK_BOX(vbox1), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
hbox = gtk_hbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(frame), hbox);

table = gtk_table_new (5,4,FALSE);
gtk_box_pack_start(GTK_BOX(hbox), table, TRUE, FALSE, 0);
label = gtk_label_new("a");
gtk_table_attach_defaults (GTK_TABLE(table), label, 1, 2, 0, 1);
label = gtk_label_new("b");
gtk_table_attach_defaults (GTK_TABLE(table), label, 2, 3, 0, 1);
label = gtk_label_new("c");
gtk_table_attach_defaults (GTK_TABLE(table), label, 3, 4, 0, 1);
label = gtk_label_new("x ");
gtk_table_attach_defaults (GTK_TABLE(table), label, 0, 1, 1, 2);
label = gtk_label_new("y ");
gtk_table_attach_defaults (GTK_TABLE(table), label, 0, 1, 2, 3);
label = gtk_label_new("z ");
gtk_table_attach_defaults (GTK_TABLE(table), label, 0, 1, 3, 4);
label = gtk_label_new(" ");
gtk_table_attach_defaults (GTK_TABLE(table), label, 0, 1, 4, 5);

button = gtk_check_button_new();
gtk_table_attach_defaults (GTK_TABLE(table), button, 1, 2, 1, 2);
if( model->moldy.strain_mask & 1)
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (strain_mask_toggle),
                         (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) 1);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);

button = gtk_check_button_new();
gtk_table_attach_defaults (GTK_TABLE(table), button, 2, 3, 1, 2);
if( model->moldy.strain_mask & 2)
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (strain_mask_toggle),
                         (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) 2);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);

button = gtk_check_button_new();
gtk_table_attach_defaults (GTK_TABLE(table), button, 3, 4, 1, 2);
if( model->moldy.strain_mask & 4)
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (strain_mask_toggle),
                         (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) 4);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);

button = gtk_check_button_new();
gtk_table_attach_defaults (GTK_TABLE(table), button, 1, 2, 2, 3);
if( model->moldy.strain_mask & 8)
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (strain_mask_toggle),
                         (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) 8);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);

button = gtk_check_button_new();
gtk_table_attach_defaults (GTK_TABLE(table), button, 2, 3, 2, 3);
if( model->moldy.strain_mask & 16)
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (strain_mask_toggle),
                         (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) 16);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);

button = gtk_check_button_new();
gtk_table_attach_defaults (GTK_TABLE(table), button, 3, 4, 2, 3);
if( model->moldy.strain_mask & 32)
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (strain_mask_toggle),
                         (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) 32);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);

button = gtk_check_button_new();
gtk_table_attach_defaults (GTK_TABLE(table), button, 1, 2, 3, 4);
if( model->moldy.strain_mask & 64)
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (strain_mask_toggle),
                         (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) 64);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);

button = gtk_check_button_new();
gtk_table_attach_defaults (GTK_TABLE(table), button, 2, 3, 3, 4);
if( model->moldy.strain_mask & 128)
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (strain_mask_toggle),
                         (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) 128);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);

button = gtk_check_button_new();
gtk_table_attach_defaults (GTK_TABLE(table), button, 3, 4, 3, 4);
if( model->moldy.strain_mask & 256)
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (strain_mask_toggle),
                         (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) 256);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);

/* Scaling options boxes */
frame = gtk_frame_new(" Scaling options ");
gtk_box_pack_start(GTK_BOX(vbox2), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);

vbox = gtk_vbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);
button = gtk_check_button_new_with_label("Individual species");
if( model->moldy.scale_options & 1)
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
gtk_box_pack_start (GTK_BOX (vbox), button, FALSE, FALSE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (scale_opt_toggle),
                         (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) 1);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);

button = gtk_check_button_new_with_label("Separate rot/trans");
if( model->moldy.scale_options & 2)
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
gtk_box_pack_start (GTK_BOX (vbox), button, FALSE, FALSE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (scale_opt_toggle),
                         (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) 2);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);

button = gtk_check_button_new_with_label("Use rolling average");
if( model->moldy.scale_options & 4)
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
gtk_box_pack_start (GTK_BOX (vbox), button, FALSE, FALSE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (scale_opt_toggle),
                         (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) 4);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);

button = gtk_check_button_new_with_label("Re-initialize");
if( model->moldy.scale_options & 8)
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
gtk_box_pack_start (GTK_BOX (vbox), button, FALSE, FALSE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (scale_opt_toggle),
                         (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) 8);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);

/* starting new page */
page = gtk_vbox_new(FALSE,0);
label = gtk_label_new(" Potentials ");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);

hbox = gtk_hbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(page), hbox);

/* split panel */
vbox1 = gtk_vbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(hbox), vbox1, FALSE, FALSE, 0);
vbox2 = gtk_vbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(hbox), vbox2, TRUE, TRUE, 0);

/* potential type */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(vbox1), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);

vbox = gtk_vbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);
label = gtk_label_new("Potential type");
gtk_box_pack_start(GTK_BOX (vbox), label, FALSE, FALSE, 5);

line->str = choose_potential(model->moldy.pot_type);
model->moldy.pot_label = gtk_entry_new_with_max_length(LINELEN);
gtk_entry_set_text(GTK_ENTRY(model->moldy.pot_label), line->str);
gtk_entry_set_editable(GTK_ENTRY(model->moldy.pot_label), FALSE);
gtk_box_pack_start (GTK_BOX (vbox), model->moldy.pot_label, FALSE, FALSE, 0);

frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(vbox1), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
 
/* Potential library */
vbox = gtk_vbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(frame), vbox);
label = gtk_label_new("Potential library");
gtk_box_pack_start (GTK_BOX (vbox), label, FALSE, FALSE, 5);

model->moldy.lib_file_entry = gtk_entry_new_with_max_length(LINELEN);
gtk_entry_set_text(GTK_ENTRY(model->moldy.lib_file_entry), model->moldy.lib_file);
gtk_entry_set_editable(GTK_ENTRY(model->moldy.lib_file_entry), FALSE);
gtk_box_pack_start (GTK_BOX (vbox), model->moldy.lib_file_entry, FALSE, FALSE, 0);

/* read potential library file button */
vbox = gtk_vbox_new(TRUE, PANEL_SPACING);
gtk_box_pack_start(GTK_BOX(vbox1), vbox, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(vbox), PANEL_SPACING);

button = gtk_button_new_with_label(" Load parameters ");
gtk_box_pack_start(GTK_BOX (vbox), button, TRUE, FALSE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                    GTK_SIGNAL_FUNC(gtk_read_potlib), (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) 0);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);

button = gtk_button_new_with_label(" Save parameters ");
gtk_box_pack_start(GTK_BOX (vbox), button, TRUE, FALSE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                    GTK_SIGNAL_FUNC(gtk_write_potlib), (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) 0);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);

button = gtk_button_new_with_label(" Clear parameters ");
gtk_box_pack_start(GTK_BOX (vbox), button, TRUE, FALSE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                    GTK_SIGNAL_FUNC(gtk_delete_pots), (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) 0);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);

/* Potential parameter list */
frame = gtk_frame_new(" Current parameter list ");
gtk_box_pack_start(GTK_BOX(vbox2), frame, FALSE, FALSE, 0);
 
vbox = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(frame), vbox);
gtk_container_set_border_width(GTK_CONTAINER(GTK_BOX(vbox)), PANEL_SPACING);

/* scrolled model pane */
scr_win = gtk_scrolled_window_new (NULL, NULL);
gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (scr_win),
                                GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
gtk_box_pack_start(GTK_BOX(vbox), scr_win, TRUE, TRUE, 0);
 
gtk_widget_set_usize(scr_win, 240, 240);

/* create the parm list */
n_cols = 4+MAX_POT_PARMS;
pot_list = gtk_clist_new_with_titles(n_cols, pot_titles);
gtk_clist_set_selection_mode(GTK_CLIST(pot_list), GTK_SELECTION_SINGLE);
/* set some properties of the list */
/* FIXME - this should probably be dependent on font size */
gtk_clist_set_column_width(GTK_CLIST(pot_list), 0, 40);
gtk_clist_set_column_width(GTK_CLIST(pot_list), 1, 40);
gtk_clist_set_column_width(GTK_CLIST(pot_list), 2, 40);
gtk_clist_set_column_width(GTK_CLIST(pot_list), 3, 40);
for(i=0;i<n_cols;i++)
  gtk_clist_set_column_width(GTK_CLIST(pot_list), 4+i, 60);

gtk_clist_set_shadow_type (GTK_CLIST(pot_list), GTK_SHADOW_OUT);
gtk_container_add(GTK_CONTAINER(scr_win), pot_list);

gtk_clist_column_titles_passive(GTK_CLIST(pot_list));
for(i=0;i<n_cols;i++)
   gtk_clist_set_column_justification(GTK_CLIST(pot_list), i, GTK_JUSTIFY_CENTER);

update_pot_dialog(model);

/* starting new page */
page = gtk_vbox_new(FALSE,0);
label = gtk_label_new(" Cutoffs ");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);

/* Cutoff settings */
frame = gtk_frame_new(" Method ");
gtk_box_pack_start(GTK_BOX(page), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
hbox = gtk_hbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(frame), hbox);

/* left vbox */
vbox1 = gtk_vbox_new(TRUE, 2);
gtk_box_pack_start (GTK_BOX (hbox), vbox1, TRUE, FALSE, 0);

/* Calculate cutoffs automatically */
button = gtk_check_button_new_with_label("Auto cutoffs");
gtk_box_pack_start(GTK_BOX(vbox1), button, TRUE, FALSE, 0);
if (model->moldy.auto_cutoff)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (auto_toggle), (gpointer) model);

/* strict cutoff toggle */
button = gtk_check_button_new_with_label("Strict cutoff");
gtk_box_pack_start(GTK_BOX(vbox1), button, TRUE, FALSE, 0);
if (model->moldy.strict)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (strict_toggle), (gpointer) model);

/* surface dipole */
button = gtk_check_button_new_with_label("Surface dipole");
gtk_box_pack_start(GTK_BOX(vbox1), button, TRUE, FALSE, 0);
if (model->moldy.surf)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (surf_toggle), (gpointer) model);

/* right vbox */
vbox2 = gtk_vbox_new(TRUE, 2);
gtk_box_pack_start (GTK_BOX (hbox), vbox2, FALSE, FALSE, 0);

/* real-space only */
button = gtk_check_button_new_with_label("Real-space only");
gtk_box_pack_start(GTK_BOX(vbox2), button, TRUE, FALSE, 0);
if (model->moldy.real_only)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (real_toggle), (gpointer) model);

/* no coulomb */
button = gtk_check_button_new_with_label("Include coulomb forces");
gtk_box_pack_start(GTK_BOX(vbox2), button, TRUE, FALSE, 0);
if (model->moldy.ewald)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (ewald_toggle), (gpointer) model);

/* Cutoff parameters */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(page), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
cutoff_box = gtk_hbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(frame), cutoff_box);

/* model->moldy.cutoff_frame = cutoff_box; */

/* left vbox */
vbox = gtk_vbox_new(TRUE, 2);
gtk_box_pack_start (GTK_BOX (cutoff_box), vbox, TRUE, TRUE, 0);
/* gtk_box_pack_start (cutoff_box, vbox, TRUE, TRUE, 0);*/

label = gtk_label_new("Real-space cutoff (A)");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);
label = gtk_label_new("Alpha parameter");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);
label = gtk_label_new("Reciprocal-space cutoff (1/A)");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);

/* right vbox */
vbox = gtk_vbox_new(TRUE, 2);
gtk_box_pack_start (GTK_BOX (cutoff_box), vbox, TRUE, TRUE, 0);
/* gtk_box_pack_start (cutoff_box, vbox, TRUE, TRUE, 0); */

entry = gtk_entry_new_with_max_length(WIDTH);
gtk_entry_set_text(GTK_ENTRY(entry),g_strdup_printf("%g",model->moldy.cutoff));
gtk_box_pack_start (GTK_BOX (vbox), entry, TRUE, TRUE, 0);
/* event */
g_signal_connect(GTK_OBJECT (entry), "changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) entry);
g_object_set_data(G_OBJECT(entry), "key", (gpointer) CUTOFF);
g_object_set_data(G_OBJECT(entry), "ptr", (gpointer) model);

entry = gtk_entry_new_with_max_length(WIDTH);
gtk_entry_set_text(GTK_ENTRY(entry),g_strdup_printf("%g", model->moldy.alpha));
gtk_box_pack_start (GTK_BOX (vbox), entry, TRUE, TRUE, 0);
/* event */
g_signal_connect(GTK_OBJECT (entry), "changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) entry);
g_object_set_data(G_OBJECT(entry), "key", (gpointer) ALPHA);
g_object_set_data(G_OBJECT(entry), "ptr", (gpointer) model);

entry = gtk_entry_new_with_max_length(WIDTH);
gtk_entry_set_text(GTK_ENTRY(entry),g_strdup_printf("%g", model->moldy.kcutoff));
gtk_box_pack_start (GTK_BOX (vbox), entry, TRUE, TRUE, 0);
/* event */
g_signal_connect(GTK_OBJECT (entry), "changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) entry);
g_object_set_data(G_OBJECT(entry), "key", (gpointer) KCUTOFF);
g_object_set_data(G_OBJECT(entry), "ptr", (gpointer) model);

/* Other parameters */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(page), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
hbox = gtk_hbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(frame), hbox);

/* left vbox */
vbox = gtk_vbox_new(TRUE, 2);
gtk_box_pack_start (GTK_BOX (hbox), vbox, TRUE, TRUE, 0);

label = gtk_label_new("Ewald accuracy");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);
label = gtk_label_new("Link cell side length (A)");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);

/* right vbox */
vbox = gtk_vbox_new(TRUE, 2);
gtk_box_pack_start (GTK_BOX (hbox), vbox, TRUE, TRUE, 0);

entry = gtk_entry_new_with_max_length(WIDTH);
gtk_entry_set_text(GTK_ENTRY(entry),g_strdup_printf("%g", model->moldy.ewald_accuracy));
gtk_box_pack_start (GTK_BOX (vbox), entry, TRUE, TRUE, 0);
/* event */
g_signal_connect(GTK_OBJECT (entry), "changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) entry);
g_object_set_data(G_OBJECT(entry), "key", (gpointer) ACCURACY);
g_object_set_data(G_OBJECT(entry), "ptr", (gpointer) model);

entry = gtk_entry_new_with_max_length(WIDTH);
gtk_entry_set_text(GTK_ENTRY(entry),g_strdup_printf("%g", model->moldy.subcell));
gtk_box_pack_start (GTK_BOX (vbox), entry, TRUE, TRUE, 0);
/* event */
g_signal_connect(GTK_OBJECT (entry), "changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) entry);
g_object_set_data(G_OBJECT(entry), "key", (gpointer) SUBCELL);
g_object_set_data(G_OBJECT(entry), "ptr", (gpointer) model);

/* starting new page */
page = gtk_vbox_new(FALSE,0);
label = gtk_label_new(" RDF ");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);

/* RDF parameters */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(page), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
hbox = gtk_hbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(frame), hbox);

/* left vbox */
vbox = gtk_vbox_new(TRUE, 2);
gtk_box_pack_start (GTK_BOX (hbox), vbox, TRUE, TRUE, 0);

label = gtk_label_new("RDF start");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);
label = gtk_label_new("RDF interval");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);
label = gtk_label_new("RDF out");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);
label = gtk_label_new("RDF limit (A)");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);
label = gtk_label_new("Binning intervals");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);

/* right vbox */
vbox = gtk_vbox_new(TRUE, 2);
gtk_box_pack_start (GTK_BOX (hbox), vbox, TRUE, TRUE, 0);

adj = (GtkAdjustment *) gtk_adjustment_new
                        (model->moldy.rdf_begin, 0, 1e6, 1, 1, 0);
button = gtk_spin_button_new(adj, 0, 0);
gtk_spin_button_set_wrap(GTK_SPIN_BUTTON(button), FALSE);
gtk_box_pack_start(GTK_BOX(vbox), button, FALSE, FALSE, 0);
gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT (adj), "value_changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) adj);
g_object_set_data(G_OBJECT(adj), "key", (gpointer) RDFSTART);
g_object_set_data(G_OBJECT(adj), "ptr", (gpointer) model);

adj = (GtkAdjustment *) gtk_adjustment_new
                        (model->moldy.rdf_int, 0, 1e6, 1, 1, 0);
button = gtk_spin_button_new(adj, 0, 0);
gtk_spin_button_set_wrap(GTK_SPIN_BUTTON(button), FALSE);
gtk_box_pack_start(GTK_BOX(vbox), button, FALSE, FALSE, 0);
gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT (adj), "value_changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) adj);
g_object_set_data(G_OBJECT(adj), "key", (gpointer) RDFINT);
g_object_set_data(G_OBJECT(adj), "ptr", (gpointer) model);

adj = (GtkAdjustment *) gtk_adjustment_new
                        (model->moldy.rdf_out, 0, 1e6, 1, 1, 0);
button = gtk_spin_button_new(adj, 0, 0);
gtk_spin_button_set_wrap(GTK_SPIN_BUTTON(button), FALSE);
gtk_box_pack_start(GTK_BOX(vbox), button, FALSE, FALSE, 0);
gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT (adj), "value_changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) adj);
g_object_set_data(G_OBJECT(adj), "key", (gpointer) RDFOUT);
g_object_set_data(G_OBJECT(adj), "ptr", (gpointer) model);

entry = gtk_entry_new_with_max_length(WIDTH);
gtk_entry_set_text(GTK_ENTRY(entry),g_strdup_printf("%g", model->moldy.rdf_limit));
gtk_box_pack_start (GTK_BOX (vbox), entry, TRUE, TRUE, 0);
/* event */
g_signal_connect(GTK_OBJECT (entry), "changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) entry);
g_object_set_data(G_OBJECT(entry), "key", (gpointer) RDFLIM);
g_object_set_data(G_OBJECT(entry), "ptr", (gpointer) model);

adj = (GtkAdjustment *) gtk_adjustment_new
                        (model->moldy.nbins, 0, 1e6, 1, 1, 0);
button = gtk_spin_button_new(adj, 0, 0);
gtk_spin_button_set_wrap(GTK_SPIN_BUTTON(button), FALSE);
gtk_box_pack_start(GTK_BOX(vbox), button, FALSE, FALSE, 0);
gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT (adj), "value_changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) adj);
g_object_set_data(G_OBJECT(adj), "key", (gpointer) NBINS);
g_object_set_data(G_OBJECT(adj), "ptr", (gpointer) model);



/* starting new page */
page = gtk_vbox_new(FALSE,0);
label = gtk_label_new(" Output ");
gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);

/* output parameters */
frame = gtk_frame_new(" Parameters ");
gtk_box_pack_start(GTK_BOX(page), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
hbox = gtk_hbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(frame), hbox);

/* left vbox */
vbox = gtk_vbox_new(TRUE, 2);
gtk_box_pack_start (GTK_BOX (hbox), vbox, TRUE, TRUE, 0);

label = gtk_label_new("Print interval");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);
label = gtk_label_new("Average start");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);
label = gtk_label_new("Average interval");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);
label = gtk_label_new("Roll interval");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);
/* Backup options */
label = gtk_label_new("Backup interval");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);

vbox = gtk_vbox_new(TRUE, 2);
gtk_box_pack_start (GTK_BOX (hbox), vbox, TRUE, TRUE, 0);

adj = (GtkAdjustment *) gtk_adjustment_new
                        (model->moldy.print_int, 0, 1e6, 1, 1, 0);
button = gtk_spin_button_new(adj, 0, 0);
gtk_spin_button_set_wrap(GTK_SPIN_BUTTON(button), FALSE);
gtk_box_pack_start(GTK_BOX(vbox), button, FALSE, FALSE, 0);
gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT (adj), "value_changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) adj);
g_object_set_data(G_OBJECT(adj), "key", (gpointer) PRINTINT);
g_object_set_data(G_OBJECT(adj), "ptr", (gpointer) model);

adj = (GtkAdjustment *) gtk_adjustment_new
                        (model->moldy.av_start, 0, 1e6, 1, 1, 0);
button = gtk_spin_button_new(adj, 0, 0);
gtk_spin_button_set_wrap(GTK_SPIN_BUTTON(button), FALSE);
gtk_box_pack_start(GTK_BOX(vbox), button, FALSE, FALSE, 0);
gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT (adj), "value_changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) adj);
g_object_set_data(G_OBJECT(adj), "key", (gpointer) AVSTART);
g_object_set_data(G_OBJECT(adj), "ptr", (gpointer) model);

adj = (GtkAdjustment *) gtk_adjustment_new
                        (model->moldy.av_int, 0, 1e6, 1, 1, 0);
button = gtk_spin_button_new(adj, 0, 0);
gtk_spin_button_set_wrap(GTK_SPIN_BUTTON(button), FALSE);
gtk_box_pack_start(GTK_BOX(vbox), button, FALSE, FALSE, 0);
gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT (adj), "value_changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) adj);
g_object_set_data(G_OBJECT(adj), "key", (gpointer) AVINT);
g_object_set_data(G_OBJECT(adj), "ptr", (gpointer) model);

adj = (GtkAdjustment *) gtk_adjustment_new
                        (model->moldy.roll_int, 0, 1e6, 1, 1, 0);
button = gtk_spin_button_new(adj, 0, 0);
gtk_spin_button_set_wrap(GTK_SPIN_BUTTON(button), FALSE);
gtk_box_pack_start(GTK_BOX(vbox), button, FALSE, FALSE, 0);
gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT (adj), "value_changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) adj);
g_object_set_data(G_OBJECT(adj), "key", (gpointer) ROLLINT);
g_object_set_data(G_OBJECT(adj), "ptr", (gpointer) model);

adj = (GtkAdjustment *) gtk_adjustment_new
                        (model->moldy.backup_int, 0, 1e6, 1, 1, 0);
button = gtk_spin_button_new(adj, 0, 0);
gtk_spin_button_set_wrap(GTK_SPIN_BUTTON(button), FALSE);
gtk_box_pack_start(GTK_BOX(vbox), button, FALSE, FALSE, 0);
gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT (adj), "value_changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) adj);
g_object_set_data(G_OBJECT(adj), "key", (gpointer) BACKINT);
g_object_set_data(G_OBJECT(adj), "ptr", (gpointer) model);

/* Dump level check radio */
frame = gtk_frame_new(" Dump contents ");
gtk_box_pack_start(GTK_BOX(page), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);

hbox = gtk_hbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(frame), hbox);

button = gtk_check_button_new_with_label("Coordinates");
if( model->moldy.dump_level & 1)
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
gtk_box_pack_start (GTK_BOX (hbox), button, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (dump_level_toggle),
                         (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) 1);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);

button = gtk_check_button_new_with_label("Velocities");
if( model->moldy.dump_level & 2)
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
gtk_box_pack_start (GTK_BOX (hbox), button, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (dump_level_toggle),
                         (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) 2);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);

button = gtk_check_button_new_with_label("Accelerations");
if( model->moldy.dump_level & 4)
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
gtk_box_pack_start (GTK_BOX (hbox), button, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (dump_level_toggle),
                         (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) 4);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);

button = gtk_check_button_new_with_label("Forces");
if( model->moldy.dump_level & 8)
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
gtk_box_pack_start (GTK_BOX (hbox), button, TRUE, TRUE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                          GTK_SIGNAL_FUNC (dump_level_toggle),
                         (gpointer) button);
g_object_set_data(G_OBJECT(button), "key", (gpointer) 8);
g_object_set_data(G_OBJECT(button), "ptr", (gpointer) model);

/* Dump files */
dump_frame = gtk_frame_new(" Dump files ");
gtk_box_pack_start(GTK_BOX(page), dump_frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(dump_frame), PANEL_SPACING);

hbox = gtk_hbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(dump_frame), hbox);

/* left vbox */
vbox1 = gtk_vbox_new(TRUE, 2);
gtk_box_pack_start (GTK_BOX (hbox), vbox1, TRUE, TRUE, 0);

label = gtk_label_new("Dump start");
gtk_box_pack_start (GTK_BOX (vbox1), label, TRUE, TRUE, 0);
label = gtk_label_new("Dump interval");
gtk_box_pack_start (GTK_BOX (vbox1), label, TRUE, TRUE, 0);
label = gtk_label_new("Dumps per file");
gtk_box_pack_start (GTK_BOX (vbox1), label, TRUE, TRUE, 0);

vbox2 = gtk_vbox_new(TRUE, 2);
gtk_box_pack_start (GTK_BOX (hbox), vbox2, TRUE, TRUE, 0);

adj = (GtkAdjustment *) gtk_adjustment_new
                        (model->moldy.dump_begin, 0, 1e6, 1, 1, 0);
button = gtk_spin_button_new(adj, 0, 0);
gtk_spin_button_set_wrap(GTK_SPIN_BUTTON(button), FALSE);
gtk_box_pack_start(GTK_BOX(vbox2), button, FALSE, FALSE, 0);
gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT (adj), "value_changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) adj);
g_object_set_data(G_OBJECT(adj), "key", (gpointer) DUMPSTART);
g_object_set_data(G_OBJECT(adj), "ptr", (gpointer) model);

adj = (GtkAdjustment *) gtk_adjustment_new
                        (model->moldy.dump_int, 0, 1e6, 1, 1, 0);
button = gtk_spin_button_new(adj, 0, 0);
gtk_spin_button_set_wrap(GTK_SPIN_BUTTON(button), FALSE);
gtk_box_pack_start(GTK_BOX(vbox2), button, FALSE, FALSE, 0);
gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT (adj), "value_changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) adj);
g_object_set_data(G_OBJECT(adj), "key", (gpointer) DUMPINT);
g_object_set_data(G_OBJECT(adj), "ptr", (gpointer) model);

adj = (GtkAdjustment *) gtk_adjustment_new
                        (model->moldy.maxdumps, 0, 1e6, 1, 1, 0);
button = gtk_spin_button_new(adj, 0, 0);
gtk_spin_button_set_wrap(GTK_SPIN_BUTTON(button), FALSE);
gtk_box_pack_start(GTK_BOX(vbox2), button, FALSE, FALSE, 0);
gtk_spin_button_set_numeric(GTK_SPIN_BUTTON(button), TRUE);
g_signal_connect(GTK_OBJECT (adj), "value_changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) adj);
g_object_set_data(G_OBJECT(adj), "key", (gpointer) NDUMPS);
g_object_set_data(G_OBJECT(adj), "ptr", (gpointer) model);

/* summary details for the model */
frame = gtk_frame_new(" Details ");
gtk_box_pack_start(GTK_BOX(box), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
hbox = gtk_hbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(frame),hbox);

/* left vbox */
vbox = gtk_vbox_new(TRUE, 2);
gtk_box_pack_start (GTK_BOX (hbox), vbox, TRUE, TRUE, 0);

label = gtk_label_new("Structure name");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);
label = gtk_label_new("Description");
gtk_box_pack_start (GTK_BOX (vbox), label, TRUE, TRUE, 0);

/* right vbox */
vbox = gtk_vbox_new(TRUE, 2);
gtk_box_pack_start (GTK_BOX (hbox), vbox, TRUE, TRUE, 0);

entry = gtk_entry_new_with_max_length(LINELEN);
gtk_entry_set_text(GTK_ENTRY(entry), model->basename);
gtk_box_pack_start (GTK_BOX (vbox), entry, TRUE, TRUE, 0);

/* event */
g_signal_connect(G_OBJECT(entry), "changed",
                   GTK_SIGNAL_FUNC(change_basename),
                   (gpointer) model);

entry = gtk_entry_new_with_max_length(WIDTH);
gtk_entry_set_text(GTK_ENTRY(entry), model->moldy.title);
gtk_box_pack_start (GTK_BOX (vbox), entry, TRUE, TRUE, 0);

/* event */
g_signal_connect(GTK_OBJECT (entry), "changed",
                          GTK_SIGNAL_FUNC (moldy_keyword),
                         (gpointer) entry);
g_object_set_data(G_OBJECT(entry), "key", (gpointer) TITLE);
g_object_set_data(G_OBJECT(entry), "ptr", (gpointer) model);


/* done */
gtk_widget_show_all(box);

/* if (!model->moldy.auto_cutoff)
  gtk_widget_set_sensitive(GTK_WIDGET(model->moldy.cutoff_frame), TRUE);
else
  gtk_widget_set_sensitive(GTK_WIDGET(model->moldy.cutoff_frame), FALSE);

if (model->moldy.dump_level)
  gtk_widget_set_sensitive(dump_frame, TRUE);
else
  gtk_widget_set_sensitive(dump_frame, FALSE); */

g_string_free(line, TRUE);
gui_model_select(model);
}
/**************************/
/* The MOLDY setup dialog */
/**************************/
void moldy_dialog(void)
{
gpointer dialog;
GtkWidget *window, *frame, *vbox;
struct model_pak *model;

model = sysenv.active_model;
if (!model)
  return;

if( model->id == MORPH || model->periodic == 2)
  return;

moldy_files_init(model);

/* request a moldy dialog */
dialog = dialog_request(MOLDY, "MOLDY", NULL, moldy_cleanup, model);
if (!dialog)
  return;
window = dialog_window(dialog);

frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(FALSE,0);
gtk_container_add(GTK_CONTAINER(frame), vbox);
gui_moldy_widget(vbox, dialog);

/* terminating buttons */
gui_stock_button(GTK_STOCK_EXECUTE, run_moldy, model,
                   GTK_DIALOG(window)->action_area);

gui_stock_button(GTK_STOCK_CLOSE, dialog_destroy, dialog,
                   GTK_DIALOG(window)->action_area);

/* done */
gtk_widget_show_all(window);
}
