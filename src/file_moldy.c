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
* FILE_MOLDY.c   Functions for reading and writing Moldy input files      *
*                                                                         *
* moldy_init()          initialize moldy variables			  *
* pots_free()      	free potential structures			  *
* site_free()      	free site structures			          *
* moldyspec_free()    	free species structures				  *
* compare_sites()       sort sites                                        *
* compare_species()     sort species by site (monatomic) or name (poly)   *
* moldy_files_init()    create valid input and output file names          *
* moldy_data_free()     remove moldy from memory			  *
* match_key()		match keywords from control file	          *
* calc_units()          calculate units				          *
* random_species()      choose from species list                          *
* skew_start()		set up configuration using skew start method	  *
* add_site_element()    add site to global element list			  *
* create_sitetype()     create site type and add to list		  *
* new_site()		create site and add to species' site list	  *	
* create_species()	create new species and add to list		  *
* create_pot()		create potential parameter entry and add to list  *
* read_moldy()		read moldy input (control and/or sys-spec) file	  *
* find_unique_sites()	identify number and type of different sites	  *
* match_molecules()	rotate molecules and compare site positions/types *
* calc_species()	identify species that molecule belongs to 	  *
* read_moldy_restart()  use Moldy util ransub to read in restart file     *
* read_moldy_potentials()	read Moldy potential parameters from file *
* write_moldy_potentials()	write potential parameters to file 	  *
* write_moldy_control() write control file for Moldy			  *
* write_moldy_sys()	write sys-spec file for Moldy			  *
*                                                                         *
* Functions for performing rotation-related calcs:                        *
* calc_centre()         calculate geometric centre of molecule            *
* calc_cofm()           calculate centre of mass of molecule              *
* moments()             calculate moments of inertia                      *
* eigen()               calculate eigenvalues and vectors from matrix     *
* calc_pfc()            calculate principal frame coordinates             *
* jacobi()              jacobi method of matrix diagonalization           *
* quat_fit()            calculate quaternion between two molecules        *
**************************************************************************/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "gdis.h"
#include "coords.h"
#include "model.h"
#include "file.h"
#include "matrix.h"
#include "interface.h"
#include "moldy.h"
/* #include "vector.h" */
#include "parse.h"
#include "task.h"
#include "space.h"
#include "quaternion.h"

#define TMPNAME ".gdistmp"
#define MATCH_KEY 1
#define FLOATING_POINT 1e-12

/* main structures */
extern struct sysenv_pak sysenv;
extern const pots_pak potspec[];
void dialog_destroy_type(gint);
void update_pot_dialog(struct model_pak *);

/***********************************/
/* MOLDY structure initialization */
/***********************************/
void moldy_init(gpointer data)
{
struct moldy_pak *moldy = data;

moldy->no_exec = TRUE;
moldy->set_framework = FALSE;
moldy->title = g_strdup("");
moldy->cpu_limit = 1e20;
moldy->temperature = 0.0;
moldy->pressure = 0.0;
moldy->mass_parm = 100.0;
moldy->deltat = 0.005;
moldy->nsteps = 0;
moldy->strict = FALSE;
moldy->surf = FALSE;
moldy->latt_start = FALSE;
moldy->density = 1.0;
moldy->const_temp = 0;
moldy->const_press = 0;
moldy->strain_mask = 200;
moldy->subcell = 0.0;
moldy->num_mols = 1; /* Used for skew start settings */
moldy->species = NULL;
moldy->pot_type = -1;
moldy->pots = NULL;
moldy->energy_unit = KJMOL;
moldy->time_unit = TIME_KJ;
moldy->charge_unit = ELCHARGE;
moldy->length_unit = ANGST;
moldy->mass_unit = AMU;
moldy->ewald = TRUE;
moldy->auto_cutoff = TRUE;
moldy->real_only = FALSE;
moldy->cutoff = 0.0;
moldy->alpha = 0.0;
moldy->kcutoff = 0.0;
moldy->ewald_accuracy = 1.013e-5;
moldy->rdf_begin = 1e6;
moldy->rdf_int = 20;
moldy->rdf_out = 5000;
moldy->rdf_limit = 10.0;
moldy->nbins = 100;
moldy->scale_int = 10;
moldy->scale_end = 1e6;
moldy->scale_options = 0;
moldy->dump_level = 0;
moldy->dump_begin = 1;
moldy->dump_int = 20;
moldy->maxdumps = 250;
moldy->av_start = 1001;
moldy->av_reset = FALSE;
moldy->av_int = 5000;
moldy->roll_int = 10;
moldy->print_int = 10;
moldy->backup_int = 500;
moldy->ttmass = 100.0;
moldy->rtmass = 100.0;
moldy->seed = 1234567;
moldy->page_width = 132;
moldy->page_length = 44;
moldy->xdr = TRUE;
moldy->text_save = FALSE;
moldy->control_file = g_strdup("");
moldy->sysspec_file = g_strdup("");
moldy->lib_file = g_strdup("");
moldy->lib_dir = g_strdup(sysenv.cwd);
moldy->out_file = g_strdup("");
moldy->save_file = g_strdup("");
moldy->restart_file = g_strdup("");
moldy->restart_dir = g_strdup("");
moldy->dump_file = g_strdup("");
moldy->backup_file = g_strdup("MDBACKUP");
moldy->temp_file = g_strdup("MDTEMPX");
}

/******************************/
/* free pots structure        */
/******************************/
void pots_free(gpointer data)
{
struct pot_pak *pot = data;

g_free(pot->site1_label);
g_free(pot->site2_label);
g_free(pot);
}

/*******************************/
/* free site structures        */
/*******************************/
void site_free(gpointer data)
{
/* struct site_pak *site = data; */

/*  g_slist_free(site->bonds); */
}
/***********************************/
/* free species type structures    */
/***********************************/
void moldyspec_free(gpointer data)
{
struct moldyspec_pak *species = data;
GSList *list;
struct site_pak *site;

g_free(species->name);

for (list=species->sites; list; list=g_slist_next(list))
  {
  site=list->data;
  site_free(site);
  }
free_slist(species->sites);

g_free(species);
}

#define DEBUG_CALC_CENTRE 0
/* Calculate geometric centre of molecule */
void calc_centre(gdouble (*sites)[3], gdouble *centre, gint n)
{
gint i,j;
gdouble min, max;

if (n < 1)
  return;

VEC3SET(centre, 0.0, 0.0, 0.0);

for (j = 3; j--;)
  {
  min = sites[0][j];
  max = sites[0][j];
  for (i = 1; i < n; i++)
    {
    min = MINIMUM(min, sites[i][j]);
    max = MAXIMUM(max, sites[i][j]);
    }
  centre[j] = (min + max)/2.0;
  }

#if DEBUG_CALC_CENTRE
  P3VEC("Centroid =",&centre[0]);
#endif
}

/* Calculate centre of mass of molecule */
#define DEBUG_CALC_COFM 0
void calc_cofm(struct model_pak *data, struct mol_pak *mol, gdouble *cofm)
{
gint i; /* Counter */
gdouble vec[3]; /* Position of current core */
gdouble mass, total_mass=0.0;
GSList *list;
struct core_pak *core;

VEC3SET(cofm, 0.0, 0.0, 0.0);

for (list=mol->cores ; list ; list=g_slist_next(list))
  {
  core = (struct core_pak *) list->data;
  if (core->status & (DELETED | HIDDEN))
    continue;
  ARR3SET(vec, core->x);
  vecmat(data->latmat, vec); /* Convert to Cartesian coords */
  mass = atom_mass(core);
  for(i=3; i--;)
    cofm[i] += mass*vec[i];
  total_mass += mass;

#if DEBUG_CALC_COFM
  printf("atom %s %f\n",core->atom_label, mass);
#endif
  }

if ( total_mass > 0.0)
  {
  VEC3MUL(cofm, 1.0/total_mass);
  }
else
  {
  ARR3SET(cofm, mol->centroid);
  vecmat(data->latmat, cofm); /* Convert to Cartesian coords */
  }

#if DEBUG_CALC_COFM
printf("C of M %f %f %f\n", cofm[0], cofm[1], cofm[2]);
#endif
}

/* Calculate moments of inertia from arrays of positions and masses */
void moments(gint n, gdouble *mass, gdouble **sites, gdouble *inertia)
{
gint i;
gdouble *x = sites[0];
gdouble *y = sites[1];
gdouble *z = sites[2];
gdouble sxx=0.0, syy=0.0, szz=0.0;
gdouble sxy=0.0, sxz=0.0, syz=0.0;

for (i=0; i<n; i++)
  {
  sxx += mass[i]*x[i]*x[i];
  syy += mass[i]*y[i]*y[i];
  szz += mass[i]*z[i]*z[i];
  sxy += mass[i]*x[i]*y[i];
  sxz += mass[i]*x[i]*z[i];
  syz += mass[i]*y[i]*z[i];
  }

inertia[0] = syy+szz; inertia[1] = -sxy; inertia[2] = sxx+szz;
inertia[3] = -sxz; inertia[4] = -syz; inertia[5] = sxx + syy;
}

/* Diagonalize a real symmetric matrix and calculate eigenvalues and vectors.
   The matrix should be in a 1-D array as the lower diagonal portion
   of the matrix, by rows. Converted from Tom Shattuck's Javascript code.
   a = matrix to diagonalize
   r = matrix of eigenvectors
   n = dimension of the square symmetric matrix
   p = factors for rotating atom positions after sorting */
void eigen(gdouble *a, gdouble *r, gint n, gdouble *p)
{
gint i, j, k;
gdouble x = 0, y = 0;
gint l = 0, m = 0, iq, ij, jq;
gint ll, il, im, mq, lq, lm, mm;
gdouble anorm = 0.0, anrmx, range = 1.0e-12;
gint redo;
gdouble thr, sinx, sinx2, cosx, cosx2, sincs;
gint pass, ia, imq, ilq, imr, ilr;

for (i=0; i<n; i++) {
  for (j=0; j<n; j++) {
    ij = i*n+j;
    if (i == j)
      r[ij] = 1.0;
    else
      r[ij] = 0.0;
  }
}
for (i=0; i<n; i++) {
  for (j=0; j<n; j++) {
    if (i != j) {
      ia = i+ (j*j+j)/2;
      anorm = anorm + a[ia]*a[ia] ;
    }
  }
}
if (anorm > 0) {
   anorm = sqrt(2)* sqrt(anorm);
   anrmx = anorm*range/n;

/* initialize indicators and compute threshold, thr */
redo = 0;
thr = anorm;
/* compare threshold with final norm */
while (thr > anrmx)
  {
  thr = thr/n;
  for (pass = 0; pass < 1000; pass++)
    {
    l=0;
/*  test for l=second from last column */
      while (l < (n-1))
        {
        m = l+1;
/*   test for m=last column */
        while (m < n)
          {
          mq = m*(m + 1)/2;
          lq = l*(l+1)/2;
          lm = l + mq;
          if (fabs(a[lm]) >= thr)
            {
            redo = 1;
            ll = l + lq;
            mm = m + mq;
            x = 0.5*(a[ll] - a[mm]);
            y = -a[lm]/ sqrt(a[lm]*a[lm] + x*x);
            if (x < 0)
              y = -1.0*y;
            sinx = y/sqrt(2.0*(1.0+(sqrt(1.0 - y*y))));
            sinx2 = sinx*sinx;
            cosx = sqrt(1.0 - sinx2);
            cosx2 = cosx*cosx;
            sincs = sinx*cosx;
/*      rotate l and m columns */
            ilq = n*l;
            imq = n*m;
            for (i=0; i<n; i++)
              {
              iq = i*(i+1)/2;
              if (i != l)
                {
                if (i != m)
	          {
                  if (i < m) { im =i + mq; } else { im = m + iq; }
                  if (i < l) { il = i + lq; } else { il = l + iq; }
                  x = a[il]*cosx - a[im]*sinx;
                  a[im] = a[il]*sinx + a[im]*cosx;
                  a[il] = x;
                  }
                }
                ilr = ilq+i;
                imr = imq+i;
                x = r[ilr]*cosx - r[imr]*sinx;
                r[imr] = r[ilr]*sinx+r[imr]*cosx;
                r[ilr] = x;
              }
            x = 2.0*a[lm]*sincs;
            y = a[ll]*cosx2 + a[mm]*sinx2 - x;
            x = a[ll]*sinx2 + a[mm]*cosx2 + x;
            a[lm] = (a[ll]-a[mm])*sincs+ a[lm]*(cosx2-sinx2);
            a[ll] = y;
            a[mm] = x;
            }

/* test for completion */
          m++;
          }
        l++;
        }
       if (! redo)
         break;
       redo = 0;
      }
    }
  }

/* sort eigen values and vectors */
  for (i = 0; i < n; i++)
    {
    ll = i+(i*i+i)/2;
    iq = n*i;
    jq = n*(i-1);
    for (j = i; j < n; j++) {
      jq += n;
      mm = j+(j*j+j)/2;
      if (a[ll] < a[mm]) {
        x = a[ll];
        a[ll] = a[mm];
        a[mm] = x;
        for (k = 0; k < n; k++) {
          ilr = iq+k;
          imr = jq+k;
          x = r[ilr];
          r[ilr] = r[imr];
          r[imr] = x;
          if (k != i && k != j)
            p[k] *= -1;
        }
      }
    }
  }
}

/* Calculate principal frame coordinates of molecule */
#define DEBUG_CALC_PFC 0
void calc_pfc(struct model_pak *data, struct mol_pak *mol, gdouble **p_f_sites, gint num_atoms)
{
gint i, isite;
gdouble mass[num_atoms];
gdouble inertia[6], rot[9];
gdouble pfs[3], cofm[3], factor[3];
GSList *list;
struct core_pak *core;

calc_cofm(data, mol, cofm);

isite=0;
for (list = mol->cores ; list ; list=g_slist_next(list))
  {
  core = (struct core_pak *) list->data;
  if (core->status & (DELETED | HIDDEN))
    continue;
  mass[isite] = atom_mass(core);
  isite++;
}

for (isite = 0; isite < num_atoms; isite++)
  {
  VEC3SET(pfs, p_f_sites[0][isite],p_f_sites[1][isite],p_f_sites[2][isite]);

  if (data->periodic == 3)
    vecmat(data->latmat, pfs); /* Conv frac coords to real coords */

  for(i=0; i < 3; i++)
    p_f_sites[i][isite] = pfs[i] - cofm[i];
  }

if (num_atoms > 1)
  {
  VEC3SET(&inertia[0], 0.0, 0.0, 0.0);
  VEC3SET(&inertia[3], 0.0, 0.0, 0.0);
  VEC3SET(&factor[0], 1.0, 1.0, 1.0);

  moments(num_atoms, mass, p_f_sites, inertia);

  eigen(inertia, rot, 3, factor);

/* rotate to principal frame */
  for (isite = 0; isite < num_atoms; isite++)
    {
    VEC3SET(pfs, p_f_sites[0][isite], p_f_sites[1][isite], p_f_sites[2][isite]);
    vecmat(rot, pfs);

    for(i=0; i < 3; i++)
      p_f_sites[i][isite] = factor[i] * pfs[i];
    }

#if DEBUG_CALC_PFC
  P3VEC("I =",&inertia[0]);
  P3VEC("   ",&inertia[3]);
  P3MAT("rotation matrix: ", rot);
#endif
  }
}

/*********************************************************
 JACOBI
 Jacobi diagonalizer with sorted output.
 Based on code by Jan Labanowski/David Heisterberg.
 a - input: matrix to diagonalize
 v - output: eigenvectors
 d - output: eigenvalues
 nrot - input: maximum number of sweeps
**********************************************************/
void jacobi(gdouble a[4][4], gdouble *d, gdouble v[4][4], gint nrot)
{
gdouble onorm, dnorm;
gdouble b, dma, q, t, c, s;
gdouble atemp, vtemp, dtemp;
gint i, j, k, l;

 for (j = 0; j < 4; j++)
   {
   for (i = 0; i < 4; i++)
     v[i][j] = 0.0;
   v[j][j] = 1.0;
   d[j] = a[j][j];
   }

 for (l = 1; l < nrot+1; l++)
   {
   dnorm = 0.0;
   onorm = 0.0;
   for (j = 0; j < 4; j++)
     {
     dnorm = dnorm + fabs(d[j]);
     for (i = 0; i < j; i++)
       onorm = onorm + fabs(a[i][j]);
     }
   if ((onorm/dnorm) <= 1.0e-12) goto Exit_now;
   for (j = 1; j < 4; j++)
     {
     for (i = 0; i < j; i++)
       {
       b = a[i][j];
       if (fabs(b) > 0.0)
         {
         dma = d[j] - d[i];
         if ((fabs(dma) + fabs(b)) <=  fabs(dma))
           {
           t = b / dma;
           }
         else
           {
           q = 0.5 * dma / b;
           t = 1.0/(fabs(q) + sqrt(1.0+q*q));
           if (q < 0.0)
             t = -t;
           }
         c = 1.0/sqrt(t * t + 1.0);
         s = t * c;
         a[i][j] = 0.0;
         for (k = 0; k < i; k++)
           {
           atemp = c * a[k][i] - s * a[k][j];
           a[k][j] = s * a[k][i] + c * a[k][j];
           a[k][i] = atemp;
           }
         for (k = i+1; k < j; k++)
           {
           atemp = c * a[i][k] - s * a[k][j];
           a[k][j] = s * a[i][k] + c * a[k][j];
           a[i][k] = atemp;
           }
         for (k = j+1; k < 4; k++)
           {
           atemp = c * a[i][k] - s * a[j][k];
           a[j][k] = s * a[i][k] + c * a[j][k];
           a[i][k] = atemp;
           }
         for (k = 0; k < 4; k++)
           {
           vtemp = c * v[k][i] - s * v[k][j];
           v[k][j] = s * v[k][i] + c * v[k][j];
           v[k][i] = vtemp;
           }
         dtemp = c*c*d[i] + s*s*d[j] - 2.0*c*s*b;
         d[j] = s*s*d[i] + c*c*d[j] +  2.0*c*s*b;
         d[i] = dtemp;
         }  /* end if */
       } /* end for i */
     } /* end for j */
   } /* end for l */

Exit_now:
 nrot = l;
 for (j = 0; j < 3; j++)
   {
   k = j;
   dtemp = d[k];
   for (i = j+1; i < 4; i++)
     {
     if (d[i] < dtemp)
       {
       k = i;
       dtemp = d[k];
       }
     }
   if (k > j)
     {
     d[k] = d[j];
     d[j] = dtemp;
     for (i = 0; i < 4; i++)
       {
       dtemp = v[i][k];
       v[i][k] = v[i][j];
       v[i][j] = dtemp;
       }
     }
   }
}
/**************************************************
 Quaternion fitting
 Based on code by Jan Labanowski/David Heisterberg.
 Find the quaternion, q, that minimizes

   |qTXq - Y| ^ 2

 This is equivalent to maximizing Re (qTXTqY).

 This is equivalent to finding the largest eigenvalue and corresponding
 eigenvector of the matrix

 [A2   AUx  AUy  AUz ]
 [AUx  Ux2  UxUy UzUx]
 [AUy  UxUy Uy2  UyUz]
 [AUz  UzUx UyUz Uz2 ]

 where

   A2   = Xx Yx + Xy Yy + Xz Yz
   Ux2  = Xx Yx - Xy Yy - Xz Yz
   Uy2  = Xy Yy - Xz Yz - Xx Yx
   Uz2  = Xz Yz - Xx Yx - Xy Yy
   AUx  = Xz Yy - Xy Yz
   AUy  = Xx Yz - Xz Yx
   AUz  = Xy Yx - Xx Yy
   UxUy = Xx Yy + Xy Yx
   UyUz = Xy Yz + Xz Yy
   UzUx = Xz Yx + Xx Yz

 INPUT
   n      - number of points
   x      - fitted molecule coordinates
   y      - reference molecule coordinates
   nr     - max number of jacobi sweeps

 OUTPUT
   q      - the best-fit quaternion
*****************************************************/
void quat_fit(struct model_pak *data, struct mol_pak *molfit,
      struct moldyspec_pak *species, gdouble *q, gint nr)
{
 gint i, j, isite=0;
 gint n = g_slist_length(species->sites);
 gdouble site_fit[3];
 gdouble centre_x[3];
 gdouble x[n][3];
 gdouble y[n][3];
 gdouble xxyx, xxyy, xxyz;
 gdouble xyyx, xyyy, xyyz;
 gdouble xzyx, xzyy, xzyz;
 gdouble c[4][4], v[4][4];
 gdouble d[4];
 GSList *list, *slist;
 struct core_pak *core;
 struct site_pak *site;

 xxyx = 0.0;
 xxyy = 0.0;
 xxyz = 0.0;
 xyyx = 0.0;
 xyyy = 0.0;
 xyyz = 0.0;
 xzyx = 0.0;
 xzyy = 0.0;
 xzyz = 0.0;

 isite=0;
 slist = species->sites;
 for (list=molfit->cores ; list ; list=g_slist_next(list))
   {
   core = (struct core_pak *) list->data;
   if (core->status & (DELETED | HIDDEN))
     continue;

   for (i = 3; i--;)
     site_fit[i] = core->x[i];
   if (data->periodic == 3)
     vecmat(data->latmat, site_fit);


   site = (struct site_pak *) slist->data;

   for (i = 3; i--;)
     {
     x[isite][i] = site_fit[i];
     y[isite][i] = site->x[i];
     }
   slist = g_slist_next(slist);
   isite++;
   }

 calc_cofm(data, molfit, centre_x);

 for(i = 0; i < n; i++)
   for (j = 0; j < 3; j++)
     x[i][j] -= centre_x[j];

/* generate the upper triangle of the quadratic form matrix */
 for (i = 0; i < n; i++)
   {
   xxyx = xxyx + x[i][0] * y[i][0];
   xxyy = xxyy + x[i][0] * y[i][1];
   xxyz = xxyz + x[i][0] * y[i][2];
   xyyx = xyyx + x[i][1] * y[i][0];
   xyyy = xyyy + x[i][1] * y[i][1];
   xyyz = xyyz + x[i][1] * y[i][2];
   xzyx = xzyx + x[i][2] * y[i][0];
   xzyy = xzyy + x[i][2] * y[i][1];
   xzyz = xzyz + x[i][2] * y[i][2];
   }

 for(i = 0; i < 4; i++)
   for(j = 0; j < 4; j++)
      c[i][j] = 0.0;

 c[0][0] = xxyx + xyyy + xzyz;

 c[0][1] = xzyy - xyyz;
 c[1][1] = xxyx - xyyy - xzyz;

 c[0][2] = xxyz - xzyx;
 c[1][2] = xxyy + xyyx;
 c[2][2] = xyyy - xzyz - xxyx;

 c[0][3] = xyyx - xxyy;
 c[1][3] = xzyx + xxyz;
 c[2][3] = xyyz + xzyy;
 c[3][3] = xzyz - xxyx - xyyy;

/* diagonalize c */
 jacobi (c, d, v, nr);

/* extract the desired quaternion */
 VEC4SET(q,v[0][3],v[1][3],v[2][3],v[3][3]);
}

#define DEBUG_COMPARE_SITES 0
gint compare_sites(gpointer ptr_site1, gpointer ptr_site2)
{
struct atomtype_pak *s1 = ptr_site1, *s2 = ptr_site2;
struct model_pak *data;
struct elem_pak elem_data;
gint code1, code2;
gdouble q1, q2, radius1, radius2;

data = sysenv.active_model;

/* Return if identical */
if (s1->id == s2->id)
  return(0);

code1 = elem_symbol_test (s1->label);
q1 = s1->charge;
get_elem_data (code1, &elem_data, data);
radius1 = elem_data.vdw;

code2 = elem_symbol_test (s2->label);
q2 = s2->charge;
get_elem_data (code2, &elem_data, data);
radius2 = elem_data.vdw;

#if DEBUG_COMPARE_SITES
puts("-------");
printf("Atom %d: %s %f %f\n", code1, s1->label, q1, radius1);
printf("Atom %d: %s %f %f\n", code2, s2->label, q2, radius2);
#endif

/* Cations before anions */
if (q1 < 0 && q2 > 0)
  return (1);
if (q1 > 0 && q2 < 0)
  return (-1);

/* Low charge before high charge */
if (fabs (q1) > fabs (q2))
  return (1);
else if (fabs (q1) < fabs (q2))
  return (-1);

/* Large ions before small ions */
if (radius1 < radius2)
  return(1);

else if (radius1 > radius2)
  return(-1);

/* Alphabetical order */
if (g_ascii_strcasecmp(s1->label, s2->label) > 0)
  return(1);
else if (g_ascii_strcasecmp(s1->label, s2->label) < 0)
  return(-1);

/* Sites must be identical */
return(0);
}
gint compare_species(gpointer ptr_species1, gpointer ptr_species2)
{
struct moldyspec_pak *s1 = ptr_species1, *s2 = ptr_species2;
struct mol_pak *mol1, *mol2;
struct core_pak *core1, *core2;

mol1  = (s1->mols)->data;
mol2  = (s2->mols)->data;
core1 = (mol1->cores)->data;
core2 = (mol2->cores)->data;

/* If monatomic, sort by atom type */
if (g_slist_length(s1->sites) == 1 &&
    g_slist_length(s2->sites) == 1)
   return(compare_atoms(core1, core2));

/* Otherwise sort by molecule */
return(compare_molecules(s1, s2));

return(0);
}

/**********************************************/
/* enforce valid input/output Moldy filenames */
/**********************************************/
void moldy_files_init(struct model_pak *model)
{
gchar *basename;

g_assert(model != NULL);
g_assert(model->basename != NULL);

/* Remove spaces from basename */
basename = g_strdup(model->basename);
parse_space_replace(basename,'_');

/* Create file names from model name (unless already set) */
if (g_ascii_strcasecmp(model->moldy.control_file, "") == 0)
  {
  g_free(model->moldy.control_file);
  model->moldy.control_file = g_strdup_printf("%s.min", basename);
  }

if (g_ascii_strcasecmp(model->moldy.sysspec_file, "") == 0)
  {
  g_free(model->moldy.sysspec_file);
  model->moldy.sysspec_file = g_strdup_printf("%s.min", basename);
  }

if (g_ascii_strcasecmp(model->moldy.out_file, "") == 0)
  {
  g_free(model->moldy.out_file);
  model->moldy.out_file = g_strdup_printf("%s.mdo", basename);
  }

if (g_ascii_strcasecmp(model->moldy.dump_file, "") == 0)
  {
  g_free(model->moldy.dump_file);
  model->moldy.dump_file = g_strdup_printf("%s%%d.dat", basename);
  }

if (g_ascii_strcasecmp(model->moldy.save_file, model->moldy.restart_file) == 0)
  {
  g_free(model->moldy.save_file);
  model->moldy.save_file = g_strdup_printf("%s.sav", basename);
  }
}

/***********************************/
/* free memory used for moldy data */
/***********************************/
#define DEBUG_FREE_MOLDY_DATA 0
void moldy_data_free(struct model_pak *model)
{
GSList *list;

g_assert(model != NULL);

#if DEBUG_FREE_MOLDY_DATA
printf("freeing moldy data...\n");
#endif
g_free(model->moldy.title);
g_free(model->moldy.control_file);
g_free(model->moldy.sysspec_file);
g_free(model->moldy.save_file);
g_free(model->moldy.restart_file);
g_free(model->moldy.restart_dir);
g_free(model->moldy.dump_file);
g_free(model->moldy.lib_file);
g_free(model->moldy.lib_dir);
g_free(model->moldy.out_file);
g_free(model->moldy.backup_file);
g_free(model->moldy.temp_file);

for (list = model->moldy.species; list; list=g_slist_next(list))
  moldyspec_free(list->data);
for (list = model->moldy.pots; list; list=g_slist_next(list))
  pots_free(list->data);

model->moldy.title = NULL;
model->moldy.control_file = NULL;
model->moldy.sysspec_file = NULL;
model->moldy.save_file = NULL;
model->moldy.restart_file = NULL;
model->moldy.restart_dir = NULL;
model->moldy.lib_file = NULL;
model->moldy.lib_dir = NULL;
model->moldy.backup_file = NULL;
model->moldy.temp_file = NULL;

#if DEBUG_FREE_MOLDY_DATA
printf("done.\n");
#endif
}
gint match_key(struct model_pak *data, gchar *name, gchar *value)
{

if (g_ascii_strcasecmp(name, "title") == 0)
  {
  data->moldy.title = g_strdup(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "sys-spec-file") == 0)
  {
  data->moldy.sysspec_file = g_strdup(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "new-sys-spec") == 0)
  {
  return(0);
  }
if (g_ascii_strcasecmp(name, "nsteps") == 0)
  {
  data->moldy.nsteps = atoi(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "step") == 0)
  {
  data->moldy.deltat = g_ascii_strtod(value, NULL);
  return(0);
  }
if (g_ascii_strcasecmp(name, "cpu-limit") == 0)
  {
  data->moldy.cpu_limit = g_ascii_strtod(value, NULL);
  return(0);
  }
if (g_ascii_strcasecmp(name, "scale-options") == 0 ||
       g_ascii_strcasecmp(name, "therm-options") == 0)
  {
  data->moldy.scale_options = atoi(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "temperature") == 0)
  {
  data->moldy.temperature=g_ascii_strtod(value, NULL);
  return(0);
  }
if (g_ascii_strcasecmp(name, "pressure") == 0)
  {
  data->moldy.pressure=g_ascii_strtod(value, NULL);
  return(0);
  }
if (g_ascii_strcasecmp(name, "w") == 0)
  {
  data->moldy.mass_parm=g_ascii_strtod(value, NULL);
  return(0);
  }
if (g_ascii_strcasecmp(name, "strict-cutoff") == 0)
  {
  data->moldy.strict=atoi(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "surface-dipole") == 0)
  {
  data->moldy.surf=atoi(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "lattice-start") == 0)
  {
  data->moldy.latt_start=atoi(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "density") == 0)
  {
  data->moldy.density=g_ascii_strtod(value, NULL);
  return(0);
  }
if (g_ascii_strcasecmp(name, "const-temp") == 0)
  {
  data->moldy.const_temp=atoi(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "const-pressure") == 0)
  {
  data->moldy.const_press=atoi(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "strain-mask") == 0)
  {
  data->moldy.strain_mask=atoi(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "subcell") == 0)
  {
  data->moldy.subcell=g_ascii_strtod(value, NULL);
  return(0);
  }
if (g_ascii_strcasecmp(name, "ewald-accuracy") == 0)
  {
  data->moldy.ewald_accuracy=g_ascii_strtod(value, NULL);
  return(0);
  }
if (g_ascii_strcasecmp(name, "time-unit") == 0)
  {
  data->moldy.time_unit=g_ascii_strtod(value, NULL);
  return(0);
  }
if (g_ascii_strcasecmp(name, "charge-unit") == 0)
  {
  data->moldy.charge_unit=g_ascii_strtod(value, NULL);
  return(0);
  }
if (g_ascii_strcasecmp(name, "length-unit") == 0)
  {
  data->moldy.length_unit=g_ascii_strtod(value, NULL);
  return(0);
  }
if (g_ascii_strcasecmp(name, "mass-unit") == 0)
  {
  data->moldy.mass_unit=g_ascii_strtod(value, NULL);
  return(0);
  }
if (g_ascii_strcasecmp(name, "cutoff") == 0)
  {
  data->moldy.cutoff=g_ascii_strtod(value, NULL);
  data->moldy.auto_cutoff = FALSE;
  return(0);
  }
if (g_ascii_strcasecmp(name, "alpha") == 0)
  {
  data->moldy.alpha=g_ascii_strtod(value, NULL);
  data->moldy.auto_cutoff = FALSE;
  return(0);
  }
if (g_ascii_strcasecmp(name, "k-cutoff") == 0)
  {
  data->moldy.kcutoff=g_ascii_strtod(value, NULL);
  data->moldy.auto_cutoff = FALSE;
  return(0);
  }
if (g_ascii_strcasecmp(name, "begin-rdf") == 0)
  {
  data->moldy.rdf_begin=atoi(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "rdf-interval") == 0)
  {
  data->moldy.rdf_int=atoi(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "rdf-out") == 0)
  {
  data->moldy.rdf_out=atoi(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "rdf-limit") == 0)
  {
  data->moldy.rdf_limit=g_ascii_strtod(value, NULL);
  return(0);
  }
if (g_ascii_strcasecmp(name, "nbins") == 0)
  {
  data->moldy.nbins=atoi(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "scale-interval") == 0)
  {
  data->moldy.scale_int=atoi(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "scale-end") == 0)
  {
  data->moldy.scale_end=atoi(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "dump-level") == 0)
  {
  data->moldy.dump_level=atoi(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "begin-dump") == 0)
  {
  data->moldy.dump_begin=atoi(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "dump-interval") == 0)
  {
  data->moldy.dump_int=atoi(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "ndumps") == 0)
  {
  data->moldy.maxdumps=atoi(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "begin-average") == 0)
  {
  data->moldy.av_start=atoi(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "reset-averages") == 0)
  {
  data->moldy.av_reset=atoi(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "average-interval") == 0)
  {
  data->moldy.av_int=atoi(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "roll-interval") == 0)
  {
  data->moldy.roll_int=atoi(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "print-interval") == 0)
  {
  data->moldy.print_int=atoi(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "backup-interval") == 0)
  {
  data->moldy.backup_int=atoi(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "ttmass") == 0)
  {
  data->moldy.ttmass=g_ascii_strtod(value, NULL);
  return(0);
  }
if (g_ascii_strcasecmp(name, "rtmass") == 0)
  {
  data->moldy.rtmass=g_ascii_strtod(value, NULL);
  return(0);
  }
if (g_ascii_strcasecmp(name, "seed") == 0)
  {
  data->moldy.seed=atoi(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "page-width") == 0)
  {
  data->moldy.page_width=abs(atoi(value));
  return(0);
  }
if (g_ascii_strcasecmp(name, "page-length") == 0)
  {
  data->moldy.page_length=abs(atoi(value));
  return(0);
  }
if (g_ascii_strcasecmp(name, "xdr") == 0)
  {
  data->moldy.xdr=atoi(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "text-mode-save") == 0)
  {
  data->moldy.text_save=atoi(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "save-file") == 0)
  {
  data->moldy.save_file=g_strdup(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "restart-file") == 0)
  {
  data->moldy.restart_file=g_strdup(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "dump-file") == 0)
  {
  data->moldy.dump_file=g_strdup(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "backup-file") == 0)
  {
  data->moldy.backup_file=g_strdup(value);
  return(0);
  }
if (g_ascii_strcasecmp(name, "temp-file") == 0)
  {
  data->moldy.temp_file=g_strdup(value);
  return(0);
  }
/* Commands not yet incorporated */
if (g_ascii_strcasecmp(name, "molecular-cutoff") == 0)
  {
  return(0);
  }

/* Command not recognized */
return(1);
}

/* Determine energy units to use */
void calc_units(struct model_pak *data)
{
gdouble e_units;

e_units = data->moldy.mass_unit * SQR(data->moldy.length_unit) /
               SQR(data->moldy.time_unit) * AVOGADRO/1e3;

if (fabs( e_units - 1) < ENERGY_TOLERANCE)
  {
  data->moldy.energy_unit = KJMOL;
  return;
  }
if (fabs( e_units - 4.184) < ENERGY_TOLERANCE)
  {
  data->moldy.energy_unit = KCAL;
  return;
  }
if (fabs( e_units - 96.4845) < ENERGY_TOLERANCE)
  {
  data->moldy.energy_unit = EV;
  return;
  }
if (fabs( e_units - 1389.3544) < ENERGY_TOLERANCE)
  {
  data->moldy.energy_unit = E2A;
  return;
  }

data->moldy.energy_unit = OTHERS;

}

/* Select species at random */
struct moldyspec_pak *random_species(GSList *species, gint nmols)
{
gint   sel = rand() % nmols;
GSList *spec_list;
struct moldyspec_pak *specdata = NULL;

spec_list = species;
while( spec_list != NULL)
  {
  specdata = (struct moldyspec_pak*) spec_list->data;
  if (sel < specdata->num_mols)
    return(specdata);

  sel -= specdata->num_mols;
  spec_list = g_slist_next(spec_list); 
  }
return(specdata);
}

/* Generate initial configuration using skew start method */
void skew_start(struct model_pak *data, GSList *spec_list)
{
gint j, ispec, imol, iatom=0;   /* Counters for species, molecules etc */
GSList *temp_list, *slist;
gdouble mass = 0.0;             /* Whole system mass                  */
gdouble       n_third = pow((gdouble)data->moldy.num_mols,1.0/3.0);
gint          nz = 1, ny = (gint)(n_third+0.5), nx = (gint)(SQR(n_third)+0.5);
gint          nmols = data->moldy.num_mols;
gint          nspec = g_slist_length(data->moldy.species);
gint          mol_select[nspec];
gdouble       delta_x = (gdouble)nx / nmols,
              delta_y = (gdouble)ny / nmols,
              delta_z = (gdouble)nz / nmols;
gchar  *label;
gdouble cart[3], quat[4], p_f_sites[3];
struct moldyspec_pak *specdata;
struct site_pak *site = NULL;
struct core_pak *core;
struct atomtype_pak *type;

srand(time(NULL)+rand()); /* Re-seed random number generator */

data->periodic = 3;

/* Loop through species */
for (temp_list = spec_list; temp_list; temp_list = g_slist_next(temp_list))
  {
  specdata = (struct moldyspec_pak*) temp_list->data;
  ispec = g_slist_index(spec_list, specdata);
  mol_select[ispec] = specdata->num_mols;
  for (slist = specdata->sites; slist; slist=g_slist_next(slist))
    {
    site = (struct site_pak*) slist->data;
    type = (struct atomtype_pak*) site->type;
    mass += type->mass * specdata->num_mols;
    }
  ispec++;
  }

/* Calculate lattice parameters */
data->pbc[0] = data->pbc[1] = data->pbc[2]  /* L = cube root of mass/density */
             = pow(mass/data->moldy.density, 1.0/3.0);
data->pbc[3] = data->pbc[4] = data->pbc[5] = PI/2.0;
matrix_lattice_init(data);

for(imol=0; imol < nmols; imol++)
  {
  do
    {
    specdata = random_species(spec_list, nmols);
    ispec = g_slist_index(spec_list, specdata);
    }
    while( mol_select[ispec] <= 0);

  mol_select[ispec]--;

  if (specdata != NULL)
    {
    if (g_slist_length(specdata->sites) > 1)
      quat_random(quat);

    /* Identify and add atoms for all sites on species */
    for(slist = specdata->sites; slist; slist = g_slist_next(slist))
      {
      site = (struct site_pak*) slist->data;
      type = (struct atomtype_pak*) site->type;
      label = g_strdup((site->type)->label);
      core = core_new(label, NULL, data);
      cart[0] = imol*delta_x;
      cart[1] = imol*delta_y;
      cart[2] = imol*delta_z;
      for(j=3;j--;)
        {
        cart[j] -= floor(cart[j]);
        p_f_sites[j] = site->x[j];
        }
      if (g_slist_length(specdata->sites) > 1)
        quat_rotate(p_f_sites, quat);
      vecmat(data->latmat, cart); /* Conv frac coords to real coords */
      for(j=3;j--;)
        cart[j] += p_f_sites[j];
      ARR3SET(core->x, cart);
      core->charge = type->charge;
      core->lookup_charge = FALSE;
      core->atom_code = elem_symbol_test(g_strstrip(label));
      core->molecule = imol;
      data->cores = g_slist_prepend(data->cores, core);
      iatom++;
      g_free(label);
      }
    }
  }
}

/* Add Moldy site to GDis' element array */
#define DEBUG_ADD_SITE_ELEM 0
void add_site_elem(struct model_pak *data, GSList *sitetype_list)
{
gint code;
GSList *list;
struct elem_pak elem_data;
struct atomtype_pak *site;

list = sitetype_list;
while( list != NULL)
  {
  site = (struct atomtype_pak *) list->data;
  code = elem_symbol_test(site->label);
  get_elem_data(code, &elem_data, data);

  elem_data.charge = site->charge;
  elem_data.weight = site->mass;

  /* NB. Necessary fudge to allow for different site types of same element */
  elem_data.number = sysenv.num_elements+site->id;

  if (!code)
    elem_data.name = g_strdup(site->label);

  put_elem_data(&elem_data, data);

#if DEBUG_ADD_SITE_ELEM
printf(" *** Site elemental data ***\n");
printf("symbol: %s\n", elem_data.symbol);
printf("  name: %s\n", elem_data.name);
printf("number: %d\n", elem_data.number);
printf("weight: %f\n", elem_data.weight);
printf("  cova: %f\n", elem_data.cova);
printf("   vdw: %f\n", elem_data.vdw);
printf("charge: %f\n", elem_data.charge);
printf("colour: %-5f %-5f %-5f\n", elem_data.colour[0],
                                   elem_data.colour[1],
                                   elem_data.colour[2]);
#endif
  list = g_slist_next(list);
  }
}

/***************************************/
/* create new site type list entry     */
/***************************************/
#define DEBUG_CREATE_SITETYPE 0
struct atomtype_pak *create_sitetype(gchar *label, gint id, gdouble mass, gdouble charge)
{
struct atomtype_pak *type_data;

type_data = g_malloc(sizeof(struct atomtype_pak));
if (!type_data)
   return(NULL);

/* remove anything but alphabetic chars */
type_data->label = g_strdup(label);
type_data->id = id;
type_data->charge = charge; 
type_data->mass = mass;

#if DEBUG_CREATE_SITETYPE
printf("site type: id %d label %s mass %f charge %f\n", type_data->id,
   type_data->label, type_data->mass, type_data->charge);
#endif
return(type_data);
}

/***************************************/
/* create new site list entry          */
/***************************************/
#define DEBUG_NEW_SITE 0
struct site_pak *new_site(gdouble *x, struct atomtype_pak *type)
{
struct site_pak *site_data;

site_data = g_malloc(sizeof(struct site_pak));
if (!site_data)
   return(NULL);

ARR3SET(site_data->x, x); 
site_data->type = type; 
site_data->bonds = NULL; 

#if DEBUG_NEW_SITE
printf("site: id %d x %f y %f z %f\n", type->id,
   site_data->x[0], site_data->x[1], site_data->x[2]);
#endif
return(site_data);
}

/***************************************/
/* create new species list entry       */
/***************************************/
#define DEBUG_CREATE_SPECIES 0
struct moldyspec_pak *create_species(gchar *name,
     GSList *sites, gdouble *cofm, gint nmols, gboolean framework)
{
gint i;
struct moldyspec_pak *specdata;
GSList *list;       /* Site info list */
struct site_pak *site, *site_data;

specdata = g_malloc(sizeof(struct moldyspec_pak));
site = g_malloc(sizeof(struct site_pak));

if (!specdata)
   return(NULL);

specdata->name = g_strdup(name);
specdata->sites = NULL;
specdata->num_mols = nmols;
specdata->framework = framework;

/* Centre species */
for (list = sites; list; list=g_slist_next(list))
  {
  site_data = (struct site_pak *) list->data;
  if (!framework)
    for(i=3; i--;)
      site_data->x[i] -= cofm[i];
  site = new_site(site_data->x, site_data->type);
  if (site_data->bonds)
    site->bonds = site_data->bonds;
  specdata->sites = g_slist_append(specdata->sites, site);
  }

specdata->mols = NULL;

#if DEBUG_CREATE_SPECIES
printf("Name: %s\nSites: %d\n", specdata->name, g_slist_length(specdata->sites));
printf("No. of mols: %d\n", nmols);
#endif

return(specdata);
}

/***************************************/
/* create potential parm list entry    */
/***************************************/
#define DEBUG_CREATE_POT 0
struct pot_pak *create_pot(gchar *label1, gdouble charge1, gchar *label2, gdouble charge2)
{
gint i;
struct pot_pak *potdata;

potdata = g_malloc(sizeof(struct pot_pak));
if (!potdata)
   return(NULL);

potdata->site1_label = g_strdup(label1);
potdata->site1_charge = charge1;
potdata->site2_label = g_strdup(label2);
potdata->site2_charge = charge2;

for(i=0; i < MAX_POT_PARMS; i++)
  potdata->parm[i] = 0.0;

return(potdata);
}

/****************************/
/* moldy parsing routine    */
/****************************/
#define DEBUG_LOAD_MOLDY 0
gint read_moldy(gchar *filename, struct model_pak *data)
{
gint i, j, idi, idj;
gint flag=0, cflag=0;
gint num_species=0, num_mols=0;             /* Cumulative no of species and molecules */
gint num_tokens, num_sites=0;
gchar **buff, line[LINELEN];
gchar *fullname;
struct moldyspec_pak *specdata=NULL;
struct atomtype_pak *typedata=NULL;
struct atomtype_pak *tempdata;
struct site_pak *site, *sitedata;
FILE *fp;
gchar name[LINELEN], value[LINELEN];
gchar *specname=NULL;
gint id, iatom=0;
gint n_items;
gdouble mass=-1;                            /* Site mass */
gdouble charge1=0.0, charge2=0.0;           /* Site charges */
gchar *label1=NULL;                         /* Site element label 1 */
gchar *label2=NULL;                         /* Site element label 2 */
gdouble p_f_sites[3];                       /* Site coords in principal frame */
gdouble quat[4];			    /* Quaternions */
gdouble cofm[3];			    /* Species centre of mass */
gint nx=1, ny=1, nz=1;                      /* Cell repeat factors */
gdouble x[3];                               /* Molecule cofm coords */
gdouble total_mass=0;                       /* Total mass of species */
GSList *species_list=NULL, *sitetype_list=NULL, *site_list=NULL;
GSList *slist, *temp_list;
struct pot_pak *pot;			    /* Potential parameter storage */
struct core_pak *core;
gboolean framework = FALSE;

/* checks */
if (data == NULL)
  return(1);
if (filename == NULL)
  return(2);

if (strlen(filename) < FILELEN)
  strcpy(data->filename, filename);
else
  {
  gui_text_show(ERROR, "File name is too long.\n");
  return(-1);
  }

fp = fopen(filename, "rt");
if (!fp)
  return(1);

#if DEBUG_LOAD_MOLDY
printf("Creating new model: %s\n", filename);
#endif

data->moldy.num_mols = 0;

/* control parameter search */
while (!fgetline(fp, line))
  {
  g_strstrip(line);
  if (strlen(line))
    {
    num_tokens = sscanf(line, " %[^= ] = %127[^#]", name, value);
    if (num_tokens == 2 || !strncasecmp(name, "title", 5))
      {
      flag = match_key(data, name, value);
      cflag++;  /* File contains control info */
      }
    if (flag)
      printf("Unknown keyword \'%s\'\n",name);

    if (!g_ascii_strcasecmp(line, "end"))
      break;
    }
  }

/* Calculate energy units */
calc_units(data);

/* Use restart filename if read from restart */
if (!g_ascii_strncasecmp(g_path_get_basename(filename), TMPNAME, strlen(TMPNAME)))
  {
  data->id = MOLDY_RES;
  filename = g_strdup(data->moldy.restart_file);
  }
else
/* Read from system specification file given in control file */
  if (strlen(data->moldy.sysspec_file) > 0)
     {
     data->id = MOLDY;
     fclose(fp);

  /* Use current working directory unless alternate path already in filename */
     if (g_ascii_strcasecmp(data->moldy.sysspec_file, g_path_get_basename(data->moldy.sysspec_file)))
       fullname = g_strdup(data->moldy.sysspec_file);
     else
       fullname = g_build_filename(sysenv.cwd, data->moldy.sysspec_file, NULL);

     g_free(data->moldy.control_file);
     if (g_ascii_strncasecmp(g_path_get_basename(filename), TMPNAME, strlen(TMPNAME)))
       data->moldy.control_file = g_path_get_basename(filename);
     else
       data->moldy.control_file = g_strdup("");
     
     filename = g_strdup(data->moldy.sysspec_file);

     fp = fopen(fullname, "rt");
     if (!fp)
       {
       gui_text_show(ERROR, "Sys-spec file not found. ");
       return(1);
       }
     }

flag = 0; /* Reset flag for counting frameworks */

if (!cflag) /* File not headed by control info */
  rewind(fp);

/* Loop, parsing 'line' for name of new species. */
while (!fgetline(fp, line))
  {
  g_strstrip(line);
  buff = tokenize(line, &num_tokens);
  g_assert(buff != NULL);

  if (num_tokens < 4)
    {
    if (num_mols > 0) /* previous species complete */
      {
      specdata = create_species(specname, site_list, cofm, num_mols, framework);
      species_list = g_slist_append(species_list, specdata);
      data->moldy.num_mols += num_mols;
      if (total_mass > 0.0)
        {
        for(i = 3; i--;)
          cofm[i] /= total_mass;
#if DEBUG_LOAD_MOLDY
printf("----------------------------------------\n");
printf("Species: %s\n", specname);
printf("Num of mols: %d\n", num_mols);
printf("Total mass: %f\n", total_mass);
printf("C of M: %f %f %f\n", cofm[0], cofm[1], cofm[2]);
printf("----------------------------------------\n");
#endif
        }
      if (site_list) free_slist(site_list);
      }
    if (g_ascii_strcasecmp(*buff, "end") == 0)
      break;
    specname = g_strdup(*buff);
    num_mols = (gint)g_ascii_strtod(*(buff+1), NULL);
    num_species++;
    num_sites = 0;
    total_mass = 0.0;

    site_list = NULL;

    VEC3SET(cofm, 0.0, 0.0, 0.0);
    framework = FALSE;
    if (num_tokens == 3)
      {
      if (g_ascii_strcasecmp(*(buff+2), "framework") == 0)
         {
         framework = TRUE;
         flag++;
         if (flag > 1)
            gui_text_show(WARNING, "More than one framework species found.\n");
         }
      }
    }
  switch (num_tokens)
    {
    case 7:
      label1 = g_strdup(*(buff+6));
    case 6:
      charge1 = g_ascii_strtod(*(buff+5), NULL);
      charge1 *= (data->moldy.charge_unit/ELCHARGE);
    case 5:
      mass = g_ascii_strtod(*(buff+4), NULL);
      mass *= (data->moldy.mass_unit/AMU);
    case 4:
      id = (gint)g_ascii_strtod(*buff, NULL);
      p_f_sites[0] = g_ascii_strtod(*(buff+1), NULL);
      p_f_sites[1] = g_ascii_strtod(*(buff+2), NULL);
      p_f_sites[2] = g_ascii_strtod(*(buff+3), NULL);
      site = g_malloc(sizeof(struct site_pak));  /* Make new list element  */
      ARR3SET(site->x, p_f_sites);
      num_sites++;

      flag = 0;
      temp_list = sitetype_list;
      while(temp_list != NULL)
        {
        typedata = (struct atomtype_pak*) temp_list->data;
        if (typedata->id == id)
           {
           /* Error checking */
           if (num_tokens > 4 && mass != typedata->mass)
             {
             printf("Error: Incorrect mass for site %d\n", id);
             return(1);
             }
           if (num_tokens > 5 && charge1 != typedata->charge)
             {
             printf("Error: Incorrect charge for site %d\n", id);
             return(1);
             }
           if (num_tokens > 6 &&  g_ascii_strcasecmp(label1, typedata->label))
             {
             printf("Error: Incorrect label for site %d\n", id);
             return(1);
             }
           /* End of error checking */
           mass = typedata->mass; charge1 = typedata->charge;
           flag++;
           break;
           }
        temp_list = g_slist_next(temp_list);
        }
      if (!flag)
        {
        if (mass == -1)
          {
          printf("Insufficient data for site %d\n", id);
          break;
          }
        else
          {
          typedata = create_sitetype(label1, abs(id), mass, charge1);
          sitetype_list = g_slist_append(sitetype_list, typedata);
          }
        }
      for(i = 3; i--;)
        cofm[i] += mass*site->x[i];
      total_mass += mass;
      site->type = typedata;
      site_list = g_slist_append(site_list, site);
      break;
    }
  g_strfreev(buff);
  };

if (data->moldy.num_mols == 0 || num_species == 0)
  {
  printf("Error: Incorrect format for sys-spec file\n");
  return(1);
  }
data->moldy.species = species_list;
data->atom_types = sitetype_list;
add_site_elem(data, sitetype_list);

/* potential parm specification */
if (fgetline(fp, line))
  return(-1);

buff = tokenize(line, &num_tokens);

for(i = 0; potspec[i].name; i++)             /* Is 'name' a known type? */
  if (g_ascii_strcasecmp(*buff, potspec[i].name) == 0)
     break;

n_items = 0;
if (!potspec[i].name)                        /* Did the loop find 'name'? */
  {
  printf("Unrecognized potential function in %s \n", filename);
  }
else
  {
  data->moldy.pot_type = i;                  /* yes */ 
  n_items = potspec[i].npar;
#if DEBUG_LOAD_MOLDY
  printf("Known potential type (%d) found\n", i);
#endif
  }

do
  {
  if (fgetline(fp, line))
     break;

  buff = tokenize(line, &num_tokens);

  if (num_tokens > 2 && n_items)
    {
    idi = atoi(*buff);
    idj = atoi(*(buff+1));
    flag = 0;
    temp_list = sitetype_list;
    while( temp_list != NULL)
      {
      tempdata = (struct atomtype_pak*) temp_list->data;
      if (tempdata->id == idi)
         {
         label1 = g_strdup(tempdata->label);
         charge1 = tempdata->charge;
         flag++;
         }
      if (tempdata->id == idj)
         {
         label2 = g_strdup(tempdata->label);
         charge2 = tempdata->charge;
         flag++;
         }
      temp_list = g_slist_next(temp_list);
      }
    if (flag == 2)
      {
      pot = create_pot(label1, charge1, label2, charge2);

      for(i=0; i < n_items; i++)
        {
        if (num_tokens > 2+i)
          pot->parm[i] = g_ascii_strtod(*(buff+2+i), NULL);
        else
          {
          printf("Insufficient potential parameters\n");
          return(-1);
          }
        }

      data->moldy.pots = g_slist_append(data->moldy.pots, pot);
      }
#if DEBUG_LOAD_MOLDY
  printf("Found pot parms for atom sites %d and %d\n", idi, idj);
#endif
    }

  } while (g_ascii_strcasecmp(*buff, "end"));

/* next line should be lattice parameters if 3D-periodic */
if (!fgetline(fp, line))
  {
  buff = tokenize(line, &num_tokens);
  if (num_tokens > 5)
    {
    data->periodic = 3;
    data->pbc[0] = g_ascii_strtod(*buff, NULL);
    data->pbc[0] *= (data->moldy.length_unit/ANGST);
    data->pbc[1] = g_ascii_strtod(*(buff+1), NULL);
    data->pbc[1] *= (data->moldy.length_unit/ANGST);
    data->pbc[2] = g_ascii_strtod(*(buff+2), NULL);
    data->pbc[2] *= (data->moldy.length_unit/ANGST);
    data->pbc[3] = D2R * g_ascii_strtod(*(buff+3), NULL);
    data->pbc[4] = D2R * g_ascii_strtod(*(buff+4), NULL);
    data->pbc[5] = D2R * g_ascii_strtod(*(buff+5), NULL);
    nx = (gint) g_ascii_strtod(*(buff+6), NULL);
    ny = (gint) g_ascii_strtod(*(buff+7), NULL);
    nz = (gint) g_ascii_strtod(*(buff+8), NULL);
    matrix_lattice_init(data);
    }

/* Check for sensible values */
  if (nx < 1 || ny < 1 || nz < 1)
    return(1);

  for (i=6; i--;)
    if (data->pbc[i] < 0.0)
      return(1);

  if (nx * ny * nz > 1000)
    {
    gui_text_show(WARNING, "Too many image cells.\n");
    return(1);
    }

  num_mols=0;
/* Read cofms of each molecule */
  do
    {
    if (fgetline(fp, line))
      break;
    n_items = sscanf(line, "%s %lf %lf %lf %lf %lf %lf %lf",
          name, &x[0], &x[1], &x[2], &quat[0], &quat[1], &quat[2], &quat[3]);
    g_strstrip(name);
    if (g_ascii_strcasecmp(name, "end") == 0)
      break;

    vecmat(data->latmat, x); /* Conv frac cofms to real coords */

/* Loop to identify species */
    temp_list = species_list;
    while( temp_list != NULL)
      {
      specdata = (struct moldyspec_pak*) temp_list->data;
      if (!g_ascii_strcasecmp(specdata->name, name))
         break;
      temp_list = g_slist_next(temp_list);
      }
/* Unknown species */
    if (temp_list == NULL)
      {
      gui_text_show(WARNING, "Unknown species ignored.\n");
      }
    else
      {
      num_mols++;
      if (g_slist_length(specdata->sites) > 1 && !specdata->framework)
        {
        if (n_items < 8)
           printf("Too few quaternions for species \"%s\" - 4 needed\n", name);
        else if (fabs(1.0 - SQR(quat[0]) - SQR(quat[1])
             - SQR(quat[2]) - SQR(quat[3])) > 1e-4)
             printf("Incorrect quaternion values found for molecule %d\n", num_mols);
        }

/* Loop through species sites, matching sites and adding atoms */
      for (slist=specdata->sites; slist; slist = g_slist_next(slist))
        {
        sitedata = (struct site_pak *) slist->data;
        typedata = sitedata->type;
        label1 = g_strdup(typedata->label);
        core = core_new(label1, NULL, data);
        ARR3SET(p_f_sites, sitedata->x);
/* coords */
        if (n_items == 8 && g_slist_length(specdata->sites) > 1 && !specdata->framework)
          quat_rotate(p_f_sites, quat); /* Rotate coords using quaternions */

        /* Add rotated atom positions to cofms */
        for(j=3; j--;)
          core->x[j] = x[j] + p_f_sites[j];
        core->charge = typedata->charge;
        core->lookup_charge = FALSE;
        core->mass = typedata->mass;
        core->lookup_mass = FALSE;
        core->atom_code = elem_symbol_test(label1);
        core->molecule = num_mols;
        printf("Hey %d\n", core->molecule);
        data->cores = g_slist_prepend(data->cores, core);
        iatom++;
        }
      }
    } while (g_ascii_strcasecmp(*buff, "end"));
    if (data->moldy.num_mols != nx*ny*nz*num_mols)
      gui_text_show(WARNING, "Wrong number of molecules in sys-spec file.\n");
  }
else /* Non-periodic system */
  {
#if DEBUG_LOAD_MOLDY
printf("Loading non-periodic system\n");
#endif
  if (num_species == 1)
    {
    data->periodic = 0;
    /* Loop through species */
    temp_list = species_list;
    while( temp_list != NULL)
      {
      specdata = (struct moldyspec_pak*) temp_list->data;

      /* Identify and add atoms for all sites on species */
      for (slist=specdata->sites; slist; slist=g_slist_next(slist))
        {
        sitedata = slist->data;
        typedata = sitedata->type;
        label1 = g_strdup(typedata->label);
        core = core_new(label1, NULL, data);
        ARR3SET(core->x, sitedata->x);
        core->charge = typedata->charge;
        core->lookup_charge = FALSE;
        core->atom_code = elem_symbol_test(label1);
        core->molecule = 1;
        data->cores = g_slist_prepend(data->cores, core);
        }
      temp_list = g_slist_next(temp_list);
      }
    }
  else
    {
    skew_start(data, species_list);
    }
  }

data->num_asym = g_slist_length(data->cores);

/* Set up filenames */
if (data->moldy.control_file == NULL || !g_ascii_strcasecmp(data->moldy.control_file, ""))
  {
  g_free(data->moldy.control_file);
  data->moldy.control_file = g_strdup_printf("%s.min", parse_strip(filename));
  }
else if (data->id == MOLDY_RES)
  if (strlen(data->moldy.restart_file) < FILELEN)
    strcpy(data->filename, data->moldy.restart_file);
  else
    {
    gui_text_show(ERROR, "File name is too long.\n");
    return(-1);
    }
else
  if (strlen(data->moldy.sysspec_file) < FILELEN)
    strcpy(data->filename, data->moldy.sysspec_file);
  else
    {
    gui_text_show(ERROR, "File name is too long.\n");
    return(-1);
    }

g_free(data->basename);
data->basename = g_strdup(parse_strip(filename));
g_free(data->moldy.sysspec_file);
data->moldy.sysspec_file = g_strdup_printf("%s.min", g_strstrip(data->basename));

/* Set periodic repeat cells */
data->image_limit[1] = nx;
data->image_limit[3] = ny;
data->image_limit[5] = nz;

/* summary */
#if DEBUG_LOAD_MOLDY
printf("-------------------------------------\n");
printf("         structure found : %s\n", data->basename);
printf("               time step : %f\n", data->moldy.deltat);
printf("        no of time steps : %d\n", data->moldy.nsteps);
printf("             temperature : %f\n", data->moldy.temperature);
printf("        potentials found : %d\n", g_slist_length(data->moldy.pots));
printf("      species data found : %d\n", g_slist_length(data->moldy.species));
printf("         molecules found : %d\n", data->moldy.num_mols);
printf("   total number of atoms : %d\n", g_slist_length(data->cores));
printf("-------------------------------------\n");
#endif

/* display init */
model_prep(data);

if (nx > 1 || ny > 1 || nz > 1)
  {
  space_make_images(CREATE, data);
  coords_init(CENT_COORDS, data);
  }

/* clean up & exit */
fclose(fp);
return(0);
}

/***************************************/
/* find all unique atom types in model */
/***************************************/
#define DEBUG_FIND_SITES 0
GSList *find_unique_sites(struct model_pak *data)
{
gint j, found;
gint num_unique_sites=0, num_atoms =  g_slist_length(data->cores);
gdouble q1, q2;  /* Total atom charges */
GSList *list, *types;
struct core_pak *core1, *core2;
struct atomtype_pak *site_type;
gint *atom_type;

#if DEBUG_FIND_SITES
  printf("%d atoms to search.\n", num_atoms);
#endif

/* destroy old list */
if (data->atom_types)
  atomtype_data_free(data->atom_types);

if (g_slist_length(data->cores) > 0)
  atom_type = g_malloc(num_atoms*sizeof(gint));
else
  return 0;

types = NULL;

/* go through all atoms */
for (list=data->cores; list; list=g_slist_next(list))
  {
  found=0;
  core1 = (struct core_pak *) list->data;

  if (core1->status & (DELETED | HIDDEN))
    continue;
  q1 = atom_charge(core1);

/* determine if current atom type is in the list */
  for (j=0; j < num_unique_sites; j++)
    {
    core2 = (struct core_pak *) g_slist_nth_data(data->cores, atom_type[j]);
    q2 = atom_charge(core2);
    if (fabs(q1-q2) < FRACTION_TOLERANCE &&
       (!g_ascii_strcasecmp(core1->atom_label, core2->atom_label)))
          {
          found++;
          break;
          }
    }
/* not in list, so add it */
  if (!found)
    {
    atom_type[j] = g_slist_index(data->cores, core1);
    num_unique_sites++;
    }
#if DEBUG_FIND_SITES
  printf("Examining atom %d - Total unique sites %d\n",
          g_slist_index(data->cores, core1)+1, num_unique_sites);
#endif
  }

for (j=0; j < num_unique_sites; j++)
  {
/* Copy data to site arrays */
  core1 = (struct core_pak *) g_slist_nth_data(data->cores, atom_type[j]);

  q1 = atom_charge(core1);
  site_type = create_sitetype(core1->atom_label, j+1, atom_mass(core1), q1);
#if DEBUG_FIND_SITES
  printf("Found: ");
  printf("Site %d Element %s Weight %f Charge %f\n", j+1, site_type->label,
    site_type->mass, site_type->charge);
#endif
  types = g_slist_prepend(types, site_type);
  }

#if DEBUG_FIND_SITES
  printf("Identified %d unique atom types\n", g_slist_length(types));
#endif

types = g_slist_reverse(types);

/* Sort sites into chemical order */
/*
types = g_slist_sort(types, (gpointer) compare_sites);
*/
g_free(atom_type);

return(types);
}

/***************************************/
/* find all unique elements in a model */
/***************************************/
#define DEBUG_FIND_SITES2 0
void find_unique_sites2(struct model_pak *data)
{
gint j, found;
gint num_unique_sites=0, num_atoms = g_slist_length(data->cores);
gdouble q1, q2;  /* Total atom charges */
GSList *list;
struct core_pak *core1, *core2;
struct atomtype_pak *site_type;
gint *atom_type;

#if DEBUG_FIND_SITES
  printf("%d atoms to search.\n", g_slist_length(data->cores));
#endif

/* destroy old list */
atomtype_data_free(data->atom_types);
data->atom_types = NULL;
atom_type = g_malloc(num_atoms*sizeof(gint));

/* go through all atoms */
for (list = data->cores; list; list=g_slist_next(list))
  {
  found=0;
  core1 = (struct core_pak *) list->data;

  if (core1->status & (DELETED | HIDDEN))
    continue;
  q1 = atom_charge(core1);
/* determine if current atom type is in the list */
  for (j=0; j < num_unique_sites; j++)
    {
    core2 = (struct core_pak *) g_slist_nth_data(data->cores, atom_type[j]);
    q2 = atom_charge(core2);
    if (fabs(q1-q2) < FRACTION_TOLERANCE &&
       (!g_ascii_strcasecmp(core1->atom_label, core2->atom_label)))
          {
          found++;
          break;
          }
    }
/* not in list, so add it */
  if (!found)
    {
    atom_type[j]=g_slist_index(data->cores, core1);
    num_unique_sites++;
    }
#if DEBUG_FIND_SITES
  printf("Examining atom %d - Total unique sites %d\n", g_slist_index(data->cores, core1)+1, num_unique_sites);
#endif
  }

for (j=0; j < num_unique_sites; j++)
  {
/* Copy data to site arrays */
  core1 = (struct core_pak *) g_slist_nth_data(data->cores, atom_type[j]);

  q1 = atom_charge(core1);
  site_type = create_sitetype(core1->atom_label, j+1, atom_mass(core1), q1);
#if DEBUG_FIND_SITES
  printf("Found: ");
  printf("Site %d Element %s Weight %f Charge %f\n", j+1, site_type->label,
    site_type->mass, site_type->charge);
#endif
  data->atom_types = g_slist_append(data->atom_types, site_type); 
  }

#if DEBUG_FIND_SITES
  printf("Identified %d unique atom types\n", num_unique_sites);
#endif

/* Sort sites into chemical order */
/*
data->atom_types = g_slist_sort(data->atom_types, (gpointer) compare_sites);
*/
g_free(atom_type);
}

/**********************************************/
/* test to see if two molecules are identical */
/**********************************************/
#define DEBUG_MATCH_MOLECULES 0
gint match_molecules(struct model_pak *data, 
         struct mol_pak *mol1, struct mol_pak *mol2, gdouble tolerance)
{
gint i, j, k, l;
gint n, num_ignore=0;
gint total_atoms = count_visible_cores(mol1->cores);
gint match=0; /* 0 = no match, 1 = match */
gint num_sites;
gint *pair; /* List of sites in molecule 2 matching those in molecule 1 */
gint nrot; /* No of rotations when attempting to match molecules */
struct core_pak *core1, *core2, *temp_core;
GSList *list1, *list2, *temp_list;
gdouble *p_f_sites1[3], *p_f_sites2[3]; /* Principle frame sites */
gdouble rot[9]; /* rotation matrix of second molecule */
gdouble vec[3];
gdouble charge1, charge2;

/* Molecule being matched */
num_sites = count_visible_cores(mol1->cores);
num_ignore = total_atoms - num_sites;

pair = g_malloc(total_atoms*sizeof(gint)); /* Record matching sites */

p_f_sites1[0] = (gdouble *) g_malloc(num_sites*sizeof(gdouble));
p_f_sites1[1] = (gdouble *) g_malloc(num_sites*sizeof(gdouble));
p_f_sites1[2] = (gdouble *) g_malloc(num_sites*sizeof(gdouble));
p_f_sites2[0] = (gdouble *) g_malloc(num_sites*sizeof(gdouble));
p_f_sites2[1] = (gdouble *) g_malloc(num_sites*sizeof(gdouble));
p_f_sites2[2] = (gdouble *) g_malloc(num_sites*sizeof(gdouble));

for(i = 0; i < 3; i++)
  for(j = 0; j < num_sites; j++)
    p_f_sites1[i][j] = p_f_sites2[i][j] = 0.0;

/* Molecule being matched */
for (list1=mol1->cores, j=0; list1; list1=g_slist_next(list1))
  {
  core1 = (struct core_pak *) list1->data;
  if (core1->status & (DELETED | HIDDEN))
    continue;
  for(i = 3; i--;)
    p_f_sites1[i][j] = core1->x[i];
  j++;
  }

/* Molecule being compared with */
for (list2 = mol2->cores, j=0; list2; list2=g_slist_next(list2))
  {
  core2 = (struct core_pak *) list2->data;
  if (core2->status & (DELETED | HIDDEN))
    continue;
  for(i = 3; i--;)
    p_f_sites2[i][j] = core2->x[i];
  j++;
  }

#if DEBUG_MATCH_MOLECULES
  printf("Molecule 1: Total atoms          %d \n", total_atoms);
  printf("            Hidden/deleted atoms %d \n", num_ignore);
#endif

calc_pfc(data, mol1, p_f_sites1, num_sites);
calc_pfc(data, mol2, p_f_sites2, num_sites);

/* Number of rotations for monatomic and polyatomic species */
if (num_sites == 1 || tolerance == 0.0)
  {
  nrot = 1;
  tolerance = FLOATING_POINT;
  }
else
  nrot = 6;

match = 1;
for(i=total_atoms;i--;) /* Set flags */
  pair[i] = -1;
i = 0;
                             /* j=index for visible atoms in mol 2 */
k = total_atoms-num_ignore;  /* k=index for hidden/deleted atoms in mol 1*/
                             /* l=index for all atoms in mol 1 */
for (list1 =  mol1->cores, l=0; list1; list1 = g_slist_next(list1), l++)
  {
  core1 = (struct core_pak *) list1->data;
  if (core1->status & (DELETED | HIDDEN))
    {
    pair[k] = l;
    k++;
    continue;
    }
  charge1 = atom_charge(core1);
  for (list2 = mol2->cores, j=0; list2; list2 = g_slist_next(list2))
    {
    core2 = (struct core_pak *) list2->data;
    if (core2->status & (DELETED | HIDDEN))
      continue;
    /* Do atom sites match? */
    charge2 = atom_charge(core2);

#if DEBUG_MATCH_MOLECULES
printf("Comparing %s %f and %s %f: Tolerance: %e\n",
      core1->atom_label, charge1, core2->atom_label, charge2, tolerance);
#endif
    if (!g_ascii_strcasecmp(core1->atom_label, core2->atom_label) &&
        fabs(charge1 - charge2) < FRACTION_TOLERANCE)
      {
      if (((fabs(p_f_sites2[0][j] - p_f_sites1[0][i]) < tolerance) &&
           (fabs(p_f_sites2[1][j] - p_f_sites1[1][i]) < tolerance) &&
           (fabs(p_f_sites2[2][j] - p_f_sites1[2][i]) < tolerance)))
        {
        pair[j] = l;
        break;
        }
      }
    j++;
    }
  i++;
  }

for(i = 0; i < num_sites; i++) /* All visible atoms paired? */
  if (pair[i] == -1)
    match = 0;

if (match && tolerance > FLOATING_POINT)
  {
  for(n = 0; n < nrot; n++)
    {
    match = 1;
    temp_list = NULL;
    for (i = 0; i < num_sites; i++)
      {
      /* Swap core list */
      temp_core = g_slist_nth_data(mol1->cores, pair[i]);
      temp_list = g_slist_append(temp_list, temp_core);
      }
      mol1->cores = temp_list;
    break;
    }

  VEC3SET(&rot[0], 1.0, 0.0, 0.0);
  VEC3SET(&rot[3], 0.0, 1.0, 0.0);
  VEC3SET(&rot[6], 0.0, 0.0, 1.0);

  switch(n)
    {
    case 0:  /* Rotate about x axis */
      rot[4] *= -1.0;
      rot[8] *= -1.0;
    break;
    case 1:  /* Rotate about y axis */
      rot[0] *= -1.0;
      rot[4] *= -1.0;
    break;
    case 2:  /* Rotate about x axis */
      rot[0] *= -1.0;
      rot[4] *= -1.0;
    break;
    case 3:  /* Rotate about z axis */
      rot[0] *= -1.0;
      rot[8] *= -1.0;
    break;
    case 4:  /* Rotate about x axis */
      rot[0] *= -1.0;
      rot[8] *= -1.0;
    break;
    case 5:
    break;
    }
/* apply rotation */
    for(j = 0; j < num_sites; j++)
      {
      vec[0] = p_f_sites2[0][j];
      vec[1] = p_f_sites2[1][j];
      vec[2] = p_f_sites2[2][j];
      vecmat(rot,vec);
      p_f_sites2[0][j] = vec[0];
      p_f_sites2[1][j] = vec[1];
      p_f_sites2[2][j] = vec[2];
      }

#if DEBUG_MATCH_MOLECULES
  printf("Match %d on rotation %d\n", match, n);
#endif
  }

return (match);
}

#define DEBUG_CALC_SPECIES 0
gint calc_species(struct model_pak *data, gdouble tolerance)
{
gint i, j, nbond;
gint flag, match, nv;
GSList *list, *slist, *mlist, *clist;
GSList *site_list, *spec_list=NULL;
struct mol_pak *mol1, *mol2;
struct elem_pak elem_data;
struct moldyspec_pak *species=NULL, *specdata=NULL;
struct site_pak *site = NULL;
struct atomtype_pak *sitetype = NULL;
struct core_pak *core, *neighbour;
gdouble *p_f_sites[3];
/* No of hidden atoms in each molecule and counter for visible atoms */
gint *num_atoms;
gchar *name = NULL;
gdouble cofm[3], charge;
gint label_no;

data->atom_types = find_unique_sites(data);
num_atoms = (gint *)g_malloc(g_slist_length(data->moles)*sizeof(gint));

#if DEBUG_CALC_SPECIES
  printf("No of molecules: %d\n", g_slist_length(data->moles));
#endif

for (mlist=data->moles; mlist; mlist=g_slist_next(mlist))
  {
  mol1 = (struct mol_pak *) mlist->data;

  site_list = NULL;  /* Reset site list for this species */
  /* Place atoms in same order */
  if (g_slist_length(mol1->cores) > 1)
    mol1->cores = g_slist_sort(mol1->cores, (gpointer) internal_atom_sort);

  i = g_slist_index(data->moles, mol1);
  num_atoms[i] = count_visible_cores(mol1->cores);

#if DEBUG_CALC_SPECIES
  printf("%d atom%sin molecule %d\n", num_atoms[i], (num_atoms[i]!=1?"s ":" "), i+1);
#endif
  match = 0;
  if (num_atoms[i] > 0)
    {
    for (list = spec_list; list; list = g_slist_next(list))
      {
      specdata = (struct moldyspec_pak *) list->data;
      mol2 = (struct mol_pak *) ((specdata->mols)->data);
      if (num_atoms[i] == num_atoms[g_slist_index(data->moles, mol2)])
        match = match_molecules(data, mol1, mol2, tolerance);

      if ( match)
        break;
      }

    if ( match)
      {
      specdata->mols = g_slist_append(specdata->mols, mol1);
      specdata->num_mols++;
#if DEBUG_CALC_SPECIES
  printf("Incrementing species %d\n", g_slist_index(spec_list, specdata)+1);
#endif
      }
    else /* Add new species */
      {
#if DEBUG_CALC_SPECIES
  printf("Assigning new species %d\n", g_slist_length(spec_list)+1);
#endif
      g_free(name);
      name = g_strdup_printf("Species%d", g_slist_length(spec_list)+1);

      VEC3SET(cofm, 0.0, 0.0, 0.0);

      for (list = mol1->cores; list; list = g_slist_next(list))
        {
        core = (struct core_pak *) list->data;
        if (core->status & (DELETED | HIDDEN))
          continue;
        get_elem_data(core->atom_code, &elem_data, data);
        charge = atom_charge(core);

        /* Use numerals in atom_label in species name if monatomic */
        label_no = -1; /* Default is no trailing numeral */
        if (num_atoms[i] == 1)
           {
           g_free(name);
           name = g_strdup(elem_data.name);
           for (j = 0; j < strlen(core->atom_label); j++)
             if (g_ascii_isdigit(core->atom_label[j]))
                {
                label_no = (gint)g_ascii_strtod(core->atom_label+j, NULL);
                break;
                }

           if (label_no >= 0)
             {
             g_free(name);
             name = g_strdup_printf("%s%d", elem_data.name, label_no);
             }

           /* Make sure name is unique */
           j=0;
           for (list = spec_list; list; list = g_slist_next(list)) 
             {
             specdata = (struct moldyspec_pak *) list->data;
             if (g_ascii_strcasecmp(name, specdata->name) == 0)
               {
               j++;
               if (j == 1)
                 specdata->name = g_strdup_printf("%s%d", name, j);
               }
             }
           
           if (j > 0)
             {
             g_free(name);
             name = g_strdup_printf("%s%d", elem_data.name, j+1);
             }
           }

        flag = 0;
        for (slist = data->atom_types; slist; slist = g_slist_next(slist))
          {
          sitetype = (struct atomtype_pak *) slist->data;
          if (!g_ascii_strcasecmp(core->atom_label, sitetype->label) &&
               fabs(sitetype->charge - charge) < FRACTION_TOLERANCE)
            {
            flag++;
            break;
            }
          }

        if (!flag)  /* Shouldn't happen */
          printf("Error: Site %d (%s) in molecule %d not identified\n",
                  g_slist_index(mol1->cores, core), core->atom_label, i); 

        site = new_site(core->x, sitetype);

#if DEBUG_CALC_SPECIES
  printf("Site %d Label %s Mass %f Charge %f Bonds %d\n",
     g_slist_index(data->atom_types, sitetype)+1, sitetype->label, sitetype->mass,
     sitetype->charge, g_slist_length(core->bonds));
#endif

        /* Create list of bonds within species */
        for (clist = connect_neighbours(core); clist; clist = clist->next)
          {
          neighbour = (struct core_pak *) clist->data;
           
          nbond = g_slist_index(mol1->cores, neighbour)+1;
          site->bonds = g_slist_prepend(site->bonds, GINT_TO_POINTER(nbond));
          }
        site->bonds = g_slist_reverse(site->bonds);
        site_list = g_slist_append(site_list, site);
        }

      if (!data->moldy.set_framework)
        {
        nv = g_slist_length(site_list);
        p_f_sites[0] = (gdouble *)g_malloc(nv*sizeof(gdouble));
        p_f_sites[1] = (gdouble *)g_malloc(nv*sizeof(gdouble));
        p_f_sites[2] = (gdouble *)g_malloc(nv*sizeof(gdouble));

        /* Convert site positions to principal frame coords */
        for (slist = site_list; slist; slist = g_slist_next(slist))
          {
          nv = g_slist_index(site_list, slist->data);
          site = (struct site_pak *) slist->data;
          for(j = 3; j--;)
            p_f_sites[j][nv] = site->x[j];
          }
        calc_pfc(data, mol1, p_f_sites, g_slist_length(site_list));
        for (slist = site_list; slist; slist = g_slist_next(slist))
          {
          nv = g_slist_index(site_list, slist->data);
          site = (struct site_pak *) slist->data;
          for(j = 3; j--;)
            site->x[j] = p_f_sites[j][nv];
          }
        }
      else
        for (slist = site_list; slist; slist = g_slist_next(slist))
          {
          site = (struct site_pak *) slist->data;
          vectmat(data->latmat, site->x);
          }

      species = create_species(name, site_list, cofm, 1, FALSE);
      species->mols = g_slist_append(species->mols, mol1);
      spec_list = g_slist_append(spec_list, species);
      free_slist(site_list);
      }
    }
  }
#if DEBUG_CALC_SPECIES
for (list = spec_list; list; list = g_slist_next(list)) 
  {
  specdata = (struct moldyspec_pak *) list->data;
  i = g_slist_length(specdata->mols);
  if (i == 1)
    printf("There is %d molecule of %s\n", i, specdata->name);
  else
    printf("There are %d molecules of %s\n", i, specdata->name);
  }
#endif

if (g_slist_length(spec_list) > 1)
  spec_list = g_slist_sort(spec_list, (gpointer) compare_species);

data->moldy.species = spec_list;

g_free(name);
return(g_slist_length(spec_list));
}

#define DEBUG_DELETE_POTS 0
gint delete_moldy_potentials(struct model_pak *model)
{

if (model == NULL)
  model = sysenv.active_model;

if (model == NULL)
  return(1);

#if DEBUG_DELETE_POTS
  printf("Deleting Moldy potential parameters.\n");
#endif

free_slist(model->moldy.pots);
g_free(model->moldy.lib_file);

model->moldy.pots = NULL;
model->moldy.lib_file = g_strdup("");
model->moldy.pot_type = -1;

update_pot_dialog(model);

return(0);
}

#define DEBUG_READ_RESTART 0
gint read_moldy_restart(gchar *filename, struct model_pak *data)
{
gint status, random;
gchar *fullname;
GString *txt, * inname;
const gchar *ctemp;

/* checks */
if (data == NULL)
  data = sysenv.active_model;

if (g_find_program_in_path("ransub") == NULL)
  {
  gui_text_show(ERROR, "You need to have Moldy's ransub utility installed for this function.\n");
  return(2);
  }

if (data == NULL)
  return(1);
if (filename == NULL)
  return(2);

if (!g_ascii_strcasecmp(filename, g_path_get_basename(filename)))
  fullname = g_build_filename(sysenv.cwd, filename, NULL);
else
  fullname = g_strdup(filename);

if (fullname == NULL)
  return(3);

if (access(fullname, R_OK))
   gui_text_show(ERROR, "Restart file not found.\n");

if (!ascii_check(fullname))
   gui_text_show(ERROR, "Restart file not in correct format.\n");

if (data->moldy.restart_file)
  g_free(data->moldy.restart_file);
data->moldy.restart_file = g_path_get_basename(filename);

if (data->moldy.restart_dir)
  g_free(data->moldy.restart_dir);
data->moldy.restart_dir = g_strdup(sysenv.cwd);

if (data->moldy.sysspec_file)
  g_free(data->moldy.sysspec_file);
data->moldy.sysspec_file = g_path_get_basename(filename);

sysenv.file_type = DATA;
dialog_destroy_type(FILE_SELECT);

srand(time(NULL)+rand()); /* Re-seed random number generator */
txt = g_string_new(NULL);
inname = g_string_new(NULL);
ctemp = g_get_home_dir();
random = (gint)(rand() % 1000);
g_string_append_printf(inname, "%s/%s%d", ctemp, TMPNAME, random); /* Attach random number */
g_string_printf(txt, "ransub -r %s -o %s >&1 /dev/null",
                      fullname, inname->str);

#if DEBUG_READ_RESTART
printf("executing: [%s]\n", txt->str);
#endif

status = system(txt->str);

#if DEBUG_READ_RESTART
printf("reading: [%s]\n", inname->str);
#endif

if (status != -1)
  {
  /* pass to moldy load routine */
  status = read_moldy(inname->str, data);

  /* delete temporary file */
  g_string_printf(txt, "rm %s", inname->str);
  system(txt->str);
  }

/* free & exit */
g_string_free(txt, TRUE);
g_string_free(inname, TRUE);
g_free(fullname);

return(status);
}
#define DEBUG_READ_POTS 0
gint read_moldy_potentials(gchar *filename, struct model_pak *data)
{
gint i=0, pflag=0;
gint n_items;
gchar line[LINELEN], pline[LINELEN];
gint line_no=0;
FILE *fp;
gint num_tokens;
gchar **buff=NULL;
gchar *fullname;
gdouble *p;
struct pot_pak *pot=NULL;

/* Atom identifiers */
gchar atom1[4], atom2[4];
gdouble charge1, charge2;

/* checks */
if (data == NULL)
  {
  data = sysenv.active_model;
  g_free(data->moldy.lib_file);
  data->moldy.lib_file = g_strdup(filename);
  }

if (data == NULL)
  return(1);
if (filename == NULL)
  return(2);

/* Remember potential file directory */
g_free(data->moldy.lib_dir);
data->moldy.lib_dir = g_strdup(sysenv.cwd);

fullname = g_build_filename(sysenv.cwd, filename, NULL);
if (fullname == NULL)
  return(3);

if (strlen(data->moldy.lib_file) > 0)
  {
/* test for library existence */
  if (access(fullname, R_OK))
    gui_text_show(WARNING, "Potential library not found.\n");
  }

fp = fopen(fullname, "rt");
if (!fp)
  return(-1);

while (!fgetline(fp,line))
  {
  buff = tokenize(line, &num_tokens);
  if (*buff == NULL)
    return(-1);

  for (i = 0; potspec[i].name; i++)            /* Is 'name' a known type? */
    if (g_ascii_strcasecmp(*buff, potspec[i].name) == 0)
      break;

  if (potspec[i].name)                        /* Did the loop find 'name'? */
    {
    pflag++;
    break;
    }
  }
if (!pflag)
  {
  printf("Unrecognized potential function in %s \n", filename);
  return(-1);
  }

data->moldy.pot_type = i;                    /* yes                       */ 
n_items = potspec[i].npar;
p = g_malloc(n_items*sizeof(gdouble));

free_slist(data->moldy.pots);
data->moldy.pots = NULL;

/* Is a valid potential library specified? */
#if DEBUG_READ_POTS
  printf("Potential name %s (type %d)\n", potspec[i].name, data->moldy.pot_type);
  printf("No. of parms %d\n",n_items);
#endif

while (!fgetline(fp,line))
  {
  buff = tokenize(line, &num_tokens);
  if (buff)
    {
    if (!g_ascii_strcasecmp(*buff, "end"))
      break;

    line_no++;
    if (sscanf(line, "%4s %lf %4s %lf %[^#]", atom1, &charge1, atom2, &charge2, pline) > 4)
      {
      for(i = 0; i < n_items; i++)
        if (*(buff+4+i) == NULL)
          {
          printf("Insufficient parameters in parameter line %d\n", line_no);
          return(-1);
          }
        else
          p[i] = g_ascii_strtod(*(buff+4+i), NULL);

#if DEBUG_READ_POTS
  printf("%s %f ",atom1, charge1);
  printf("%s %f ",atom2, charge2);
  for(i=0; i<n_items; i++)
    printf("%f ",p[i]);
  printf("\n");
#endif
      /* Check validity of atoms */ 
      if (elem_symbol_test(g_strstrip(atom1)) && elem_symbol_test(g_strstrip(atom2)))
        {   
        pot = create_pot(atom1, charge1,
                         atom2, charge2);
        for(i = 0; i < n_items; i++)
          pot->parm[i] = (gdouble) *(p+i);

        data->moldy.pots = g_slist_append(data->moldy.pots, pot);
        }
      }
    }
  }
g_strfreev(buff);
g_free(p);
fclose(fp);

data->moldy.pots = g_slist_reverse(data->moldy.pots);

/* finished - only now we destroy the file dialog */
sysenv.file_type = DATA;
dialog_destroy_type(FILE_SELECT);

update_pot_dialog(data);
return(0);
}

#define DEBUG_WRITE_POTS 0
gint write_moldy_potentials(gchar *filename, struct model_pak *model)
{
gint i, ptype;
gchar *fullname;
FILE *fp;
struct pot_pak *pot;
GSList *plist, *slist1, *slist2;
struct atomtype_pak *site1=NULL, *site2=NULL;
struct elem_pak elem_data;

/* Atom identifiers */
gint code1, code2, code3, code4; /* Element codes */
gchar *atom1=NULL, *atom2=NULL;
gdouble charge1=0.0, charge2=0.0;
gboolean elem1_flag, elem2_flag; /* Flags to indicate label = element symbol */

/* checks */
if (!model)
  model = sysenv.active_model;

if (model == NULL)
  return(1);
if (filename == NULL)
  return(2);

if (!model->moldy.pots)
    return(-1);

model->atom_types = find_unique_sites(model);
ptype = model->moldy.pot_type;

/* Remember potential file directory */
g_free(model->moldy.lib_dir);
model->moldy.lib_dir = g_strdup(sysenv.cwd);

filename = parse_extension_set(filename, "pot");
fullname = g_build_filename(sysenv.cwd, filename, NULL);
if (fullname == NULL)
  return(3);

fp = fopen(fullname, "wt");
if (!fp)
  return(3);

/* Write file header */
  fprintf(fp, "# Potential parameters from structure \"%s\",", model->basename);
  fprintf(fp, " written by GDIS\n");
/* Write potential type */
  fprintf(fp, "%s\n", potspec[ptype].name);

/* Write each pair potential */
  for (plist = model->moldy.pots; plist; plist = g_slist_next(plist))
    {
    pot = plist->data;

    /* Is 1st symbol in pot file specific or general? */
    code1 = elem_symbol_test(g_strstrip(pot->site1_label));
    get_elem_data(code1, &elem_data, model);
    if (!g_ascii_strcasecmp(pot->site1_label, elem_data.symbol))
      elem1_flag = TRUE;
    else
      elem1_flag = FALSE;

    /* Is 2nd symbol in pot file specific or general? */
    code2 = elem_symbol_test(g_strstrip(pot->site2_label));
    get_elem_data(code2, &elem_data, model);
    if (!g_ascii_strcasecmp(pot->site2_label, elem_data.symbol))
      elem2_flag = TRUE;
    else
      elem2_flag = FALSE;

    for (slist1 = model->atom_types; slist1; slist1 = g_slist_next(slist1))
      {
      site1 = (struct atomtype_pak *) slist1->data;
      code3 = elem_symbol_test(g_strstrip(site1->label));

      for (slist2 = slist1; slist2; slist2 = g_slist_next(slist2))
        {
        site2 = (struct atomtype_pak *) slist2->data;
        code4 = elem_symbol_test(g_strstrip(site2->label));

        /* Reset identifiers */
        charge1 = 0.0;
        charge2 = 0.0;
        g_free(atom1);
        g_free(atom2);

        if (code1 == code3 && code2 == code4)
          {
          charge1 = site1->charge;
          charge2 = site2->charge;

          if (elem1_flag)
            {
            get_elem_data(code3, &elem_data, model);
            atom1 = g_strdup(elem_data.symbol);
            }
          else
            atom1 = g_strdup(site1->label);

          if (elem2_flag)
            {
            get_elem_data(code4, &elem_data, model);
            atom2 = g_strdup(elem_data.symbol);
            }
          else
            atom2 = g_strdup(site2->label);
          }
        else if (code1 == code4 && code2 == code3)
          {
          charge1 = site2->charge;
          charge2 = site1->charge;

          if (elem1_flag)
            {
            get_elem_data(code4, &elem_data, model);
            atom1 = g_strdup(elem_data.symbol);
            }
          else
            atom1 = g_strdup(site2->label);

          if (elem2_flag)
            {
            get_elem_data(code3, &elem_data, model);
            atom2 = g_strdup(elem_data.symbol);
            }
          else
            atom2 = g_strdup(site1->label);
          }
        else
          {
          atom1 = g_strdup("");
          atom2 = g_strdup("");
          }
        }
      }
  if (((!g_ascii_strcasecmp(pot->site1_label, atom1)
    && fabs(pot->site1_charge - charge1) < FRACTION_TOLERANCE)
    && (!g_ascii_strcasecmp(pot->site2_label, atom2)
    && fabs(pot->site2_charge - charge2) < FRACTION_TOLERANCE)))
    {
    fprintf(fp, "%s %f %s %f ", site1->label, site1->charge,
         site2->label, site2->charge);
    for(i=0; i<potspec[ptype].npar; i++)
      fprintf(fp, "%f ", pot->parm[i]);
    fprintf(fp, "\n");
    }

#if DEBUG_WRITE_POTS
  printf("%s %f %s %f", site1->label, site1->charge,
          site2->label, site2->charge);
  for(i=0; i<potspec[ptype].npar; i++)
    printf(" %f ", pot->parm[i]);
  printf("\n");
#endif
    }
  fprintf(fp, "end\n");
  fclose(fp);

/* finished - only now we destroy the file dialog */
gui_text_show(INFO, "Moldy potential file written successfully.\n");
dialog_destroy_type(FILE_SELECT);

g_free(atom1);
g_free(atom2);
g_free(model->moldy.lib_file);
model->moldy.lib_file = g_strdup(filename);

update_pot_dialog(model);
sysenv.file_type = DATA;

return(0);
}

/*******************************************/
/* write a MOLDY control file              */
/*******************************************/
#define DEBUG_WRITE_CONTROL 0
gint write_moldy_control(gchar *filename, struct model_pak *data)
{
struct moldy_pak moldy;
time_t tp;
FILE *fp;

if (data == NULL)
  return(1);
if (filename == NULL || strlen(filename) == 0)
  return(2);

moldy = data->moldy;

fp = fopen(filename, "wt");
if (!fp)
  {
  printf("Failed to open %s.\n", filename);
  return(2);
  }

#if DEBUG_WRITE_CONTROL
  printf("Writing control parameters to %s\n", filename);
#endif

tp = time(NULL);
fprintf(fp, "# Control file written by GDIS on %s", ctime((&tp)));
if (g_ascii_strcasecmp(moldy.title, ""))
  fprintf(fp, "title=%s\n", moldy.title);

if (strlen(moldy.restart_file) > 0)
    fprintf(fp, "restart-file=%s\n", moldy.restart_file);
else
  {
  if (data->periodic == 3)
    fprintf(fp, "lattice-start=%d\n", 1);
  else
    fprintf(fp, "density=%f\n", moldy.density);

  if (g_ascii_strcasecmp(moldy.control_file, moldy.sysspec_file))
    fprintf(fp, "sys-spec-file=%s\n", moldy.sysspec_file);
  }
/* Simulation parameters */
fprintf(fp, "nsteps=%d\n", moldy.nsteps);
fprintf(fp, "step=%g\n", moldy.deltat);
fprintf(fp, "cpu-limit=%g\n", moldy.cpu_limit);
fprintf(fp, "temperature=%g\n", moldy.temperature);

/* Save files */
if (moldy.text_save)
  fprintf(fp, "text-mode-save=%d\n", moldy.text_save);
if (strlen(moldy.save_file) > 0)
  fprintf(fp, "save-file=%s\n", moldy.save_file);
if (moldy.dump_level > 0 && strlen(moldy.dump_file) > 0)
  {
  fprintf(fp, "dump-file=%s\n", moldy.dump_file);
  fprintf(fp, "begin-dump=%d\n", moldy.dump_begin);
  fprintf(fp, "dump-interval=%d\n", moldy.dump_int);
  fprintf(fp, "dump-level=%d\n", moldy.dump_level);
  fprintf(fp, "ndumps=%d\n", moldy.maxdumps);
  }
if (g_ascii_strcasecmp(moldy.temp_file, "MDTEMPX"))
  fprintf(fp, "temp-file=%s\n", moldy.temp_file);
if (!moldy.xdr)
  fprintf(fp, "xdr=%d\n", moldy.xdr);

/* Temperature options */
fprintf(fp, "const-temp=%d\n",moldy.const_temp);
if (moldy.const_temp == NOSEHOOVER)
  {
  fprintf(fp, "ttmass=%g\n", moldy.ttmass);
  fprintf(fp, "rtmass=%g\n", moldy.rtmass);
  }
fprintf(fp, "scale-interval=%d\n", moldy.scale_int);
fprintf(fp, "scale-end=%d\n", moldy.scale_end);
fprintf(fp, "scale-options=%d\n", moldy.scale_options);

/* Constant pressure options */
fprintf(fp, "const-pressure=%d\n", moldy.const_press);
fprintf(fp, "pressure=%g\n", moldy.pressure);
fprintf(fp, "strain-mask=%d\n", moldy.strain_mask);
fprintf(fp, "w=%g\n", moldy.mass_parm);

/* Units */
fprintf(fp, "time-unit=%-.15g\n", moldy.time_unit);
fprintf(fp, "mass-unit=%-.15g\n", moldy.mass_unit);
fprintf(fp, "charge-unit=%-.15g\n", moldy.charge_unit);
fprintf(fp, "length-unit=%-.15g\n", moldy.length_unit);

/* Cutoff and Ewald sum */
if (moldy.auto_cutoff)
  {
  if (moldy.cutoff > 0)
     fprintf(fp, "cutoff=%g\n", moldy.cutoff);
  fputs("alpha=0\n", fp);
  if (moldy.kcutoff > 0)
     fprintf(fp, "k-cutoff=%g\n", moldy.kcutoff);
  }
else
  {
  if (moldy.real_only)
    {
    if (!moldy.ewald)
      fputs("alpha=-1\n", fp);
    else
      {
      fputs("alpha=1e-8\n", fp);
      }
    if (moldy.cutoff > 0)
      fprintf(fp, "cutoff=%g\n", moldy.cutoff);
    else
      fputs("cutoff=10.0\n", fp);
    }
  else
    {
    if (moldy.cutoff > 0)
      fprintf(fp, "cutoff=%g\n", moldy.cutoff);
    if (moldy.alpha > 0)
      fprintf(fp, "alpha=%g\n", moldy.alpha);
    if (moldy.kcutoff > 0)
      fprintf(fp, "k-cutoff=%g\n", moldy.kcutoff);
    }
  }
/* fprintf(fp, "ewald-accuracy=%g\n", moldy.ewald_accuracy); */

/* Averaging and backup options */
fprintf(fp, "roll-interval=%d\n", moldy.roll_int);
fprintf(fp, "print-interval=%d\n", moldy.print_int);
fprintf(fp, "begin-average=%d\n", moldy.av_start);
fprintf(fp, "average-interval=%d\n", moldy.av_int);
if (moldy.av_reset)
  fprintf(fp, "reset-averages=%d\n", moldy.av_reset);
fprintf(fp, "backup-interval=%d\n", moldy.backup_int);

/* Other options */
if (moldy.strict)
  fprintf(fp, "strict-cutoff=%d\n", 1);
if (moldy.surf)
  fprintf(fp, "surface-dipole=%d\n", 1);
if (moldy.subcell)
  fprintf(fp, "subcell=%g\n", moldy.subcell);

/* RDF settings */
fprintf(fp, "begin-rdf=%d\n", moldy.rdf_begin);
fprintf(fp, "rdf-interval=%d\n", moldy.rdf_int);
fprintf(fp, "rdf-out=%d\n", moldy.rdf_out);
fprintf(fp, "rdf-limit=%g\n", moldy.rdf_limit);
fprintf(fp, "nbins=%d\n", moldy.nbins);

/* Formatting, etc */
fprintf(fp, "seed=%f\n", moldy.seed);
fprintf(fp, "page-width=%d\n", moldy.page_width);
fprintf(fp, "page-length=%d\n", moldy.page_length);

fprintf(fp, "end\n");

/* done */
fclose(fp);

gui_text_show(INFO, "Moldy control file written successfully.\n");

#if DEBUG_WRITE_CONTROL
  printf("Moldy control file written successfully.\n");
#endif
return(0);
}
/*******************************************/
/* write a MOLDY system specification file */
/*******************************************/
#define DEBUG_WRITE_SYSSPEC 0
gint write_moldy_sys(gchar *filename, struct model_pak *data)
{
gint i, j, k;
gint ptype = data->moldy.pot_type;
gint num_species;
gint num_unique;
gint repeat, repeat_x, repeat_y, repeat_z;
gdouble x[3], quat[4];
time_t tp;
FILE *fp;
gint *pflag;
struct pot_pak *pdata;
struct moldyspec_pak *species = NULL;
struct mol_pak *mole;
struct site_pak *site;
struct atomtype_pak *type;
GSList *potlist = NULL, *mlist, *clist, *slist;
GSList *type_list;

if (data == NULL)
  return(1);

if (strlen(data->moldy.restart_file) > 0)
  return(1);

/* Break up molecules if bonds turned off */
if (!data->show_bonds)
  {
  data->build_molecules = FALSE;
  wipe_bonds(data);
  connect_refresh(data);
#if DEBUG_WRITE_SYSSPEC
puts("Bonds wiped.\n");
#endif
  }

repeat_x = data->image_limit[1]+data->image_limit[0];
repeat_y = data->image_limit[3]+data->image_limit[2];
repeat_z = data->image_limit[5]+data->image_limit[4];
repeat = repeat_x * repeat_y * repeat_z;

if (!g_ascii_strcasecmp(data->moldy.control_file, data->moldy.sysspec_file)) /* Append to control file */
  fp = fopen(filename, "a");
else
  fp = fopen(filename, "wt");

if (!fp)
  {
  printf("Failed to open %s.\n", filename);
  return(2);
  }

#if DEBUG_WRITE_SYSSPEC
if (!g_ascii_strcasecmp(data->moldy.control_file, data->moldy.sysspec_file)) /* Append to control file */
  printf("Appending sys-spec info to %s\n", filename);
else
  printf("Writing sys-spec info to %s\n", filename);
#endif

/* Determine no and type of different species */
/* data->atom_types = find_unique_sites(data); */
num_species = calc_species(data, POSITION_TOLERANCE);
num_unique = g_slist_length(data->atom_types);
pflag = g_malloc0(num_unique*num_unique*sizeof(gint));

if (data->periodic != 3)
  repeat *= data->moldy.num_mols;

tp = time(NULL);
fprintf(fp, "# System specification file written by GDIS on %s", ctime((&tp)));

fprintf(fp, "# %s\n", data->basename);
if (g_ascii_strcasecmp(data->moldy.title, ""))
  fprintf(fp, "# %s\n", data->moldy.title);

for (slist = data->moldy.species; slist; slist = g_slist_next(slist))
  {
  species = (struct moldyspec_pak *) slist->data;
  fprintf(fp, "%-14s %d", species->name, g_slist_length(species->mols)*repeat);
  if (species->framework || data->moldy.set_framework)
    fputs(" framework", fp);
  fputs("\n", fp);
  for (mlist = species->sites; mlist; mlist = g_slist_next(mlist))
    {
    site = (struct site_pak *) mlist->data;
    type = site->type;
    fprintf(fp, "%-5d %12.6f %12.6f %12.6f %12.4f %12.6f %s\n",
        type->id,
        site->x[0], site->x[1], site->x[2],
        type->mass*(AMU/data->moldy.mass_unit),
        type->charge*(ELCHARGE/data->moldy.charge_unit),
        type->label);
#if DEBUG_WRITE_SYSSPEC
  printf("Species %d: %s Site: id %d\n",
      g_slist_index(data->moldy.species, species), species->name, type->id);
#endif 
    }
  }
fprintf(fp, "end\n");

#if DEBUG_WRITE_SYSSPEC
  puts("Writing potential parameters.");
#endif 

/* If no potentials specified, assume same no of parameters as LJ */
if (ptype < 0 || ptype > 5)
  ptype = 0;

fprintf(fp, "%s\n", potspec[ptype].name);

/*  if (g_slist_length(data->moldy.pots) == 0)
    gui_text_show(WARNING, "No potentials specified for Moldy.\n");
*/
potlist = data->moldy.pots;
while(potlist != NULL)
  {
  i = j = -1;
  pdata = potlist->data;
  for (type_list = data->atom_types; type_list; type_list = g_slist_next(type_list))
    {
    type = type_list->data;
    if (!g_ascii_strcasecmp(pdata->site1_label, type->label) && pdata->site1_charge == type->charge)
       i = type->id;
    if (!g_ascii_strcasecmp(pdata->site2_label, type->label) && pdata->site2_charge == type->charge)
       j = type->id;
    }

  if (i > 0 && j > 0)  
    {
    fprintf(fp, "  %d   %d ", i, j);
    for(k=0; k < potspec[ptype].npar; k++)
      fprintf(fp, "%16.6f  ", pdata->parm[k]);
    pflag[(i-1)+(j-1)*num_unique]++;
    fputs("\n", fp);
    }
  potlist = g_slist_next(potlist);
  }
for(i=0; i < num_unique; i++)
  for(j=i; j < num_unique; j++)
    if (pflag[i+j*num_unique] == 0 && pflag[j+i*num_unique] == 0)
      {
      /* printf("Warning: No parameters found between sites %d and %d\n",
             i+1, j+1); */
      /* Set parameters to 0 if not found in parameter list */
      fprintf(fp, "  %d   %d ", i+1, j+1);
      for(k=0; k < potspec[ptype].npar; k++)
        fprintf(fp, "%16.6f  ", 0.0);
      fputs("\n", fp);
      }
fprintf(fp, "end\n");

if (data->periodic == 3)
  {
#if DEBUG_WRITE_SYSSPEC
  puts("Writing molecule coordinates.");
#endif 
  if (num_species > 0)
    {
    fprintf(fp, "%9.4f %9.4f %9.4f %9.4f  %9.4f  %9.4f  ",
       data->pbc[0], data->pbc[1], data->pbc[2],
          R2D*fabs(data->pbc[3]), R2D*fabs(data->pbc[4]), R2D*fabs(data->pbc[5]));

    fprintf(fp, "%d %d %d\n", repeat_x, repeat_y, repeat_z);

    for (mlist=data->moldy.species; mlist; mlist=g_slist_next(mlist))
      {
      species = (struct moldyspec_pak *) mlist->data;

      for (clist=species->mols; clist; clist=g_slist_next(clist))
        {
        mole = (struct mol_pak *) clist->data;

        if ( !data->moldy.set_framework)
          {
          calc_cofm(data, mole, x);
          vecmat(data->ilatmat, x);  /* Convert cofm to fractional coords */

        /* Ensure that coords are between [0,1) */
          for (j=3; j--;)
            {
            x[j] -= floor(x[j]);
            if (fabs(fmod(x[j],1)) < 5*FRACTION_TOLERANCE
                || fabs(x[j]-1.0)  < 5*FRACTION_TOLERANCE) 
              x[j] = 0.0;
            }
          }
        else
          VEC3SET(x, 0.0, 0.0, 0.0);

        fprintf(fp, "%-14s %10.5f %10.5f %10.5f ", species->name, x[0], x[1], x[2]);
        if (g_slist_length(species->sites) > 1 && !data->moldy.set_framework)
          {
          quat_fit(data, mole, species, quat, NJACOBI);
          fprintf(fp, "%10.5f %10.5f %10.5f %10.5f",
            quat[0], quat[1], quat[2], quat[3]);
          }
        fputs("\n", fp);
        }
      }
    }
  fprintf(fp, "end\n");
  }

/* done */
fclose(fp);

gui_text_show(INFO, "Moldy sys-spec file written successfully.\n");

if (!data->show_bonds)
  {
  data->build_molecules = TRUE;
  connect_refresh(data);
#if DEBUG_WRITE_SYSSPEC
puts("Bonds re-joined.\n");
#endif
  }

#if DEBUG_WRITE_SYSSPEC
  printf("Moldy sys-spec file written successfully.\n");
#endif
return(0);
}
