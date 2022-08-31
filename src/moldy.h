/*
Copyright (C) 2003 by Sean David Fleming
Copyright (C) 2015 by Craig Andrew James Fisher

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

/* Site specification */
enum {BY_ATOM_LABEL, BY_ELEM_SYMBOL};

/* Pressure control */
enum {NOPRESS, PRSTRESS, PRPRESS, WCSTRESS, ANDERSEN};

/* Temperature control */
enum {SCALING, NOSEHOOVER, GAUSSIAN};

/* Ewald summation */
enum {EWALD_INC, REAL_INC};

#define MAX_POT_PARMS 8 /* Maximum no of potential parms */
#define ENERGY_TOLERANCE 1e-3

#define NJACOBI 30  /* Max number of Jacobi sweeps */

typedef struct      /* Holds moldy potential function information */
{
   gchar *name;
   gint  npar;
} pots_pak;

GSList *find_unique_sites (struct model_pak *);
gint calc_species (struct model_pak *, gdouble);

void moldy_load_dialog(gchar *, gint, gpointer, struct model_pak *, gint);
gint delete_moldy_potentials(struct model_pak *);
void moldy_files_init(struct model_pak *);
void connect_file_bonds(GSList *, struct model_pak *);
void update_pot_dialog(struct model_pak *data);
