#include<stdio.h>
#include<math.h>
#include "general.h"
#include "vertex.h"
#include "bond.h"
#include "triangle.h"
#include "cell.h"
#include "vesicle.h"
#include "io.h"
#include "initial_distribution.h"
#include "frame.h"
#include "timestep.h"

/** Entrance function to the program
  * @param argv is a number of parameters used in program call (including the program name
  * @param argc is a pointer to strings (character arrays) which holds the arguments
  * @returns returns 0 on success, any other number on fail.
*/

int main(int argv, char *argc[]){
ts_bool retval;
    ts_vertex_list *vlist=init_vertex_list(5);
ts_vertex_list *vlist1;
ts_bond_list *blist=init_bond_list();
ts_triangle_list *tlist=init_triangle_list();
ts_cell_list *clist=init_cell_list(3,3,3,0.3);

retval=vtx_add_cneighbour(blist,vlist->vtx[1],vlist->vtx[0]);
if(retval==TS_FAIL) printf("1. already a member or vertex is null!\n");

retval=vtx_add_neighbour(vlist->vtx[0],vlist->vtx[1]);
if(retval==TS_FAIL) printf("2. already a member or vertex is null!\n");
fprintf(stderr,"Was here");
retval=vertex_list_remove_vtx(vlist->vtx[1]->neigh,vlist->vtx[0]);
retval=vertex_list_remove_vtx(vlist->vtx[0]->neigh,vlist->vtx[1]);
vtx_add_neighbour(vlist->vtx[0],vlist->vtx[1]);
fprintf(stderr,"Was here too!\n");

vlist->vtx[0]->x=1.0;
vlist->vtx[0]->x=1.1;
vlist1=vertex_list_copy(vlist);
bond_add(blist, vlist->vtx[1],vlist->vtx[0]);
triangle_add(tlist,vlist->vtx[1],vlist->vtx[2],vlist->vtx[3]);

triangle_add(tlist,vlist->vtx[1],vlist->vtx[2],vlist->vtx[3]);

printf("Cell idx=1 has vertices=%u\n",clist->cell[0]->nvertex);
cell_add_vertex(clist->cell[0], vlist->vtx[0]);
printf("Cell idx=1 has vertices=%u\n",clist->cell[0]->nvertex);
printf("Cell idx=1 has vertex[0] has x coordinate=%e\n",clist->cell[0]->vertex[0]->x);
cell_list_cell_occupation_clear(clist);
printf("Cell idx=1 has vertices=%u\n",clist->cell[0]->nvertex);
cell_add_vertex(clist->cell[0], vlist->vtx[0]);


triangle_list_free(tlist);
bond_list_free(blist);
vtx_list_free(vlist);
cell_list_free(clist);

vtx_list_free(vlist1);
printf("Tests complete.\n");

return 0; //program finished perfectly ok. We return 0.
} 
