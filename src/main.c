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
#include "poly.h"

/** Entrance function to the program
  * @param argv is a number of parameters used in program call (including the program name
  * @param argc is a pointer to strings (character arrays) which holds the arguments
  * @returns returns 0 on success, any other number on fail.
*/

int main(int argv, char *argc[]){
ts_uint inititer,mcsweeps, iterations;
ts_vesicle *vesicle;
/* THIS SHOULD GO INTO UNIT TEST
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
retval=vtx_remove_neighbour(vlist->vtx[1],vlist->vtx[0]);
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
*/
vesicle=parsetape(&mcsweeps, &inititer, &iterations);

/*Testing */
//vesicle->poly_list=init_poly_list(1400,20,vesicle->vlist);

//poly_list_free(vesicle->poly_list);
/*End testing*/
printf("Vertex %d ima x komponento %e\n",123,vesicle->vlist->vtx[123]->x);
dump_state(vesicle);
vesicle->vlist->vtx[123]->x=123.0;
printf("Vertex %d ima x komponento %e\n",123,vesicle->vlist->vtx[123]->x);
write_vertex_xml_file(vesicle,0);
vesicle=restore_state();
printf("Stevilo vertexov je %d\n",vesicle->vlist->n);
printf("Vertex %d ima x komponento %e\n",123,vesicle->vlist->vtx[123]->x);
write_vertex_xml_file(vesicle,1);
write_master_xml_file("test.pvd");
return 0;

run_simulation(vesicle, mcsweeps, inititer, iterations);

write_master_xml_file("test.pvd");
write_dout_fcompat_file(vesicle,"dout");
vesicle_free(vesicle);

return 0; //program finished perfectly ok. We return 0.
} 
