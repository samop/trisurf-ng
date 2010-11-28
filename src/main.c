#include<stdio.h>
#include<math.h>
#include "general.h"
#include "vertex.h"
#include "bond.h"
#include "triangle.h"
//#include "io.h"
//#include "initial_timestep.h"

/** Entrance function to the program
  * @param argv is a number of parameters used in program call (including the program name
  * @param argc is a pointer to strings (character arrays) which holds the arguments
  * @returns returns 0 on success, any other number on fail.
*/

int main(int argv, char *argc[]){
ts_bool retval;
ts_vertex_list *vlist=init_vertex_list(5);
ts_bond_list *blist=init_bond_list();
ts_triangle_list *tlist=init_triangle_list();



retval=vtx_add_cneighbour(blist,vlist->vtx[1],vlist->vtx[0]);
if(retval==TS_FAIL) printf("1. already a member or vertex is null!\n");

retval=vtx_add_cneighbour(blist,vlist->vtx[0],vlist->vtx[1]);
if(retval==TS_FAIL) printf("2. already a member or vertex is null!\n");

retval=vtx_remove_neighbour(vlist->vtx[0],vlist->vtx[1]);
vtx_add_neighbour(vlist->vtx[0],vlist->vtx[1]);

vlist->vtx[0]->data->x=1.0;
vlist->vtx[0]->data->x=1.1;

bond_add(blist, vlist->vtx[1],vlist->vtx[0]);
triangle_add(tlist,vlist->vtx[1],vlist->vtx[2],vlist->vtx[3]);

triangle_add(tlist,vlist->vtx[1],vlist->vtx[2],vlist->vtx[3]);

triangle_list_free(tlist);
bond_list_free(blist);
vtx_list_free(vlist);
printf("Done.\n");
return 0; //program finished perfectly ok. We return 0.
} 
