#include<stdio.h>
#include<math.h>
#include "general.h"
#include "vertex.h"
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


retval=vtx_add_neighbour(VTX(1),VTX(0));
if(retval==TS_FAIL) printf("1. already a member or vertex is null!\n");
retval=vtx_add_neighbour(VTX(0),VTX(1));
if(retval==TS_FAIL) printf("2. already a member or vertex is null!\n");
VTX_DATA(1)->x=1.0;
vtx_list_free(vlist);
printf("Done.\n");
return 0; //program finished perfectly ok. We return 0.
} 
