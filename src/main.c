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
ts_int i;
ts_bool retval;
ts_vertex **vlist=init_vertex_list(5);


retval=vtx_add_neighbour(&vlist[0],&vlist[1]);
if(retval==TS_FAIL) printf("1. already a member or vertex is null!\n");
retval=vtx_add_neighbour(&vlist[1],&vlist[0]);
if(retval==TS_FAIL) printf("2. already a member or vertex is null!\n");

for(i=0;i<5;i++){
vtx_free(&vlist[i]);
}
free(vlist);
printf("Done.\n");
return 0; //program finished perfectly ok. We return 0.
} 
