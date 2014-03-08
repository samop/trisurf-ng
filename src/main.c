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
	ts_vesicle *vesicle;
	ts_tape *tape;
	ts_uint start_iteration=0;
	parse_args(argv,argc);
	ts_fprintf(stdout,"\nStarting program...\n\n");
if(force_from_tape){
ts_fprintf(stdout,"****************************************************\n");
ts_fprintf(stdout,"**** Reinitializing initial geometry from tape *****\n");
ts_fprintf(stdout,"****************************************************\n\n");
tape=parsetape("tape");
vesicle=create_vesicle_from_tape(tape);
} else {

ts_fprintf(stdout,"**********************************************************************\n");
ts_fprintf(stdout,"**** Recreating vesicle from dump file and continuing simulation *****\n");
ts_fprintf(stdout,"**********************************************************************\n\n");
tape=parsetape("tape");
vesicle=restore_state(&start_iteration);
}

run_simulation(vesicle, tape->mcsweeps, tape->inititer, tape->iterations);
write_master_xml_file("test.pvd");
write_dout_fcompat_file(vesicle,"dout");
vesicle_free(vesicle);
tape_free(tape);
return 0; //program finished perfectly ok. We return 0.
} 
