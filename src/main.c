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
	force_from_tape=0;
	parse_args(argv,argc); // sets global variable command_line_args (defined in io.h)
	ts_fprintf(stdout,"Starting program...\n\n");
	if(command_line_args.force_from_tape){
		ts_fprintf(stdout,"************************************************\n");
		ts_fprintf(stdout,"**** Generating initial geometry from tape *****\n");
		ts_fprintf(stdout,"************************************************\n\n");
		tape=parsetape("tape");
		vesicle=create_vesicle_from_tape(tape);
	} else {

		ts_fprintf(stdout,"**********************************************************************\n");
		ts_fprintf(stdout,"**** Recreating vesicle from dump file and continuing simulation *****\n");
		ts_fprintf(stdout,"**********************************************************************\n\n");
		tape=parsetape("tape");
		vesicle=restore_state(&start_iteration);
		// nove vrednosti iz tapea...
/*	vesicle->bending_stiffness=tape->xk0;
		set_global_values(tape->xk0);
		vesicle->pressure=tape->pressure;
*/
		if(command_line_args.reset_iteration_count) start_iteration=tape->inititer+1;
		else start_iteration++;

		if(start_iteration>=tape->iterations){
			ts_fprintf(stdout, "Simulation already completed. if you want to rerun it try with --force-from-tape or --reset-iteration-count\n\n");
			return 0;
		}
	}

	run_simulation(vesicle, tape->mcsweeps, tape->inititer, tape->iterations, start_iteration);
	write_master_xml_file("test.pvd");
	write_dout_fcompat_file(vesicle,"dout");
	vesicle_free(vesicle);
	tape_free(tape);
	return 0; //program finished perfectly ok. We return 0.
} 
