#include<stdlib.h>
#include<stdio.h>
#include<math.h>
//#include "io.h"
#include "general.h"
#include "timestep.h"
#include "vertexmove.h"
#include "bondflip.h"
#include "frame.h"
#include "io.h"
ts_bool run_simulation(ts_vesicle *vesicle, ts_uint mcsweeps, ts_uint inititer, ts_uint iterations){
	ts_uint i, j;

	centermass(vesicle);
	cell_occupation(vesicle);
	ts_fprintf(stdout, "Starting simulation (first %d x %d MC sweeps will not be recorded on disk)\n", inititer, mcsweeps);
	for(i=0;i<inititer+iterations;i++){
		for(j=0;j<mcsweeps;j++){
			single_timestep(vesicle);
		}
		centermass(vesicle);
		cell_occupation(vesicle);
		if(i>inititer){
			write_vertex_xml_file(vesicle,i-inititer);
		}
	}
	return TS_SUCCESS;
}

ts_bool single_timestep(ts_vesicle *vesicle){
    ts_bool retval;
    ts_double rnvec[3];
    ts_uint i, b;
    for(i=0;i<vesicle->vlist->n;i++){
        rnvec[0]=drand48();
        rnvec[1]=drand48();
        rnvec[2]=drand48();
        retval=single_verticle_timestep(vesicle,vesicle->vlist->vtx[i],rnvec);
    }

//	ts_int cnt=0;
    for(i=0;i<vesicle->vlist->n;i++){
//why is rnvec needed in bondflip?
/*        rnvec[0]=drand48();
        rnvec[1]=drand48();
        rnvec[2]=drand48();
*/ 
	b=rand() % vesicle->blist->n;
        //find a bond and return a pointer to a bond...
        //call single_bondflip_timestep...
        retval=single_bondflip_timestep(vesicle,vesicle->blist->bond[b],rnvec);
//	if(retval==TS_SUCCESS) cnt++;        
    } 
//	printf("Bondflip success rate in one sweep: %d/%d=%e\n", cnt,vesicle->blist->n,(double)cnt/(double)vesicle->blist->n);
	if(retval);
    return TS_SUCCESS;
}



