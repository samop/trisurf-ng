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
#include "stats.h"
#include "sh.h"
#include "vesicle.h"

ts_bool run_simulation(ts_vesicle *vesicle, ts_uint mcsweeps, ts_uint inititer, ts_uint iterations, ts_uint start_iteration){
	ts_uint i, j;
	ts_double r0;
	ts_double l1,l2,l3,volume=0.0,area=0.0,vmsr,bfsr, vmsrt, bfsrt;
	ts_ulong epochtime;
// 	char filename[255];
	FILE *fd=fopen("statistics.csv","w");
	if(fd==NULL){
		fatal("Cannot open statistics.csv file for writing",1);
	}
	fprintf(fd, "Epoch OuterLoop VertexMoveSucessRate BondFlipSuccessRate Volume Area lamdba1 lambda2 lmabda3\n");
	centermass(vesicle);
	cell_occupation(vesicle);
	if(start_iteration<inititer) ts_fprintf(stdout, "Starting simulation (first %d x %d MC sweeps will not be recorded on disk)\n", inititer, mcsweeps);
	for(i=start_iteration;i<inititer+iterations;i++){
		vmsr=0.0;
		bfsr=0.0;
		for(j=0;j<mcsweeps;j++){
			single_timestep(vesicle, &vmsrt, &bfsrt);
			vmsr+=vmsrt;
			bfsr+=bfsrt;
		}
		vmsr/=(ts_double)mcsweeps;
		bfsr/=(ts_double)mcsweeps;
		centermass(vesicle);
		cell_occupation(vesicle);
		ts_fprintf(stdout,"Done %d out of %d iterations (x %d MC sweeps).\n",i+1,inititer+iterations,mcsweeps);
            dump_state(vesicle,i);
		if(i>=inititer){
			write_vertex_xml_file(vesicle,i-inititer);
			write_master_xml_file("test.pvd");
			epochtime=get_epoch();			
			gyration_eigen(vesicle, &l1, &l2, &l3);
			get_area_volume(vesicle, &area,&volume);
			vesicle_volume(vesicle);
			r0=getR0(vesicle);
			preparationSh(vesicle,r0);
			calculateYlmi(vesicle);
			calculateUlm(vesicle);
			storeUlm2(vesicle);
			saveAvgUlm2(vesicle);

			fprintf(fd, "%lu %u %e %e %e %e %e %e %e\n",epochtime,i,vmsr,bfsr,volume, area,l1,l2,l3);
			
		//	sprintf(filename,"timestep-%05d.pov",i-inititer);
		//	write_pov_file(vesicle,filename);
		}
	}
	fclose(fd);
	return TS_SUCCESS;
}

ts_bool single_timestep(ts_vesicle *vesicle,ts_double *vmsr, ts_double *bfsr){
    ts_bool retval;
    ts_double rnvec[3];
    ts_uint i,j,b;
    ts_uint vmsrcnt=0;
    for(i=0;i<vesicle->vlist->n;i++){
        rnvec[0]=drand48();
        rnvec[1]=drand48();
        rnvec[2]=drand48();
        retval=single_verticle_timestep(vesicle,vesicle->vlist->vtx[i],rnvec);
	if(retval==TS_SUCCESS) vmsrcnt++;        
    }

	ts_int bfsrcnt=0;
    for(i=0;i<3*vesicle->vlist->n;i++){
	b=rand() % vesicle->blist->n;
        //find a bond and return a pointer to a bond...
        //call single_bondflip_timestep...
        retval=single_bondflip_timestep(vesicle,vesicle->blist->bond[b],rnvec);
	if(retval==TS_SUCCESS) bfsrcnt++;        
    }

	for(i=0;i<vesicle->poly_list->n;i++){
		for(j=0;j<vesicle->poly_list->poly[i]->vlist->n;j++){
			rnvec[0]=drand48();
			rnvec[1]=drand48();
			rnvec[2]=drand48();
			retval=single_poly_vertex_move(vesicle,vesicle->poly_list->poly[i],vesicle->poly_list->poly[i]->vlist->vtx[j],rnvec);	
		}
	}


	for(i=0;i<vesicle->filament_list->n;i++){
		for(j=0;j<vesicle->filament_list->poly[i]->vlist->n;j++){
			rnvec[0]=drand48();
			rnvec[1]=drand48();
			rnvec[2]=drand48();
			retval=single_filament_vertex_move(vesicle,vesicle->filament_list->poly[i],vesicle->filament_list->poly[i]->vlist->vtx[j],rnvec);	
		}
	}
 

//	printf("Bondflip success rate in one sweep: %d/%d=%e\n", cnt,3*vesicle->blist->n,(double)cnt/(double)vesicle->blist->n/3.0);
	*vmsr=(ts_double)vmsrcnt/(ts_double)vesicle->vlist->n;
	*bfsr=(ts_double)bfsrcnt/(ts_double)vesicle->vlist->n/3.0;
    return TS_SUCCESS;
}



