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
#include "shcomplex.h"
#include "vesicle.h"
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>


ts_bool run_simulation(ts_vesicle *vesicle, ts_uint mcsweeps, ts_uint inititer, ts_uint iterations, ts_uint start_iteration){
	ts_uint i, j,k,l,m;
	ts_double r0,kc1,kc2,kc3,kc4;
	ts_double l1,l2,l3,volume=0.0,area=0.0,vmsr,bfsr, vmsrt, bfsrt;
	ts_ulong epochtime;
	FILE *fd1,*fd2=NULL;
// 	char filename[255];
	FILE *fd=fopen("statistics.csv","w");
	if(fd==NULL){
		fatal("Cannot open statistics.csv file for writing",1);
	}
	fprintf(fd, "Epoch OuterLoop VertexMoveSucessRate BondFlipSuccessRate Volume Area lamdba1 lambda2 lambda3 Kc(2-9) Kc(6-9) Kc(2-end) Kc(3-6)\n");

	 if(vesicle->sphHarmonics!=NULL){
		fd2=fopen("ulm2.csv","w");
		if(fd2==NULL){
			fatal("Cannot open ulm2.csv file for writing",1);
		}
		fprintf(fd2, "Timestep u_00^2 u_10^2 u_11^2 u_20^2 ...\n");	

	}



	centermass(vesicle);
	cell_occupation(vesicle);
	vesicle_volume(vesicle); //needed for constant volume at this moment
	if(start_iteration<inititer) ts_fprintf(stdout, "Starting simulation (first %d x %d MC sweeps will not be recorded on disk)\n", inititer, mcsweeps);
	for(i=start_iteration;i<inititer+iterations;i++){
		vmsr=0.0;
		bfsr=0.0;
/*    vesicle_volume(vesicle);
    fprintf(stderr,"Volume before TS=%1.16e\n", vesicle->volume); */
		for(j=0;j<mcsweeps;j++){
			single_timestep(vesicle, &vmsrt, &bfsrt);
			vmsr+=vmsrt;
			bfsr+=bfsrt;
		}
/*
    vesicle_volume(vesicle);
    fprintf(stderr,"Volume after TS=%1.16e\n", vesicle->volume); */
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
			vesicle_volume(vesicle); //calculates just volume. Area is not added to ts_vesicle yet!
			get_area_volume(vesicle, &area,&volume); //that's why I must recalculate area (and volume for no particular reason).
			r0=getR0(vesicle);
            if(vesicle->sphHarmonics!=NULL){
			    preparationSh(vesicle,r0);
			    //calculateYlmi(vesicle);
			    calculateUlmComplex(vesicle);
			    storeUlmComplex2(vesicle);
			    saveAvgUlm2(vesicle);
                kc1=calculateKc(vesicle, 2,9);
                kc2=calculateKc(vesicle, 6,9);
                kc3=calculateKc(vesicle, 2,vesicle->sphHarmonics->l);
                kc4=calculateKc(vesicle, 3,6);
            
				fd1=fopen("state.dat","w");
				fprintf(fd1,"%e %e\n",vesicle->volume, getR0(vesicle));
				for(k=0;k<vesicle->vlist->n;k++){
					fprintf(fd1,"%e %e %e %e %e\n",
						vesicle->vlist->vtx[k]->x,
						vesicle->vlist->vtx[k]->y,
						vesicle->vlist->vtx[k]->z,
						vesicle->vlist->vtx[k]->solAngle,
						vesicle->vlist->vtx[k]->relR
					);
				}
				fclose(fd1);
		
			fprintf(fd2,"%u ", i);
			for(l=0;l<vesicle->sphHarmonics->l;l++){
				for(m=l;m<2*l+1;m++){
					fprintf(fd2,"%e ", gsl_complex_abs2(vesicle->sphHarmonics->ulmComplex[l][m]) );
				}
			}
				fprintf(fd2,"\n");
	
		    	fflush(fd2);	

            }

			fprintf(fd, "%lu %u %e %e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e\n",epochtime,i,vmsr,bfsr,volume, area,l1,l2,l3,kc1, kc2, kc3,kc4);

		    fflush(fd);	
		//	sprintf(filename,"timestep-%05d.pov",i-inititer);
		//	write_pov_file(vesicle,filename);
		}
	}
	fclose(fd);
	if(fd2!=NULL) fclose(fd2);
	return TS_SUCCESS;
}

ts_bool single_timestep(ts_vesicle *vesicle,ts_double *vmsr, ts_double *bfsr){
//    vesicle_volume(vesicle);
//    fprintf(stderr,"Volume before TS=%1.16e\n", vesicle->volume);
    ts_bool retval;
    ts_double rnvec[3];
    ts_uint i,j, b;
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
       //     b++; retval=TS_FAIL;
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
//    vesicle_volume(vesicle);
//    fprintf(stderr,"Volume after TS=%1.16e\n", vesicle->volume);
    return TS_SUCCESS;
}



