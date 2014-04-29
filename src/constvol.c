#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include "general.h"
#include "constvol.h"
#include "triangle.h"
#include "energy.h"
#include "vertex.h"
#include "cell.h"

ts_bool constvolume(ts_vesicle *vesicle, ts_vertex *vtx_avoid, ts_double Vol, ts_double *retEnergy, ts_vertex **vtx_moved_retval, ts_vertex **vtx_backup){
    ts_vertex *vtx_moved;
    ts_uint vtxind,i,j;
    ts_uint Ntries=20;
	ts_vertex *backupvtx;
    ts_double Rv, dh, dvol, voldiff, oenergy,delta_energy;
    backupvtx=(ts_vertex *)calloc(sizeof(ts_vertex),10);
    ts_double l0 = (1.0 + sqrt(vesicle->dmax))/2.0; //make this a global constant if necessary
    for(i=0;i<Ntries;i++){
        vtxind=rand() % vesicle->vlist->n;
        vtx_moved=vesicle->vlist->vtx[vtxind];
        if(vtx_moved==vtx_avoid) continue;

        for(j=0;j<vtx_moved->neigh_no;j++){
            if(vtx_moved->neigh[j]==vtx_avoid) continue;
/*            for(k=0;k<vtx_moved->neigh[j]->neigh_no;k++){
                if(vtx_moved->neigh[j]->neigh[k]==vtx_avoid) continue;
            }   
*/

        }
         
	    memcpy((void *)&backupvtx[0],(void *)vtx_moved,sizeof(ts_vertex));
        //move vertex in specified direction. first try, test move!

        Rv=sqrt(pow(vtx_moved->x,2)+pow(vtx_moved->y,2)+pow(vtx_moved->z,2));
        dh=2.0*Vol/(sqrt(3.0)*l0*l0);
//        fprintf(stderr,"Prej (x,y,z)=(%e,%e,%e).\n",vtx_moved->x,vtx_moved->y,vtx_moved->z);
	    vtx_moved->x=vtx_moved->x*(1.0-dh/Rv);
	    vtx_moved->y=vtx_moved->y*(1.0-dh/Rv);
	    vtx_moved->z=vtx_moved->z*(1.0-dh/Rv);
//        fprintf(stderr,"Potem (x,y,z)=(%e,%e,%e). Vol=%e\n",vtx_moved->x,vtx_moved->y,vtx_moved->z,Vol);

        //check for constraints
          if(constvolConstraintCheck(vesicle, vtx_moved)==TS_FAIL){
		    vtx_moved=memcpy((void *)vtx_moved,(void *)&backupvtx[0],sizeof(ts_vertex));
            continue;
        }
//        fprintf(stderr,"Sprejet.\n");

        // All checks OK!
            fprintf(stderr, "Step 1 success\n");

        for(j=0;j<vtx_moved->neigh_no;j++){
        	memcpy((void *)&backupvtx[j+1],(void *)vtx_moved->neigh[j],sizeof(ts_vertex));
	    }
        dvol=0.0;
        for(j=0;j<vtx_moved->tristar_no;j++){
            dvol-=vtx_moved->tristar[j]->volume;
            triangle_normal_vector(vtx_moved->tristar[j]);
            dvol+=vtx_moved->tristar[j]->volume;
        }

        voldiff=dvol-Vol;

        if(fabs(voldiff)/vesicle->volume < vesicle->tape->constvolprecision){
            //calculate energy, return change in energy...
             oenergy=vtx_moved->energy;
            energy_vertex(vtx_moved);
            delta_energy=vtx_moved->xk*(vtx_moved->energy - oenergy);
            //the same is done for neighbouring vertices
            for(i=0;i<vtx_moved->neigh_no;i++){
                oenergy=vtx_moved->neigh[i]->energy;
                energy_vertex(vtx_moved->neigh[i]);
                delta_energy+=vtx_moved->neigh[i]->xk*(vtx_moved->neigh[i]->energy-oenergy);
            }
            *retEnergy=delta_energy;
            *vtx_backup=backupvtx;
            *vtx_moved_retval=vtx_moved;
            fprintf(stderr, "Preliminary success\n");
            return TS_SUCCESS;
        }        
            fprintf(stderr, "Step 2 success\n");
        //do it again ;)
        dh=Vol*dh/dvol;
		vtx_moved=memcpy((void *)vtx_moved,(void *)&backupvtx[0],sizeof(ts_vertex));
        vtx_moved->x=vtx_moved->x*(1-dh/Rv);
	    vtx_moved->y=vtx_moved->y*(1-dh/Rv);
	    vtx_moved->z=vtx_moved->z*(1-dh/Rv);
        //check for constraints
        if(constvolConstraintCheck(vesicle, vtx_moved)==TS_FAIL){
	        for(j=0;j<vtx_moved->neigh_no;j++){
        	    memcpy((void *)vtx_moved->neigh[j],(void *)&backupvtx[j+1],sizeof(ts_vertex));
	        }
	        vtx_moved=memcpy((void *)vtx_moved,(void *)&backupvtx[0],sizeof(ts_vertex));
            continue;
        }

        dvol=0.0;
        for(j=0;j<vtx_moved->tristar_no;j++){
            dvol-=vtx_moved->tristar[j]->volume;
            triangle_normal_vector(vtx_moved->tristar[j]);
            dvol+=vtx_moved->tristar[j]->volume;
        }

            fprintf(stderr, "Step 3a success voldiff=%e\n",voldiff);
        voldiff=dvol-Vol;
            fprintf(stderr, "Step 3b success voldiff=%e\n",voldiff);
        if(fabs(voldiff)/vesicle->volume < vesicle->tape->constvolprecision){
            //calculate energy, return change in energy...
            oenergy=vtx_moved->energy;
            energy_vertex(vtx_moved);
            delta_energy=vtx_moved->xk*(vtx_moved->energy - oenergy);
            //the same is done for neighbouring vertices
            for(i=0;i<vtx_moved->neigh_no;i++){
                oenergy=vtx_moved->neigh[i]->energy;
                energy_vertex(vtx_moved->neigh[i]);
                delta_energy+=vtx_moved->neigh[i]->xk*(vtx_moved->neigh[i]->energy-oenergy);
            }
            *retEnergy=delta_energy;
            *vtx_backup=backupvtx;
            *vtx_moved_retval=vtx_moved;
            fprintf(stderr, "DVOL=%e\n",voldiff);
            return TS_SUCCESS;
        }        


    }
    free(backupvtx);
            fprintf(stderr, "fail\n");
    return TS_FAIL;
}


ts_bool constvolConstraintCheck(ts_vesicle *vesicle, ts_vertex *vtx){ 
        ts_uint i;
        ts_double dist;
        ts_uint cellidx;
    	//distance with neighbours check
        for(i=0;i<vtx->neigh_no;i++){
            dist=vtx_distance_sq(vtx,vtx->neigh[i]);
            if(dist<1.0 || dist>vesicle->dmax) {
		    return TS_FAIL;
		    }
        }
        // Distance with grafted poly-vertex check:	
	    if(vtx->grafted_poly!=NULL){
		    dist=vtx_distance_sq(vtx,vtx->grafted_poly->vlist->vtx[0]);
            if(dist<1.0 || dist>vesicle->dmax) {
		    return TS_FAIL;
		    }
	    }

        // Nucleus penetration check:
	    if (vtx->x*vtx->x + vtx->y*vtx->y + vtx->z*vtx->z < vesicle->R_nucleus){
		    return TS_FAIL;
	    }

        //self avoidance check with distant vertices
	    cellidx=vertex_self_avoidance(vesicle, vtx);
    	//check occupation number
	    return cell_occupation_number_and_internal_proximity(vesicle->clist,cellidx,vtx);
}



ts_bool constvolumerestore(ts_vertex *vtx_moved,ts_vertex *vtx_backup){
    ts_uint j;
	 memcpy((void *)vtx_moved,(void *)&vtx_backup[0],sizeof(ts_vertex));
     for(j=0;j<vtx_moved->neigh_no;j++){
        	    memcpy((void *)vtx_moved->neigh[j],(void *)&vtx_backup[j+1],sizeof(ts_vertex));
	}
    free(vtx_backup);
    return TS_SUCCESS;
}

ts_bool constvolumeaccept(ts_vesicle *vesicle,ts_vertex *vtx_moved, ts_vertex *vtx_backup){
    ts_bool retval;
	ts_uint cellidx=vertex_self_avoidance(vesicle, vtx_moved);
    if(vtx_moved->cell!=vesicle->clist->cell[cellidx]){
		retval=cell_add_vertex(vesicle->clist->cell[cellidx],vtx_moved);
		if(retval==TS_SUCCESS) cell_remove_vertex(vtx_backup[0].cell,vtx_moved);
		
	}
    free(vtx_backup);

    return TS_SUCCESS;
}
