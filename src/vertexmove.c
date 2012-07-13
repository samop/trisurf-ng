#include<stdlib.h>
#include<math.h>
#include "general.h"
#include "vertex.h"
#include "bond.h"
#include "triangle.h"
#include "vesicle.h"
#include "energy.h"
#include "timestep.h"
#include "cell.h"
//#include "io.h"
#include<stdio.h>
#include "vertexmove.h"
#include <string.h>

ts_bool single_verticle_timestep(ts_vesicle *vesicle,ts_vertex *vtx,ts_double
*rn){
    ts_uint i;
    ts_double dist;
    ts_bool retval; 
    ts_uint cellidx; 
    ts_double delta_energy,oenergy;
	//This will hold all the information of vtx and its neighbours
	ts_vertex backupvtx[20];
	memcpy((void *)&backupvtx[0],(void *)vtx,sizeof(ts_vertex));

	//Some stupid tests for debugging cell occupation!
/*     	cellidx=vertex_self_avoidance(vesicle, vtx);
	if(vesicle->clist->cell[cellidx]==vtx->cell){
		fprintf(stderr,"Idx match!\n");
	} else {
		fprintf(stderr,"***** Idx don't match!\n");
		fatal("ENding.",1);
	}
*/

    	//temporarly moving the vertex
	vtx->x=vtx->x+vesicle->stepsize*(2.0*rn[0]-1.0);
    	vtx->y=vtx->y+vesicle->stepsize*(2.0*rn[1]-1.0);
    	vtx->z=vtx->z+vesicle->stepsize*(2.0*rn[2]-1.0);
    	//distance with neighbours check
    for(i=0;i<vtx->neigh_no;i++){
        dist=vtx_distance_sq(vtx,vtx->neigh[i]);
        if(dist<1.0 || dist>vesicle->dmax) {
		vtx=memcpy((void *)vtx,(void *)&backupvtx[0],sizeof(ts_vertex));
		return TS_FAIL;
		}
    }
    //self avoidance check with distant vertices
     cellidx=vertex_self_avoidance(vesicle, vtx);
    //check occupation number
     retval=cell_occupation_number_and_internal_proximity(vesicle->clist,cellidx,vtx);
	
    if(retval==TS_FAIL){
		vtx=memcpy((void *)vtx,(void *)&backupvtx[0],sizeof(ts_vertex));
        return TS_FAIL;
    } 
   
 
    //if all the tests are successful, then energy for vtx and neighbours is calculated
	for(i=0;i<vtx->neigh_no;i++){
	memcpy((void *)&backupvtx[i+1],(void *)vtx->neigh[i],sizeof(ts_vertex));
	}



    delta_energy=0;
    //update the normals of triangles that share bead i.
    for(i=0;i<vtx->tristar_no;i++) triangle_normal_vector(vtx->tristar[i]);
	oenergy=vtx->energy;
    energy_vertex(vtx);
    delta_energy=vtx->xk*(vtx->energy - oenergy);
    //the same is done for neighbouring vertices
    for(i=0;i<vtx->neigh_no;i++){
        oenergy=vtx->neigh[i]->energy;
        energy_vertex(vtx->neigh[i]);
        delta_energy+=vtx->neigh[i]->xk*(vtx->neigh[i]->energy-oenergy);
    }
//   fprintf(stderr, "DE=%f\n",delta_energy);
    //MONTE CARLOOOOOOOO
    if(delta_energy>=0){
#ifdef TS_DOUBLE_DOUBLE
        if(exp(-delta_energy)< drand48() )
#endif
#ifdef TS_DOUBLE_FLOAT
        if(expf(-delta_energy)< (ts_float)drand48())
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
        if(expl(-delta_energy)< (ts_ldouble)drand48())
#endif
    {
    //not accepted, reverting changes
	vtx=memcpy((void *)vtx,(void *)&backupvtx[0],sizeof(ts_vertex));
	for(i=0;i<vtx->neigh_no;i++){
		vtx->neigh[i]=memcpy((void *)vtx->neigh[i],(void *)&backupvtx[i+1],sizeof(ts_vertex));
	}
	
    //update the normals of triangles that share bead i.
   for(i=0;i<vtx->tristar_no;i++) triangle_normal_vector(vtx->tristar[i]);

    return TS_FAIL; 
    }
}
		
//	oldcellidx=vertex_self_avoidance(vesicle, &backupvtx[0]);
	if(vtx->cell!=vesicle->clist->cell[cellidx]){
		retval=cell_add_vertex(vesicle->clist->cell[cellidx],vtx);
//		if(retval==TS_SUCCESS) cell_remove_vertex(vesicle->clist->cell[oldcellidx],vtx);
		if(retval==TS_SUCCESS) cell_remove_vertex(backupvtx[0].cell,vtx);
		
	}
//	if(oldcellidx);
    //END MONTE CARLOOOOOOO
    return TS_SUCCESS;
}

