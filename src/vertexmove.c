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

ts_bool single_verticle_timestep(ts_vesicle *vesicle,ts_vertex *vtx,ts_double
*rn){
    ts_uint i;
    ts_double dist;
    ts_bool retval; 
    ts_uint cellidx; 
    ts_double xold,yold,zold;
    ts_double delta_energy,oenergy;
    ts_vertex *ovtx;
    ts_vertex *tvtx=(ts_vertex *)calloc(1,sizeof(ts_vertex));
    ts_double costheta,sintheta,phi,r;

    //randomly we move the temporary vertex
//    tvtx->x=vtx->x+vesicle->stepsize*(2.0*rn[0]-1.0);
//    tvtx->y=vtx->y+vesicle->stepsize*(2.0*rn[1]-1.0);
//    tvtx->z=vtx->z+vesicle->stepsize*(2.0*rn[2]-1.0);
	//random move in a sphere with radius stepsize:
	r=vesicle->stepsize*rn[0];
	phi=rn[1]*2*M_PI;
	costheta=2*rn[2]-1;
	sintheta=sqrt(1-pow(costheta,2));
	tvtx->x=vtx->x+r*sintheta*cos(phi);
	tvtx->y=vtx->y+r*sintheta*sin(phi);
	tvtx->z=vtx->z+r*costheta;
    //check we if some length to neighbours are too much
    for(i=0;i<vtx->neigh_no;i++){
        dist=vtx_distance_sq(tvtx,vtx->neigh[i]);
        if(dist<1.0 || dist>vesicle->dmax) {
		vtx_free(tvtx);
//	fprintf(stderr,"Fail 1, dist=%f, vesicle->dmax=%f\n", dist, vesicle->dmax);
		return TS_FAIL;
		}
    }
    //self avoidance check with distant vertices
     cellidx=vertex_self_avoidance(vesicle, tvtx);
    //check occupation number
     retval=cell_occupation_number_and_internal_proximity(vesicle->clist,cellidx,vtx,tvtx);
    if(retval==TS_FAIL){
	vtx_free(tvtx);
//	fprintf(stderr,"Fail 2\n");
        return TS_FAIL;
    } 
    
    //if all the tests are successful, then we update the vertex position
    xold=vtx->x;
    yold=vtx->y;
    zold=vtx->z;
    ovtx=malloc(sizeof(ts_vertex));
    vtx_copy(ovtx,vtx);
    vtx->x=tvtx->x;
    vtx->y=tvtx->y;
    vtx->z=tvtx->z;

    delta_energy=0;
    //update the normals of triangles that share bead i.
    for(i=0;i<vtx->tristar_no;i++) triangle_normal_vector(vtx->tristar[i]);
    //energy and curvature
    energy_vertex(vtx);
    delta_energy=vtx->xk*(vtx->energy - ovtx->energy);
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
    vtx->x=xold;
    vtx->y=yold;
    vtx->z=zold;
    //update the normals of triangles that share bead i.
    for(i=0;i<vtx->tristar_no;i++) triangle_normal_vector(vtx->tristar[i]);
    //energy and curvature
    energy_vertex(vtx);
    //the same is done for neighbouring vertices
	for(i=0;i<vtx->neigh_no;i++) energy_vertex(vtx->neigh[i]);
	free(ovtx->bond_length);
    free(ovtx->bond_length_dual);
    free(ovtx);
    vtx_free(tvtx);
    return TS_FAIL; 
    }
}
    //END MONTE CARLOOOOOOO

    //TODO: change cell occupation if necessary!
//	fprintf(stderr,"Success!!\n");
    free(ovtx->bond_length);
    free(ovtx->bond_length_dual);
    free(ovtx);
    vtx_free(tvtx);
    return TS_SUCCESS;
}

