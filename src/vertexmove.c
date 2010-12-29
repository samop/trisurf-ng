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
    ts_vertex *tvtx=(ts_vertex *)malloc(sizeof(ts_vertex));
	tvtx->data=init_vertex_data();
    ts_bool retval; 
    ts_uint cellidx; 
    ts_double xold,yold,zold;
    ts_double delta_energy,oenergy;
    ts_vertex *ovtx;

    //randomly we move the temporary vertex
    tvtx->data->x=vtx->data->x+vesicle->stepsize*(2.0*rn[0]-1.0);
    tvtx->data->y=vtx->data->y+vesicle->stepsize*(2.0*rn[1]-1.0);
    tvtx->data->z=vtx->data->z+vesicle->stepsize*(2.0*rn[2]-1.0);
    //check we if some length to neighbours are too much
    for(i=0;i<vtx->data->neigh_no;i++){
        dist=vtx_distance_sq(tvtx,vtx->data->neigh[i]);
        if(dist<1.0 || dist>vesicle->dmax) return TS_FAIL;
    }
fprintf(stderr,"Was here!\n");
    //self avoidance check with distant vertices
     cellidx=vertex_self_avoidance(vesicle, tvtx);
    //check occupation number
     retval=cell_occupation_number_and_internal_proximity(vesicle->clist,cellidx,vtx,tvtx);
    if(retval==TS_FAIL){
        return TS_FAIL;
    } 
    
    //if all the tests are successful, then we update the vertex position
    xold=vtx->data->x;
    yold=vtx->data->y;
    zold=vtx->data->z;
    ovtx=malloc(sizeof(ts_vertex));
    vtx_copy(ovtx,vtx);
    vtx->data->x=tvtx->data->x;
    vtx->data->y=tvtx->data->y;
    vtx->data->z=tvtx->data->z;

    delta_energy=0;
    //update the normals of triangles that share bead i.
    for(i=0;i<vtx->data->tristar_no;i++) triangle_normal_vector(vtx->data->tristar[i]);
    //energy and curvature
    energy_vertex(vtx);
    delta_energy=vtx->data->xk*(vtx->data->energy - ovtx->data->energy);
    //the same is done for neighbouring vertices
    for(i=0;i<vtx->data->neigh_no;i++){
        oenergy=vtx->data->neigh[i]->data->energy;
        energy_vertex(vtx->data->neigh[i]);
        delta_energy+=vtx->data->neigh[i]->data->xk*(vtx->data->neigh[i]->data->energy-oenergy);
    }
    fprintf(stderr, "DE=%f\n",delta_energy);
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
    vtx->data->x=xold;
    vtx->data->y=yold;
    vtx->data->z=zold;
    //update the normals of triangles that share bead i.
    for(i=0;i<vtx->data->tristar_no;i++) triangle_normal_vector(vtx->data->tristar[i]);
    //energy and curvature
    energy_vertex(vtx);
    //the same is done for neighbouring vertices
    for(i=0;i<vtx->data->neigh_no;i++) energy_vertex(vtx->data->neigh[i]);
  free(ovtx->data->bond_length);
    free(ovtx->data->bond_length_dual);
    free(ovtx);
    return TS_FAIL; 
    }
}
    //END MONTE CARLOOOOOOO

    //TODO: change cell occupation if necessary!

    free(ovtx->data->bond_length);
    free(ovtx->data->bond_length_dual);
    free(ovtx);
    return TS_SUCCESS;
}

