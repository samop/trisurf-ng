#include<stdlib.h>
#include<stdio.h>
#include<math.h>
//#include "io.h"
#include "general.h"
#include "timestep.h"
#include "vertexmove.h"
#include "bondflip.h"

ts_bool single_timestep(ts_vesicle *vesicle){
    ts_bool retval;
    ts_double rnvec[3];
    ts_uint i;
    for(i=0;i<vesicle->vlist->n;i++){
        rnvec[0]=drand48();
        rnvec[1]=drand48();
        rnvec[2]=drand48();
        retval=single_verticle_timestep(vesicle,vesicle->vlist->vtx[i],rnvec);
    }

//	ts_int cnt=0;
    for(i=0;i<vesicle->blist->n;i++){
        rnvec[0]=drand48();
        rnvec[1]=drand48();
        rnvec[2]=drand48();
        //find a bond and return a pointer to a bond...
        //call single_bondflip_timestep...
        retval=single_bondflip_timestep(vesicle,vesicle->blist->bond[i],rnvec);
//	if(retval==TS_SUCCESS) cnt++;        
    } 
//	printf("Bondflip success rate in one sweep: %d/%d=%e\n", cnt,vesicle->blist->n,(double)cnt/(double)vesicle->blist->n);
	if(retval);
    return TS_SUCCESS;
}



