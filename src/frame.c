#include<stdlib.h>
#include "general.h"
#include "cell.h"

ts_bool centermass(ts_vesicle *vesicle){
    ts_uint i;
    vesicle->cm[0]=0;
    vesicle->cm[1]=0;
    vesicle->cm[2]=0;
    for(i=0;i<vesicle->vlist.n;i++){
        vesicle->cm[0]+=vesicle->vlist.vertex[i].x;
        vesicle->cm[1]+=vesicle->vlist.vertex[i].y;
        vesicle->cm[2]+=vesicle->vlist.vertex[i].z; 
    } 
    vesicle->cm[0]/=(float)vesicle->vlist.n;
    vesicle->cm[1]/=(float)vesicle->vlist.n;
    vesicle->cm[2]/=(float)vesicle->vlist.n;
    return TS_SUCCESS;
}

ts_bool cell_ocupation(ts_vesicle *vesicle){
    ts_uint i,j,cellidx;
    ts_double shift;
    ts_double dcell;
    shift=(ts_double) vesicle->clist.ncmax[0]/2;
    dcell=1.0/(1.0 + vesicle->stepsize);
    ts_uint ncx, ncy,ncz;

    cell_list_cell_ocupation_clear(&vesicle->clist);
    for(i=0;i<vesicle->vlist.n;i++){
  
    cellidx=vertex_self_avoidance(vesicle, &vesicle->vlist.vertex[i]);
    vesicle->vlist.vertex[i].cell=&(vesicle->clist.cell[cellidx]);

    cell_add_vertex(&vesicle->clist.cell[cellidx],&vesicle->vlist.vertex[i]);

  //  if(ncx > vesicle->clist.ncmax[0]) vesicle->clist.ncmax[0]=ncx;
  //  if(ncy > vesicle->clist.ncmax[1]) vesicle->clist.ncmax[1]=ncy;
  //  if(ncz > vesicle->clist.ncmax[2]) vesicle->clist.ncmax[2]=ncz;
    }


/* This was already done in previous for loop.... Have I gained some time? 


    for(i=0;i<vesicle->clist.ncmax[0]*vesicle->clist.ncmax[1]*vesicle->clist.ncmax[2];i++){
        vesicle->clist.cell[i].nvertex=0;
        for(j=0;j<vesicle->vlist.n;j++){
            //add_vertextocell;
            //add_vertextomonomer;
        }
    }

*/
    return TS_SUCCESS;
}
