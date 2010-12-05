#include<stdlib.h>
#include "general.h"
#include "cell.h"
#include "frame.h"
ts_bool centermass(ts_vesicle *vesicle){
    ts_uint i, n=vesicle->vlist->n;
    ts_vertex **vtx=vesicle->vlist->vtx;
    vesicle->cm[0]=0;
    vesicle->cm[1]=0;
    vesicle->cm[2]=0;
    for(i=0;i<n;i++){
        vesicle->cm[0]+=vtx[i]->data->x;
        vesicle->cm[1]+=vtx[i]->data->y;
        vesicle->cm[2]+=vtx[i]->data->z; 
    } 
    vesicle->cm[0]/=(float)n;
    vesicle->cm[1]/=(float)n;
    vesicle->cm[2]/=(float)n;
    return TS_SUCCESS;
}

ts_bool cell_occupation(ts_vesicle *vesicle){
    ts_uint i,cellidx, n=vesicle->vlist->n;
    ts_double shift;
    ts_double dcell;
    shift=(ts_double) vesicle->clist->ncmax[0]/2;
    dcell=1.0/(1.0 + vesicle->stepsize);

    cell_list_cell_occupation_clear(vesicle->clist);
    for(i=0;i<n;i++){
  
    cellidx=vertex_self_avoidance(vesicle, vesicle->vlist->vtx[i]);
    vesicle->vlist->vtx[i]->data->cell=vesicle->clist->cell[cellidx];

    cell_add_vertex(vesicle->clist->cell[cellidx],vesicle->vlist->vtx[i]);

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
