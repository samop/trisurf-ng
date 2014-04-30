#include<stdlib.h>
#include "general.h"
#include "cell.h"
#include "frame.h"


#include "triangle.h"
ts_bool centermass(ts_vesicle *vesicle){
    ts_uint i,j, n=vesicle->vlist->n;
    ts_vertex **vtx=vesicle->vlist->vtx;
    vesicle->cm[0]=0;
    vesicle->cm[1]=0;
    vesicle->cm[2]=0;
    for(i=0;i<n;i++){
        vesicle->cm[0]+=vtx[i]->x;
        vesicle->cm[1]+=vtx[i]->y;
        vesicle->cm[2]+=vtx[i]->z; 
    } 
    vesicle->cm[0]/=(ts_float)n;
    vesicle->cm[1]/=(ts_float)n;
    vesicle->cm[2]/=(ts_float)n;

    for(i=0;i<n;i++){
        vtx[i]->x-=vesicle->cm[0];
        vtx[i]->y-=vesicle->cm[1];
        vtx[i]->z-=vesicle->cm[2]; 
    } 
//move polymers for the same vector as we moved vesicle
	for(i=0;i<vesicle->poly_list->n;i++){
		for(j=0;j<vesicle->poly_list->poly[i]->vlist->n;j++){
			vesicle->poly_list->poly[i]->vlist->vtx[j]->x-=vesicle->cm[0];
			vesicle->poly_list->poly[i]->vlist->vtx[j]->y-=vesicle->cm[1];
			vesicle->poly_list->poly[i]->vlist->vtx[j]->z-=vesicle->cm[2];
		}
    }
//move filaments for the same vector as we moved vesicle
	for(i=0;i<vesicle->filament_list->n;i++){
		for(j=0;j<vesicle->filament_list->poly[i]->vlist->n;j++){
			vesicle->filament_list->poly[i]->vlist->vtx[j]->x-=vesicle->cm[0];
			vesicle->filament_list->poly[i]->vlist->vtx[j]->y-=vesicle->cm[1];
			vesicle->filament_list->poly[i]->vlist->vtx[j]->z-=vesicle->cm[2];
		}
    }

    vesicle->cm[0]=0.0;
    vesicle->cm[1]=0.0;
    vesicle->cm[2]=0.0;

    for(i=0;i<vesicle->tlist->n;i++){
        triangle_normal_vector(vesicle->tlist->tria[i]);
    }


    return TS_SUCCESS;
}

ts_bool cell_occupation(ts_vesicle *vesicle){
    ts_uint i,j,cellidx, n=vesicle->vlist->n;

    cell_list_cell_occupation_clear(vesicle->clist);
    for(i=0;i<n;i++){
    cellidx=vertex_self_avoidance(vesicle, vesicle->vlist->vtx[i]);
//	already done in cell_add_vertex
//    vesicle->vlist->vtx[i]->cell=vesicle->clist->cell[cellidx];

    cell_add_vertex(vesicle->clist->cell[cellidx],vesicle->vlist->vtx[i]);
    }

//Add all polymers to cells
    for(i=0;i<vesicle->poly_list->n;i++){
	for(j=0;j<vesicle->poly_list->poly[i]->vlist->n;j++){
    	cellidx=vertex_self_avoidance(vesicle, vesicle->poly_list->poly[i]->vlist->vtx[j]);
    	cell_add_vertex(vesicle->clist->cell[cellidx],vesicle->poly_list->poly[i]->vlist->vtx[j]);
	}
    }
//Add all filaments to cells
     for(i=0;i<vesicle->filament_list->n;i++){
	for(j=0;j<vesicle->filament_list->poly[i]->vlist->n;j++){
    	cellidx=vertex_self_avoidance(vesicle, vesicle->filament_list->poly[i]->vlist->vtx[j]);
    	cell_add_vertex(vesicle->clist->cell[cellidx],vesicle->filament_list->poly[i]->vlist->vtx[j]);
	}
    }
   

    return TS_SUCCESS;
}
