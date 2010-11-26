#include<stdlib.h>
#include "general.h"
#include "vertex.h"
#include <stdio.h>
ts_bool init_cell_list(ts_cell_list *clist, ts_double stepsize){
    ts_uint i;
    ts_uint nocells=clist->ncmax[0]*clist->ncmax[1]*clist->ncmax[2];
    clist->cell=malloc(nocells*sizeof(ts_cell));
    clist->dcell=1.0/(1.0 + stepsize);
    clist->shift=(ts_double) clist->ncmax[0]/2;
    clist->cellno=nocells;
    for(i=0;i<nocells;i++){
        clist->cell[i].idx=i+1; // We enumerate cells! Probably never required!
        clist->cell[i].nvertex=0;
        clist->cell[i].vertex=NULL;
    }
    return TS_SUCCESS;
}

ts_bool cell_list_free(ts_cell_list *clist){
    ts_uint i;
     ts_uint nocells=clist->ncmax[0]*clist->ncmax[1]*clist->ncmax[2];

    for(i=0;i<nocells;i++)
         if(clist->cell->vertex != NULL) free(clist->cell->vertex);
    free(clist->cell);

    return TS_SUCCESS;
}

inline ts_uint vertex_self_avoidance(ts_vesicle *vesicle, ts_vertex *vtx){
    ts_uint i,cellidx;
    ts_uint ncx, ncy,ncz;

    ncx=(ts_uint)((vtx->x-vesicle->cm[0])*vesicle->clist.dcell+vesicle->clist.shift);
    ncy=(ts_uint)((vtx->y-vesicle->cm[1])*vesicle->clist.dcell+vesicle->clist.shift);
    ncz=(ts_uint)((vtx->z-vesicle->cm[2])*vesicle->clist.dcell+vesicle->clist.shift);
//    fprintf(stderr,"(ncx,ncy,ncz)=(%i %i %i)\t",ncx,ncy,ncz);
//    fprintf(stderr,"(ncxmax,ncymax,nczmax)=(%i %i %i)\n",vesicle->clist.ncmax[0], vesicle->clist.ncmax[1], vesicle->clist.ncmax[2]);
    if(ncx == vesicle->clist.ncmax[0]-1 || ncx == 2){
        fatal("Vesicle is positioned outside the cell covered area. Coordinate x is the problem.",1500);
    }
    if(ncy == vesicle->clist.ncmax[1]-1 || ncy == 2){
        fatal("Vesicle is positioned outside the cell covered area. Coordinate y is the problem.",1500);
    }
    if(ncz == vesicle->clist.ncmax[2]-1 || ncz == 2){
        fatal("Vesicle is positioned outside the cell covered area. Coordinate z is the problem.",1500);
    }
    cellidx=ncz+(ncy-1)*vesicle->clist.ncmax[2] +
(ncx-1)*vesicle->clist.ncmax[2]* vesicle->clist.ncmax[1] - 1; // -1 is because of 0 based indexing
    return cellidx;
}



ts_bool cell_add_vertex(ts_cell *cell, ts_vertex *vtx){
    
    cell->nvertex++;
	cell->vertex=realloc(cell->vertex,cell->nvertex*sizeof(ts_vertex *));
		if(vtx->neigh == NULL){
			fatal("Reallocation of memory failed during insertion of vertex neighbour in vertex_add_neighbour",3);
        }
    cell->vertex[cell->nvertex-1]=vtx;

    return TS_SUCCESS;
}

ts_bool cell_list_cell_ocupation_clear(ts_cell_list *clist){
    ts_uint i;
    for(i=0;i<clist->cellno;i++){
        if(clist->cell[i].vertex != NULL){
            free(clist->cell[i].vertex);
            clist->cell[i].vertex=NULL;
        }
        clist->cell[i].nvertex=0;
    }
    return TS_SUCCESS;
}

ts_bool cell_occupation_number_and_internal_proximity(ts_cell_list *clist, ts_uint cellidx, ts_vertex *vtx, ts_vertex *tvtx){
    ts_uint ncx,ncy,ncz,remainder,cell_occupation;
    ts_uint i,j,k,l,neigh_cidx,mcell;
    ts_double dist;
    ncx=(cellidx+1)/(clist->ncmax[2]*clist->ncmax[1])+1;
    remainder=(cellidx+1)%(clist->ncmax[2]*clist->ncmax[1]);
    ncy=remainder/clist->ncmax[2]+1;
    ncz=remainder%clist->ncmax[2];
//    fprintf(stderr,"here are ncx,ncy,ncz=%i,%i,%i\n",ncx,ncy,ncz);

    for(i=ncx-1;i<=ncx+1;i++){
        for(j=ncy-1;j<=ncy+1;j++){
            for(k=ncz-1;k<=ncz+1;k++){
                neigh_cidx=k+(j-1)*clist->ncmax[2]+(i-1)*clist->ncmax[2]*clist->ncmax[1] -1;
          //      fprintf(stderr,"neigh_cell_index=%i\n",neigh_cidx);
                cell_occupation=clist->cell[neigh_cidx].nvertex;
          //      fprintf(stderr, "cell_occupation=%i\n",cell_occupation);
                if(cell_occupation>clist->max_occupancy){
                    fatal("Neighbouring cell occupation more than set max_occupancy value.",2000);
                }
// Now we check whether we didn't come close to some other vertices in the same
// cell!
                if(cell_occupation>1){
                    for(l=0;l<cell_occupation;l++){
                        if(clist->cell[neigh_cidx].vertex[l]!=vtx){
                    //        fprintf(stderr,"calling dist on vertex %i\n",l);
                           dist=vertex_distance_sq(clist->cell[neigh_cidx].vertex[l],tvtx);
                    //        fprintf(stderr,"dist was %f\n",dist);
                            if(dist<1) return TS_FAIL;
                        }
                    }
                }
            }
        }
    }
    
    return TS_SUCCESS;
}
