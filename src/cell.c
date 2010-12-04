#include<stdlib.h>
#include "general.h"
#include "vertex.h"

ts_cell_list  *init_cell_list(ts_uint ncmax1, ts_uint ncmax2, ts_uint ncmax3, ts_double stepsize){
    ts_uint i;
    ts_uint nocells=ncmax1*ncmax2*ncmax3;
    ts_cell_list *clist=(ts_cell_list *)malloc(sizeof(ts_cell_list));
    if(clist==NULL) fatal("Error while allocating memory for cell list!",100);

    clist->ncmax[0]=ncmax1;
    clist->ncmax[1]=ncmax2;
    clist->ncmax[2]=ncmax3;
    clist->cellno=nocells;
    clist->dcell=1.0/(1.0 + stepsize);
    clist->shift=(ts_double) clist->ncmax[0]/2;

    clist->cell=(ts_cell **)malloc(nocells*sizeof(ts_cell *));
    if(clist->cell==NULL) fatal("Error while allocating memory for cell list! ncmax too large?",101);

    for(i=0;i<nocells;i++){
        clist->cell[i]=(ts_cell *)malloc(sizeof(ts_cell));
        if(clist->cell[i]==NULL) fatal("Error while allocating memory for cell list! ncmax too large?",102);
        clist->cell[i]->idx=i+1; // We enumerate cells! Probably never required!
        clist->cell[i]->data=(ts_cell_data *)calloc(1,sizeof(ts_cell_data));
    }
    return clist;
}

ts_bool cell_free(ts_cell* cell){
    if(cell->data!=NULL){
        if(cell->data->vertex!=NULL) free(cell->data->vertex);
        free(cell->data);
    }
    free(cell);
    return TS_SUCCESS;
}

ts_bool cell_list_free(ts_cell_list *clist){
    ts_uint i;
    if(clist==NULL) return TS_FAIL;
    ts_uint nocells=clist->cellno;
    for(i=0;i<nocells;i++)
         if(clist->cell[i] != NULL) cell_free(clist->cell[i]);
    free(clist->cell);
    free(clist);
    return TS_SUCCESS;
}


//TODO: not debugged at all!
inline ts_uint vertex_self_avoidance(ts_vesicle *vesicle, ts_vertex *vtx){
    ts_uint cellidx;
    ts_uint ncx, ncy,ncz;
    ts_cell_list *clist=vesicle->clist;
    ncx=(ts_uint)((vtx->data->x-vesicle->cm[0])*clist->dcell+clist->shift);
    ncy=(ts_uint)((vtx->data->y-vesicle->cm[1])*clist->dcell+clist->shift);
    ncz=(ts_uint)((vtx->data->z-vesicle->cm[2])*clist->dcell+clist->shift);

    if(ncx == clist->ncmax[0]-1 || ncx == 2){
        fatal("Vesicle is positioned outside the cell covered area. Coordinate x is the problem.",1500);
    }
    if(ncy == clist->ncmax[1]-1 || ncy == 2){
        fatal("Vesicle is positioned outside the cell covered area. Coordinate y is the problem.",1500);
    }
    if(ncz == clist->ncmax[2]-1 || ncz == 2){
        fatal("Vesicle is positioned outside the cell covered area. Coordinate z is the problem.",1500);
    }
    cellidx=ncz+(ncy-1)*clist->ncmax[2] + (ncx-1)*clist->ncmax[2]* 
                                    clist->ncmax[1] - 1; // -1 is because of 0 based indexing
        return cellidx;
}


//TODO: looks ok, but debug anyway in the future
ts_bool cell_add_vertex(ts_cell *cell, ts_vertex *vtx){
    cell->data->nvertex++;
	cell->data->vertex=(ts_vertex **)realloc(cell->data->vertex,cell->data->nvertex*sizeof(ts_vertex *));
		if(vtx->data->neigh == NULL){
			fatal("Reallocation of memory failed during insertion of vertex neighbour in vertex_add_neighbour",3);
        }
    cell->data->vertex[cell->data->nvertex-1]=vtx;
    return TS_SUCCESS;
}


ts_bool cell_list_cell_occupation_clear(ts_cell_list *clist){
    ts_uint i;
    for(i=0;i<clist->cellno;i++){
        if(clist->cell[i]->data->vertex != NULL){
            free(clist->cell[i]->data->vertex);
            clist->cell[i]->data->vertex=NULL;
        }
        clist->cell[i]->data->nvertex=0;
    }
    return TS_SUCCESS;
}

// TODO: compiles ok, but it is completely untested and undebugged. It was debugged before rewrite, but this was long time ago.
ts_bool cell_occupation_number_and_internal_proximity(ts_cell_list *clist, ts_uint cellidx, ts_vertex *vtx, ts_vertex *tvtx){
    ts_uint ncx,ncy,ncz,remainder,cell_occupation;
    ts_uint i,j,k,l,neigh_cidx;
    ts_double dist;
    ncx=(cellidx+1)/(clist->ncmax[2]*clist->ncmax[1])+1; //+1 because of zero indexing.
    remainder=(cellidx+1)%(clist->ncmax[2]*clist->ncmax[1]);
    ncy=remainder/clist->ncmax[2]+1;
    ncz=remainder%clist->ncmax[2];
//    fprintf(stderr,"here are ncx,ncy,ncz=%i,%i,%i\n",ncx,ncy,ncz);

    for(i=ncx-1;i<=ncx+1;i++){
        for(j=ncy-1;j<=ncy+1;j++){
            for(k=ncz-1;k<=ncz+1;k++){
                neigh_cidx=k+(j-1)*clist->ncmax[2]+(i-1)*clist->ncmax[2]*clist->ncmax[1] -1;
          //      fprintf(stderr,"neigh_cell_index=%i\n",neigh_cidx);
                cell_occupation=clist->cell[neigh_cidx]->data->nvertex;
          //      fprintf(stderr, "cell_occupation=%i\n",cell_occupation);
                if(cell_occupation>clist->max_occupancy){
                    fatal("Neighbouring cell occupation more than set max_occupancy value.",2000);
                }
// Now we check whether we didn't come close to some other vertices in the same
// cell!
                if(cell_occupation>1){
                    for(l=0;l<cell_occupation;l++){
                        if(clist->cell[neigh_cidx]->data->vertex[l]!=vtx){
                    //        fprintf(stderr,"calling dist on vertex %i\n",l);
                           dist=vtx_distance_sq(clist->cell[neigh_cidx]->data->vertex[l],tvtx);
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
