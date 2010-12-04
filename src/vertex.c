#include<stdlib.h>
#include<math.h>
#include<string.h>
#include "general.h"
#include "vertex.h"
#include "bond.h"
#include<stdio.h>

ts_vertex_list *init_vertex_list(ts_uint N){	
	ts_int i;
    ts_vertex *tlist;
    ts_vertex_list *vlist=(ts_vertex_list *)malloc(sizeof(ts_vertex_list));
    
	if(N==0){
		err("Initialized vertex list with zero elements. Pointer set to NULL");
        vlist->n=0;
		vlist->vtx=NULL;
		return vlist;
	}
	
    vlist->vtx=(ts_vertex **)malloc(N*sizeof(ts_vertex *));
    tlist=(ts_vertex *)malloc(N*sizeof(ts_vertex));
    if(vlist->vtx==NULL || tlist==NULL)
        fatal("Fatal error reserving memory space for vertex list! Could number of requsted vertices be too large?", 100);
    for(i=0;i<N;i++) {
        vlist->vtx[i]=&tlist[i];
        vlist->vtx[i]->data=init_vertex_data();
        vlist->vtx[i]->idx=i;
    }
    vlist->n=N;
	return vlist;
}

ts_vertex_data *init_vertex_data(){
    ts_vertex_data *data;
    data=(ts_vertex_data *)calloc(1,sizeof(ts_vertex_data));
    if(data==NULL)
        fatal("Fatal error reserving memory space for ts_vertex! Memory full?", 100);
    return data;
}

/*
ts_bool vtx_set_global_values(ts_vertex **vlist, ts_vesicle *vesicle){
    ts_double xk=vesicle->bending_rigidity;
    ts_uint i;
    for(i=0;i<vesicle->vlist.n;i++){
        vlist[i]->xk=xk;
    }
    return TS_SUCCESS;
}
*/



ts_bool vtx_add_neighbour(ts_vertex *vtx, ts_vertex *nvtx){
    ts_uint i;
    /* no neighbour can be null! */
    if(vtx==NULL || nvtx==NULL) return TS_FAIL;
    
    /*if it is already a neighbour don't add it to the list */
    for(i=0; i<vtx->data->neigh_no;i++){
        if(vtx->data->neigh[i]==nvtx) return TS_FAIL;
    }
    ts_uint nn=++vtx->data->neigh_no;
    vtx->data->neigh=(ts_vertex **)realloc(vtx->data->neigh, nn*sizeof(ts_vertex *));
    vtx->data->neigh[nn-1]=nvtx;

    /* pa se sosedu dodamo vertex */
    /*if it is already a neighbour don't add it to the list */
    for(i=0; i<nvtx->data->neigh_no;i++){
        if(nvtx->data->neigh[i]==vtx) return TS_FAIL;
    } 
    nn=++nvtx->data->neigh_no;
    nvtx->data->neigh=(ts_vertex **)realloc(nvtx->data->neigh, nn*sizeof(ts_vertex *));
    nvtx->data->neigh[nn-1]=vtx;


    return TS_SUCCESS;
}

/* TODO: optimize this. test this. */
ts_bool vtx_remove_neighbour(ts_vertex *vtx, ts_vertex *nvtx){
/* find a neighbour */
/* remove it from the list while shifting remaining neighbours up */
    ts_uint i,j=0;
    for(i=0;i<vtx->data->neigh_no;i++){
        if(vtx->data->neigh[i]!=nvtx){
            vtx->data->neigh[j]=vtx->data->neigh[i];
            j++;
        }
    }
/* resize memory. potentionally time consuming */
    vtx->data->neigh_no--;
    vtx->data->neigh=(ts_vertex **)realloc(vtx->data->neigh,vtx->data->neigh_no*sizeof(ts_vertex *));
    if(vtx->data->neigh == NULL && vtx->data->neigh_no!=0)
        fatal("Reallocation of memory failed during removal of vertex neighbour in vtx_remove_neighbour",100);

/* repeat for the neighbour */
/* find a neighbour */
/* remove it from the list while shifting remaining neighbours up */
    for(i=0;i<nvtx->data->neigh_no;i++){
        if(nvtx->data->neigh[i]!=vtx){
            nvtx->data->neigh[j]=nvtx->data->neigh[i];
            j++;
        }
    }
/* resize memory. potentionally time consuming. */
    nvtx->data->neigh_no--;
    nvtx->data->neigh=(ts_vertex **)realloc(nvtx->data->neigh,nvtx->data->neigh_no*sizeof(ts_vertex *));
    if(nvtx->data->neigh == NULL && nvtx->data->neigh_no!=0)
        fatal("Reallocation of memory failed during removal of vertex neighbour in vtx_remove_neighbour",100);

    return TS_SUCCESS;
}


ts_bool vtx_add_bond(ts_bond_list *blist,ts_vertex *vtx1,ts_vertex *vtx2){
    ts_bond *bond;
    bond=bond_add(blist,vtx1,vtx2);
    if(bond==NULL) return TS_FAIL;
    vtx1->data->bond_no++;
    vtx2->data->bond_no++;

    vtx1->data->bond=(ts_bond **)realloc(vtx1->data->bond, vtx1->data->bond_no*sizeof(ts_bond *)); 
    vtx2->data->bond=(ts_bond **)realloc(vtx2->data->bond, vtx2->data->bond_no*sizeof(ts_bond *)); 
    vtx1->data->bond[vtx1->data->bond_no-1]=bond;
    vtx2->data->bond[vtx2->data->bond_no-1]=bond;
    return TS_SUCCESS;
}

ts_bool vtx_add_cneighbour(ts_bond_list *blist, ts_vertex *vtx1, ts_vertex *vtx2){
    ts_bool retval;
    retval=vtx_add_neighbour(vtx1,vtx2);
    if(retval==TS_SUCCESS)
    retval=vtx_add_bond(blist,vtx1,vtx2); 
    return retval;
}

/*TODO: write and optimize this urgently before use! */
ts_bool vtx_remove_cneighbour(ts_bond_list *blist, ts_vertex *vtx1, ts_vertex
*vtx2){
//    ts_bool retval;
/* remove the bond */
//retval=vtx_remove_bond(blist,vtx1,vtx2);
/* remove the vertices */
    return TS_SUCCESS;
}



ts_bool vtx_data_free(ts_vertex_data *data){
    if(data->neigh!=NULL)   free(data->neigh);
    if(data->tristar!=NULL) free(data->tristar);
    if(data->bond!=NULL)    free(data->bond);
    if(data->cell!=NULL)    free(data->cell);
    free(data);
    return TS_SUCCESS;
}

/*not usable. can be deleted */
ts_bool vtx_free(ts_vertex  *vtx){
    vtx_data_free(vtx->data);
    free(vtx);
    return TS_SUCCESS;
}

ts_bool vtx_list_free(ts_vertex_list *vlist){
    int i;
    for(i=0;i<vlist->n;i++){
        vtx_data_free(vlist->vtx[i]->data);
    }
    free(*(vlist->vtx));
    free(vlist->vtx);
    free(vlist);
    return TS_SUCCESS;
}

inline ts_double vtx_distance_sq(ts_vertex *vtx1, ts_vertex *vtx2){
    ts_double dist;
    ts_vertex_data *vd1=vtx1->data, *vd2=vtx2->data;
#ifdef TS_DOUBLE_DOUBLE
    dist=pow(vd1->x-vd2->x,2) + pow(vd1->y-vd2->y,2) + pow(vd1->z-vd2->z,2);
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
    dist=powl(vd1->x-vd2->x,2) + powl(vd1->y-vd2->y,2) + powl(vd1->z-vd2->z,2);
#endif
#ifdef TS_DOUBLE_FLOAT
    dist=powf(vd1->x-vd2->x,2) + powf(vd1->y-vd2->y,2) + powf(vd1->z-vd2->z,2);
#endif
    return(dist);
}

