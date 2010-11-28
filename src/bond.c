#include<stdlib.h>
#include "general.h"
#include "vertex.h"

ts_bond_list *init_bond_list(){
    ts_bond_list *blist=(ts_bond_list *)malloc(sizeof(ts_bond_list));
	blist->n=0;
	blist->bond=NULL;
	return blist;
}

ts_bond  *bond_add(ts_bond_list *blist, ts_vertex *vtx1, ts_vertex *vtx2){
    
    /* no vertices must be null! */
    if(vtx1==NULL || vtx2==NULL) return NULL;
    /* TODO: Verify if the bond already exists... Don't do multiple bonds */
	blist->n++;
	blist->bond=(ts_bond **)realloc(blist->bond,blist->n*sizeof(ts_bond *));
	if(blist->bond==NULL) fatal("Cannot reallocate memory for additional **ts_bond.",100);
    blist->bond[blist->n-1]=(ts_bond *)malloc(sizeof(ts_bond));
    if(blist->bond[blist->n-1]==NULL) fatal("Cannot allocate memory for additional *ts_bond.",100);
    blist->bond[blist->n-1]->data=(ts_bond_data *)malloc(sizeof(ts_bond_data));
    
	//NOW insert vertices into data!	
	blist->bond[blist->n - 1]->data->vtx1=vtx1;	
	blist->bond[blist->n - 1]->data->vtx2=vtx2;

    //Should we calculate bond length NOW?
	
	return blist->bond[blist->n-1];
}

ts_bool bond_list_free(ts_bond_list *blist){
    ts_uint i;
    for(i=0;i<blist->n;i++){
    free(blist->bond[i]->data);
    free(blist->bond[i]);
    }
    free(blist->bond);
    free(blist);
    return TS_SUCCESS;
}
