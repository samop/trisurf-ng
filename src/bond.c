#include<stdlib.h>
#include "general.h"
#include<stdio.h>

ts_bool init_bond_list(ts_bond_list *blist){
	blist->n=0;
	blist->bond=NULL;
	return TS_SUCCESS;
}

ts_bool bond_add(ts_bond_list *blist, ts_vertex *vtx1, ts_vertex *vtx2){
	blist->n++;
	blist->bond=realloc(blist->bond,blist->n*sizeof(ts_bond));
	if(blist->bond==NULL) fatal("Cannot reallocate memory for additional *ts_bond.",5);
	//NOW insert vertices!	
	blist->bond[blist->n - 1].vtx1=vtx1;	
	blist->bond[blist->n - 1].vtx2=vtx2;	
	return TS_SUCCESS;
}

ts_bool bond_list_free(ts_bond_list *blist){
    free(blist->bond);
    return TS_SUCCESS;
}
