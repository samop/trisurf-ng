#include<general.h>
#include "vesicle.h"
ts_bool vesicle_translate(ts_vesicle *vesicle,ts_double x, ts_double y, ts_double z){
	ts_uint i;
	ts_vertex *vtx=vesicle->vlist.vertex;
	ts_uint nn=vesicle->vlist.n;
	for(i=0;i<nn;i++){
		vtx[i].x+=x;
		vtx[i].y+=y;
		vtx[i].z+=z;
	}
	return TS_SUCCESS;
}

ts_bool vesicle_free(ts_vesicle *vesicle){
    vertex_list_free(&vesicle->vlist);
    bond_list_free(&vesicle->blist);
    triangle_list_free(&vesicle->tlist);
    cell_list_free(&vesicle->clist);
    return TS_SUCCESS;
}
