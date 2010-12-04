#include<general.h>
#include "vesicle.h"
#include "vertex.h"
#include "triangle.h"
#include "bond.h"
#include "cell.h"


ts_vesicle *init_vesicle(ts_uint N, ts_uint ncmax1, ts_uint ncmax2, ts_uint
ncmax3, ts_double stepsize){
    ts_vesicle *vesicle;
    vesicle->vlist=init_vertex_list(N);
    vesicle->blist=init_bond_list();
    vesicle->tlist=init_triangle_list();
    vesicle->clist=init_cell_list(ncmax1, ncmax2, ncmax3, stepsize);
    return TS_SUCCESS;
}

ts_bool vesicle_translate(ts_vesicle *vesicle,ts_double x, ts_double y, ts_double z){
	ts_uint i;
	ts_vertex *vtx=vesicle->vlist->vertex;
	ts_uint nn=vesicle->vlist->n;
	for(i=0;i<nn;i++){
		vtx[i]->data->x+=x;
		vtx[i]->data->y+=y;
		vtx[i]->data->z+=z;
	}
	return TS_SUCCESS;
}

ts_bool vesicle_free(ts_vesicle *vesicle){
    vertex_list_free(vesicle->vlist);
    bond_list_free(vesicle->blist);
    triangle_list_free(vesicle->tlist);
    cell_list_free(vesicle->clist);
    return TS_SUCCESS;
}
