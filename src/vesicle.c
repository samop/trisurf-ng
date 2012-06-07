#include<general.h>
#include "vesicle.h"
#include "vertex.h"
#include "triangle.h"
#include "bond.h"
#include "cell.h"
#include "stdlib.h"

ts_vesicle *init_vesicle(ts_uint N, ts_uint ncmax1, ts_uint ncmax2, ts_uint
ncmax3, ts_double stepsize){
    ts_vesicle *vesicle=(ts_vesicle *)malloc(sizeof(ts_vesicle));
    vesicle->vlist=init_vertex_list(N);
    vesicle->blist=init_bond_list();
    vesicle->tlist=init_triangle_list();
    vesicle->clist=init_cell_list(ncmax1, ncmax2, ncmax3, stepsize);
    return vesicle;
}

ts_bool vesicle_translate(ts_vesicle *vesicle,ts_double x, ts_double y, ts_double z){
	ts_uint i;
	ts_vertex **vtx=vesicle->vlist->vtx;
	ts_uint nn=vesicle->vlist->n;
	for(i=0;i<nn;i++){
		vtx[i]->x+=x;
		vtx[i]->y+=y;
		vtx[i]->z+=z;
	}
	return TS_SUCCESS;
}

ts_bool vesicle_free(ts_vesicle *vesicle){
    vtx_list_free(vesicle->vlist);
    bond_list_free(vesicle->blist);
    triangle_list_free(vesicle->tlist);
    cell_list_free(vesicle->clist);
    free(vesicle);
    return TS_SUCCESS;
}

ts_bool vesicle_volume(ts_vesicle *vesicle){
    ts_double volume;
    ts_double vol;
    ts_uint i;
    ts_triangle **tria=vesicle->tlist->tria;
    volume=0;
    for(i=0; i<vesicle->tlist->n;i++){
        vol=(tria[i]->vertex[0]->x+ tria[i]->vertex[1]->x + tria[i]->vertex[2]->x) * tria[i]->xnorm + 
       (tria[i]->vertex[0]->y+ tria[i]->vertex[1]->y + tria[i]->vertex[2]->y) * tria[i]->ynorm + 
    (tria[i]->vertex[0]->z+ tria[i]->vertex[1]->z + tria[i]->vertex[2]->z) *
tria[i]->znorm;
    volume=volume-vol/18.0;
    }

    vesicle->volume=volume;
    return TS_SUCCESS;
}
