#include"general.h"
#include"poly.h"
#include<stdlib.h>
#include"vertex.h"
#include"bond.h"
#include<math.h>


ts_poly	*init_poly(ts_uint n, ts_vertex *grafted_vtx){
	ts_poly	*poly=(ts_poly *)calloc(1,sizeof(ts_poly));
	poly->vlist = init_vertex_list(n);
	poly->blist = init_bond_list();
	poly->grafted_vtx = grafted_vtx;
	grafted_vtx->grafted_poly = poly;

	ts_uint i;
	for(i=0;i<n-1;i++){
		vtx_add_cneighbour(poly->blist, poly->vlist->vtx[i], poly->vlist->vtx[i+1]);
		vtx_add_neighbour(poly->vlist->vtx[i+1], poly->vlist->vtx[i]);
	}

	return poly;
}


ts_poly_list *init_poly_list(ts_uint n_poly, ts_uint n_mono, ts_vertex_list *vlist){
	ts_poly_list *poly_list=(ts_poly_list *)calloc(1,sizeof(ts_poly_list));
	poly_list->poly	= (ts_poly **)calloc(n_poly,sizeof(ts_poly *));
	ts_uint i=0,j=0;
	ts_uint gvtxi;
	ts_double xnorm,ynorm,znorm,normlength;

	if (n_poly > vlist->n){fatal("Number of polymers larger then numbero f vertices on a vesicle.",310);}
	
	while(i<n_poly){
		gvtxi = rand() % vlist->n;
		if (vlist->vtx[gvtxi]->grafted_poly == NULL){
		poly_list->poly[i] = init_poly(n_mono, vlist->vtx[gvtxi]);
		i++;
		}
	}
	
	poly_list->n = n_poly;

/* Make straight poylmers normal to membrane. Dist. between poly vertices put to 1*/
	for (i=0;i<poly_list->n;i++){

		xnorm=0.0;
		ynorm=0.0;
		znorm=0.0;
		for (j=0;j<poly_list->poly[i]->grafted_vtx->tristar_no;j++){
			xnorm+=poly_list->poly[i]->grafted_vtx->tristar[j]->xnorm;
			ynorm+=poly_list->poly[i]->grafted_vtx->tristar[j]->ynorm;
			znorm+=poly_list->poly[i]->grafted_vtx->tristar[j]->znorm;	
		}
		normlength=sqrt(xnorm*xnorm+ynorm*ynorm+znorm*znorm);
		xnorm=xnorm/normlength;
		ynorm=ynorm/normlength;
		znorm=znorm/normlength;

		for (j=0;j<poly_list->poly[i]->vlist->n;j++){
			poly_list->poly[i]->vlist->vtx[j]->x = poly_list->poly[i]->grafted_vtx->x + xnorm*(ts_double)(j+1);
			poly_list->poly[i]->vlist->vtx[j]->y = poly_list->poly[i]->grafted_vtx->y + ynorm*(ts_double)(j+1);
			poly_list->poly[i]->vlist->vtx[j]->z = poly_list->poly[i]->grafted_vtx->z + znorm*(ts_double)(j+1);
		}
	}

	return poly_list;
}


ts_bool poly_free(ts_poly *poly){

	if (poly->grafted_vtx!=NULL){
		poly->grafted_vtx->grafted_poly=NULL;
	}
	vtx_list_free(poly->vlist);
	bond_list_free(poly->blist);
	free(poly);

	return TS_SUCCESS;
}

ts_bool poly_list_free(ts_poly_list *poly_list){
	ts_uint i;

	for(i=0;i<poly_list->n;i++){
		poly_free(poly_list->poly[i]);
	}
	free(poly_list->poly);
	free(poly_list);
	
	return TS_SUCCESS;
}