ts_poly	*init_poly(ts_uint n, ts_vertex *grafted_vtx){
	ts_poly	*poly=(ts_poly *)calloc(1,sizeof(ts_poly));
	poly->vlist = init_vertex_list(n);
	poly->blist = init_bond_list();
	poly->grafted_vtx = grafted_vtx;
	grafted_vtx->grafted_poly = poly;

	ts_uint i;
	for(i=0,i<n-1,i++){
		vtx_add_cneighbour(poly->blist, poly->vlist->vtx[i], poly->vlist->vtx[i+1]);
	}

	return poly;
}


ts_poly_list *init_poly_list(ts_uint n_poly, ts_uint n_mono, ts_vertex_list vlist){
	ts_poly_list *poly_list=(ts_poly_list *)calloc(1,sizeof(ts_poly_list));

	if (n_poly > vlist->n){fatal("Number of polymers larger then numbero f vertices on a vesicle.",310);}
	
	ts_uint i=0;
	ts_uint gvtxi;
	while(i<n_poly){
		gvtxi = rand() % vlist->n
		if (vlist->vtx[gvtxi]->grafted_poly == NULL){
		poly_list->poly = init_poly(n_mono, vlist->vtx[gvtxi]);
		i++;
		}
	}
	
	poly_list->n = n_poly;

	return poly_list;
}

