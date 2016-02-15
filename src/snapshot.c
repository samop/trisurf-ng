#include<stdio.h>
#include<general.h>
#include<snapshot.h>
ts_bool xml_trisurf_data(FILE *fh, ts_vesicle *vesicle){
	xml_trisurf_header(fh, vesicle);
	xml_trisurf_tria(fh,vesicle->tlist);
	xml_trisurf_tria_neigh(fh,vesicle->tlist);
	xml_trisurf_vtx_neigh(fh,vesicle->vlist);
	xml_trisurf_vtx_tristar(fh,vesicle->vlist);
	xml_trisurf_footer(fh);
	return TS_SUCCESS;
}

ts_bool xml_trisurf_header(FILE *fh, ts_vesicle *vesicle){
	fprintf(fh, "<trisurf nvtx=\"%u\" npoly=\"%u\" nfono=\"%u\">\n", vesicle->vlist->n, vesicle->poly_list->n, vesicle->poly_list->poly[0]->vlist->n);
	return TS_SUCCESS;
}

ts_bool xml_trisurf_footer(FILE *fh){
	fprintf(fh, "</trisurf>\n");
	return TS_SUCCESS;
}

ts_bool xml_trisurf_tria(FILE *fh, ts_triangle_list *tlist){
	ts_uint i;
	fprintf(fh,"<tria>\n");
	for(i=0; i<tlist->n;i++){
		fprintf(fh,"%u %u %u\n",tlist->tria[i]->vertex[0]->idx, tlist->tria[i]->vertex[1]->idx, tlist->tria[i]->vertex[2]->idx);
	}
	fprintf(fh,"</tria>\n");
	return TS_SUCCESS;
}

ts_bool xml_trisurf_tria_neigh(FILE *fh, ts_triangle_list *tlist){
	ts_uint i;
	fprintf(fh,"<trianeigh>\n");
	for(i=0; i<tlist->n;i++){
		fprintf(fh,"%u %u %u\n",tlist->tria[i]->neigh[0]->idx, tlist->tria[i]->neigh[1]->idx, tlist->tria[i]->neigh[2]->idx);
	}
	fprintf(fh,"</trianeigh>\n");
	return TS_SUCCESS;
}

ts_bool xml_trisurf_vtx_neigh(FILE *fh, ts_vertex_list *vlist){
	ts_uint i,j;
	for(i=0;i<vlist->n;i++){
		fprintf(fh,"<vtxn idx=\"%u\">",vlist->vtx[i]->idx);
		for(j=0;j<vlist->vtx[i]->neigh_no;j++){
			fprintf(fh,"%u ",vlist->vtx[i]->neigh[j]->idx);
		}
		fprintf(fh, "</vtxn>\n");
	}
	return TS_SUCCESS;
}

ts_bool xml_trisurf_vtx_tristar(FILE *fh, ts_vertex_list *vlist){
	ts_uint i,j;
	for(i=0;i<vlist->n;i++){
		fprintf(fh,"<tristar idx=\"%u\">",vlist->vtx[i]->idx);
		for(j=0;j<vlist->vtx[i]->tristar_no;j++){
			fprintf(fh,"%u ",vlist->vtx[i]->tristar[j]->idx);
		}
		fprintf(fh, "</tristar>\n");
	}
	return TS_SUCCESS;
}

