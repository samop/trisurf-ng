#include<stdio.h>
#include<general.h>
#include<snapshot.h>
#include<stdlib.h>
#include<stdarg.h>
#include <zlib.h>
ts_uint ts_sprintf(ts_string *str, char *fmt, ...){
	va_list ap;
	va_start(ap,fmt);
	ts_uint n=vsprintf(&(str->string[str->beg]),fmt,ap);
	va_end(ap);
	str->beg+=n;
	return n;
}


ts_bool xml_trisurf_data(FILE *fh, ts_vesicle *vesicle){

	ts_string *data=(ts_string *)malloc(sizeof(ts_sprintf));
	data->string=(char *)malloc(512000*sizeof(char)); /*TODO: warning, can break if the string is to long */
	data->beg=0;
	
	xml_trisurf_header(fh, vesicle);
	xml_trisurf_tria(data,vesicle->tlist);
	xml_trisurf_tria_neigh(data,vesicle->tlist);
	xml_trisurf_vtx_neigh(data,vesicle->vlist);	
	xml_trisurf_vtx_tristar(data,vesicle->vlist);
#ifdef COMPRESSION
	z_stream defstream;
	defstream.zalloc = Z_NULL;
	defstream.zfree = Z_NULL;
	defstream.opaque = Z_NULL;
	defstream.avail_in = data->beg+1;
	defstream.next_in = (unsigned char *)data->string;	
	char *compr=(char *)malloc((data->beg+1)*sizeof(char *));
	defstream.avail_out = data->beg+1;
	defstream.next_out = (unsigned char *)compr;
	deflateInit(&defstream, Z_BEST_COMPRESSION);
    	deflate(&defstream, Z_FINISH);
    	deflateEnd(&defstream);
	fwrite(compr, sizeof(unsigned char), defstream.total_out, fh);
//	fprintf(fh,"%s",compr);
//printf("Uncompressed size is: %lu\n", defstream.total_out);
	free(compr);
#else
	fprintf(fh,"%s", data->string);
#endif
	free(data->string);
	free(data);
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

ts_bool xml_trisurf_tria(ts_string *data, ts_triangle_list *tlist){
	ts_uint i;
	ts_sprintf(data,"<tria>\n");
	for(i=0; i<tlist->n;i++){
		ts_sprintf(data,"%u %u %u\n",tlist->tria[i]->vertex[0]->idx, tlist->tria[i]->vertex[1]->idx, tlist->tria[i]->vertex[2]->idx);
	}
	ts_sprintf(data,"</tria>\n");
	return TS_SUCCESS;
}

ts_bool xml_trisurf_tria_neigh(ts_string *data, ts_triangle_list *tlist){
	ts_uint i;
	ts_sprintf(data,"<trianeigh>\n");
	for(i=0; i<tlist->n;i++){
		ts_sprintf(data,"%u %u %u\n",tlist->tria[i]->neigh[0]->idx, tlist->tria[i]->neigh[1]->idx, tlist->tria[i]->neigh[2]->idx);
	}
	ts_sprintf(data,"</trianeigh>\n");
	return TS_SUCCESS;
}

ts_bool xml_trisurf_vtx_neigh(ts_string *data, ts_vertex_list *vlist){
	ts_uint i,j;
	for(i=0;i<vlist->n;i++){
		ts_sprintf(data,"<vtxn idx=\"%u\">",vlist->vtx[i]->idx);
		for(j=0;j<vlist->vtx[i]->neigh_no;j++){
			ts_sprintf(data,"%u ",vlist->vtx[i]->neigh[j]->idx);
		}
		ts_sprintf(data, "</vtxn>\n");
	}
	return TS_SUCCESS;
}

ts_bool xml_trisurf_vtx_tristar(ts_string *data, ts_vertex_list *vlist){
	ts_uint i,j;
	for(i=0;i<vlist->n;i++){
		ts_sprintf(data,"<tristar idx=\"%u\">",vlist->vtx[i]->idx);
		for(j=0;j<vlist->vtx[i]->tristar_no;j++){
			ts_sprintf(data,"%u ",vlist->vtx[i]->tristar[j]->idx);
		}
		ts_sprintf(data, "</tristar>\n");
	}
	return TS_SUCCESS;
}

