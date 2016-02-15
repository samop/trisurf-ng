#ifndef _H_SNAPSHOT
#define _H_SNAPSHOT



typedef struct{
	char *string;
	ts_uint beg;
} ts_string;

#define COMPRESSION

ts_bool xml_trisurf_data(FILE *fh, ts_vesicle *vesicle);
ts_bool xml_trisurf_header(FILE *fh, ts_vesicle *vesicle);
ts_bool xml_trisurf_footer(FILE *fh);
ts_bool xml_trisurf_tria(ts_string *data, ts_triangle_list *tlist);
ts_bool xml_trisurf_tria_neigh(ts_string *data, ts_triangle_list *tlist);
ts_bool xml_trisurf_vtx_neigh(ts_string *data, ts_vertex_list *vlist);
ts_bool xml_trisurf_vtx_tristar(ts_string *data, ts_vertex_list *vlist);
#endif
