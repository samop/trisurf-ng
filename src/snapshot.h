#ifndef _H_SNAPSHOT
#define _H_SNAPSHOT

ts_bool xml_trisurf_data(FILE *fh, ts_vesicle *vesicle);
ts_bool xml_trisurf_header(FILE *fh, ts_vesicle *vesicle);
ts_bool xml_trisurf_footer(FILE *fh);
ts_bool xml_trisurf_tria(FILE *fh, ts_triangle_list *tlist);
ts_bool xml_trisurf_tria_neigh(FILE *fh, ts_triangle_list *tlist);
ts_bool xml_trisurf_vtx_neigh(FILE *fh, ts_vertex_list *vlist);
ts_bool xml_trisurf_vtx_tristar(FILE *fh, ts_vertex_list *vlist);
#endif
