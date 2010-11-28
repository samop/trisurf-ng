#ifndef _TRIANGLE_H
#define _TRIANGLE_H

ts_triangle_list *init_triangle_list(void);
ts_triangle *triangle_add(ts_triangle_list *tlist, ts_vertex *vtx1, ts_vertex *vtx2, ts_vertex *vtx3);

/** Adds a neighbouring triangle in a list
 *	@param *tria is a pointer to the triangle, to which additional member want to be added
 *	@param *ntria is a pointer to neighbouring triangle
 *	@returns TS_SUCCESS on success, TS_FAIL otherwise. If memory cannot be alloccated, this is considered as serious error and the execution is immediately terminated with error code returned to the underlying operating system
 */
ts_bool triangle_add_neighbour(ts_triangle *tria, ts_triangle *ntria);
ts_bool triangle_normal_vector(ts_triangle *tria);
ts_bool triangle_data_free(ts_triangle_data *triang);
ts_bool triangle_list_free(ts_triangle_list *tlist);

#endif
