#ifndef _H_CELL
#define _H_CELL
inline ts_uint vertex_self_avoidance(ts_vesicle *vesicle, ts_vertex *vtx);
ts_bool cell_occupation_number_and_internal_proximity(ts_cell_list *clist,
ts_uint cellidx, ts_vertex *vtx, ts_vertex *tvtx);
#endif
