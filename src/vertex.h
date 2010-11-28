#ifndef _VERTEX_H
#define _VERTEX_H

/** @brief Creates initial vertex list
 *  
 *  Allocates memory and initializes the vertices.
 *	@param vertex is a structure holding information about 
 *      vertices
 *	@param N is a number of vertices that are used in simulation
 *	@param zero_them is boolean value. 0 skip setting zeros to idx 
 *      and (x,y,z) coordinates for each points, 1 means to zero all 
 *      information on points > 1 requests zeroing of coordinates and 
 *      indexing the vertexes 0..N-1.
 *	@returns ts_bool value 1 on success, 0 otherwise
*/
ts_vertex_list *init_vertex_list(ts_uint N);
ts_vertex_data *init_vertex_data(void);
ts_bool vtx_add_neighbour(ts_vertex *vtx, ts_vertex *nvtx);
ts_bool vtx_add_cneighbour(ts_bond_list *blist,ts_vertex *vtx1,ts_vertex *vtx2);
ts_bool vtx_add_bond(ts_bond_list *blist,ts_vertex *vtx1,ts_vertex *vtx2);
ts_bool vtx_data_free(ts_vertex_data *data);
ts_bool vtx_free(ts_vertex *vtx);
ts_bool vtx_list_free(ts_vertex_list *vlist);
#endif
