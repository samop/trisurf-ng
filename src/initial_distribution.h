/** @brief initial bond length */
#define A0 1.2 


/** @brief Creates initial distribution of vertices
  *
  * @param *vlist is the pointer to vertex list
  * @param N is the number of vertices to be initialized
  * @returns TS_SUCCESS on success, TS_FAIL otherwise. If allocation fails, the execution is terminated, reporting error code to the underlying operating system.
*/
ts_vesicle *initial_distribution_dipyramid(ts_uint nshell, ts_uint ncmax1, ts_uint ncmax2, ts_uint ncmax3, ts_double stepsize);


/** Sets the initial position of the vertexes to dipyramid
 *
 *      @param *vlist is a pointer to list of vertices
 *      @returns TS_SUCCESS on success, TS_FAIL otherwise
 */
ts_bool pentagonal_dipyramid_vertex_distribution(ts_vertex_list *vlist);

/** Finds the neighbouring vertices and add them to a list of each vertex
 *
 *      @param *vlist is a pointer to a ts_vertex_list
 *      @returns TS_SUCCESS if successful, TS_FAIL otherwise
 */
ts_bool init_vertex_neighbours(ts_vertex_list *vlist);
//ts_bool init_vertex_neighbours(ts_vertex_list *vlist);

/** interior sites and their neighbours in circ. order + the triangles they are holding together */
ts_bool init_sort_neighbours(ts_vertex_list *vlist);
ts_bool init_vesicle_bonds(ts_vesicle *vesicle);
ts_bool init_triangles(ts_vesicle *vesicle);
ts_bool init_triangle_neighbours(ts_vesicle *vesicle);
ts_bool init_common_vertex_triangle_neighbours(ts_vesicle *vesicle);
ts_bool init_normal_vectors(ts_triangle_list *tlist);

