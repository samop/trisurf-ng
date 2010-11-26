#ifndef _IO_H
#define _IO_H

/** @brief Prints the position of vertices for the whole list
 *  
 *  The function is meant more or less as a debug tool, but can be used in production
 *  environment aswell.
 *  the output is in form of idx: x y z
 *  	@param *vlist is a structure holding information on vertex list.
 *	@returns TS_SUCCESS on successful execution, TS_FAIL otherwise.
*/

ts_bool print_vertex_list(ts_vertex_list *vlist);

/** @brief Prints the neighbours of all the vertices
 *  
 *  The function is meant more or less as a debug tool, but can be used in production
 *  environment aswell.
 *  the output is in form of idx(number of neighbours): (x1,y1,z1) (x2,y2,z2) ...
 *  	@param *vlist is a structure holding information on vertex list.
 *	@returns TS_SUCCESS on successful execution, TS_FAIL otherwise.
*/
ts_bool print_vertex_neighbours(ts_vertex_list *vlist);


/** @brief Function outputs the vetex list file to comply to old fortran format
 *  
 *  	@param *vlist is a list of vertices
 *  	@param *filename is a name of the output file to be created (note that if the file already
 *  exists it will be overwritten.
 */
ts_bool write_vertex_fcompat_file(ts_vertex_list *vlist,ts_char *filename);
ts_bool fprint_vertex_list(FILE *fh,ts_vertex_list *vlist);
ts_bool fprint_tristar(FILE *fh, ts_vesicle *vesicle);
ts_bool fprint_triangle_list(FILE *fh, ts_vesicle *vesicle);
ts_bool fprint_vertex_data(FILE *fh,ts_vertex_list *vlist);
ts_bool fprint_bonds(FILE *fh,ts_vesicle *vesicle);
ts_bool write_dout_fcompat_file(ts_vesicle *vesicle, ts_char *filename);
ts_bool read_tape_fcompat_file(ts_vesicle *vesicle, ts_char *filename);


/** @brief Outputs file in vtk format, compatible with paraview.
 *
 *	@param *vlist is a list of vertices
 *	@param *filename is a name of the output file. If exists, it will be overwritten
 *	@param *text is a description line (max. 255 characters) to be included in the file
 */
ts_bool write_vertex_vtk_file(ts_vesicle *vesicle,ts_char *filename, ts_char *text);

ts_bool parsetape(ts_vesicle *vesicle,ts_uint *iterations);


#endif
