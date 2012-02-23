#ifndef _GENERAL_H
#define _GENERAL_H

#include<stdarg.h>
#include<stdio.h>

/* @brief This is a header file, defining general constants and structures.
  * @file header.h
  * @author Samo Penic
  * @date 5.3.2001
  *
  * Header file for general inclusion in all the code, defining data structures
  * and general constans. All datatypes used in the code is also defined here.
  *
  */

/* Defines */
/** @brief Return value of type bz_bool that indiceates successful function finish 
  *
  * Function usualy return some value, which are the result of certain operation. Functions that don't
  * return any parameters can return value, that indicates if the function call finished successfully.
  * In case of successful function run, the functions should return TS_SUCCESS to the caller. This define
  * is set here to get uniformity among all the functions used in program.
  *
  * Example of usage:
  *		ts_boot somefunction(ts_int param1, ....){
  *			...
  *			return TS_SUCCESS;
  *		}
  */
#define TS_SUCCESS 0

/** @brief Return value of type bz_bool that indicates unsuccessful function finish 
  *
  * Function usualy return some value, which are the result of certain operation. Functions that don't
  * return any parameters can return value, that indicates if the function call finished successfully.
  * In case of unsuccessful function run, the functions should return TS_FAIL to the caller. This define
  * is set here to get uniformity among all the functions used in program.
  *
  * Example of usage:
  *
  *		ts_boot somefunction(ts_int param1, ....){
  *			...
  *			return TS_FAIL;
  *		}
  */
#define TS_FAIL 1

/* CONSTANTS */

/* DATA TYPES */
/** @brief Sets the default datatype for ts_double
 *
 * Requred for some functions to work, like "pow" from math.h. If ts_double is defined as
 * float program must run with "powf". Where type dependant function is used it checks this
 * define directive to decide which version to compile in. Available options
 *
 *	TS_DOUBLE_FLOAT
 *	TS_DOUBLE_DOUBLE
 *	TS_DOUBLE_LONGDOUBLE
*/
#define TS_DOUBLE_DOUBLE

/** For the purpose of greater flexibility all data types used in the program
 *  shouldn't use standard C types, but should use types defined here.
 *	ts_int (uses int)
 */
typedef int ts_int;
/** For the purpose of greater flexibility all data types used in the program
 *  shouldn't use standard C types, but should use types defined here.
 *	ts_uint (uses unsigned int)
 */
typedef unsigned int ts_uint;
/** For the purpose of greater flexibility all data types used in the program
 *  shouldn't use standard C types, but should use types defined here.
 *	ts_long (uses long)
 */
typedef long ts_long;
/** For the purpose of greater flexibility all data types used in the program
 *  shouldn't use standard C types, but should use types defined here.
 *	ts_ulong (uses unsigned long)
 */
typedef unsigned long ts_ulong;
/** For the purpose of greater flexibility all data types used in the program
 *  shouldn't use standard C types, but should use types defined here.
 *	ts_float (uses float)
 */
typedef float ts_float;
/** For the purpose of greater flexibility all data types used in the program
 *  shouldn't use standard C types, but should use types defined here.
 *	ts_double (uses double)
 */
typedef double ts_double;
/** For the purpose of greater flexibility all data types used in the program
 *  shouldn't use standard C types, but should use types defined here.
 *	ts_char (uses char)
 */
typedef char ts_char;
/** For the purpose of greater flexibility all data types used in the program
 *  shouldn't use standard C types, but should use types defined here.
 *	ts_uchar (uses unsigned char)
 */
typedef unsigned char ts_uchar;
/** For the purpose of greater flexibility all data types used in the program
 *  shouldn't use standard C types, but should use types defined here.
 *	ts_bool (uses char)
 */
typedef char ts_bool;


/* STRUCTURES */

/** @brief Data structure of all data connected to a vertex
 *
 *  ts_vertex_data holds the data for one single point (bead, vertex). To understand how to use it
 *  here is a detailed description of the fields in the data structure. */
struct ts_vertex_data {
        ts_uint idx; /**< Represents index of the vertex point. Should become obsolete, since it is also present in ts_vertex structure. */        
        ts_double x; /**< The x coordinate of vertex. */
        ts_double y; /**< The y coordinate of vertex. */
        ts_double z; /**< The z coordinate of vertex. */
        ts_uint neigh_no; /**< The number of neighbours. */
        struct ts_vertex **neigh; /**< The pointer that holds neigh_no pointers to this structure. */
        ts_double *bond_length; /**< Obsolete! The bond lenght is moved to ts_bond */
        ts_double *bond_length_dual; /**< Obsolete! Bond length in dual lattice is moved to ts_bond! */
        ts_double curvature;
        ts_double energy;
        ts_double energy_h;
        ts_uint tristar_no;
        struct ts_triangle **tristar; /**< The list of triangles this vertex belongs to. This is an array of pointers to ts_triangle structure of tristar_no length */
        ts_uint bond_no;
        struct ts_bond **bond; /**< Array of pointers of lenght bond_no that stores information on bonds. */
        struct ts_cell *cell; /**< Which cell do we belong to? */
        ts_double xk;
        ts_double c;
        ts_uint id;
};
typedef struct ts_vertex_data ts_vertex_data;

struct ts_vertex {
        ts_uint idx;
        ts_vertex_data *data;
};
typedef struct ts_vertex ts_vertex;

typedef struct {
    ts_uint n;
    ts_vertex **vtx;

} ts_vertex_list;

struct ts_bond {
    ts_uint idx;
	ts_vertex *vtx1;
	ts_vertex *vtx2;
    ts_double bond_length;
    ts_double bond_length_dual;
};
typedef struct ts_bond ts_bond;

struct ts_bond_list {
    ts_uint n;
    ts_bond **bond;
};
typedef struct ts_bond_list ts_bond_list;

struct ts_triangle {
    ts_uint idx;
	ts_vertex *vertex[3];
	ts_uint neigh_no;
	struct ts_triangle **neigh;
	ts_double xnorm;
	ts_double ynorm;
	ts_double znorm;
};
typedef struct ts_triangle ts_triangle;

struct ts_triangle_list{
    ts_uint n;
    ts_triangle **tria;
};
typedef struct ts_triangle_list ts_triangle_list;

typedef struct ts_cell_data {
    ts_vertex **vertex;
    ts_uint nvertex;
} ts_cell_data;

typedef struct ts_cell {
    ts_uint idx;
    ts_cell_data *data;
} ts_cell; 

typedef struct ts_cell_list{
    ts_uint ncmax[3];
    ts_uint cellno;
    ts_cell **cell;
    ts_double dcell;
    ts_double shift;
    ts_double max_occupancy;
} ts_cell_list;


typedef struct {
	ts_vertex_list *vlist;
	ts_bond_list *blist;
	ts_triangle_list *tlist;
    ts_cell_list *clist;
	ts_uint nshell;
    ts_double bending_rigidity;
    ts_double dmax;
    ts_double stepsize;
    ts_double cm[3];
} ts_vesicle;


/* GLOBAL VARIABLES */

int quiet;


/* FUNCTIONS */

/** Non-fatal error function handler:
 *      @param text is a description of an error
 *      @returns doesn't return anything
*/
void err(char *text);

/** Fatal error function handler:
 *      @param text is a description of an error
 *      @param errcode is a (non-zero) error code
 *      @returns terminates the execution of program with errcode set
*/
void fatal(char *text, ts_int errcode);

ts_uint ts_fprintf(FILE *fd, char *fmt, ...);

#define VTX(n) &(vlist->vtx[n])
#define VTX_DATA(n) vlist->vtx[n].data

#endif
