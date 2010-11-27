#ifndef _GENERAL_H
#define _GENERAL_H

#include<stdarg.h>

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
 *  ts_vertex holds the data for one single point (bead, vertex) in the space. To understand how to use it
 *  here is a detailed description of the fields in the data structure. */
struct ts_vertex_data {
        ts_uint idx; /**< Represents index of the vertex point. Should become obsolete in C. */        
        ts_double x; /**< The x coordinate of vertex. */        
        ts_double y; /**< The y coordinate of vertex. */
        ts_double z; /**< The z coordinate of vertex. */
        ts_uint neigh_no; /**< The number of neighbours. */
        struct ts_vertex **neigh; /**< The pointer that holds neigh_no pointers to this structure. Careful when using pointers to pointers! Also developers do mistakes here.  */
        ts_double *bond_length;
        ts_double *bond_length_dual;
        ts_double curvature;
        ts_double energy;
        ts_double energy_h;
        ts_uint tristar_no;
        struct ts_triangle **tristar;
        struct ts_bond **bond;
        struct ts_cell *cell;
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
    ts_vertex *vtx;

} ts_vertex_list;


/** ts_bond is a structure that describes a bond */
typedef struct {
	ts_vertex *vtx1;
	ts_vertex *vtx2;
    ts_double bond_length;
    ts_double bond_length_dual;
} ts_bond;

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

typedef struct ts_cell {
    ts_uint idx;
    ts_vertex **vertex;
    ts_uint nvertex;
} ts_cell;

typedef struct {
	ts_vertex **vlist;
	ts_bond **blist;
	ts_triangle **tlist;
    ts_cell **clist;
	ts_uint nshell;
    ts_uint nvertex;
    ts_uint nbond;
    ts_uint ntria;
    ts_cell ncell;
    ts_double dcell;
    ts_double shift;
    ts_double max_occupancy;
    ts_uint ncmax[3];
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

//ts_uint ts_fprintf(FILE *fd, char *fmt, va_list ap);

#define VTX(n) &(vlist->vtx[n])
#define VTX_DATA(n) vlist->vtx[n].data

#endif
