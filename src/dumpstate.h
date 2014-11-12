#ifndef _H_DUMPSTATE
#define _H_DUMPSTATE

#include <libxml/parser.h>
#include <libxml/tree.h>



typedef struct {
	long int npoints;
    long int ncells;
    long int idx;
    ts_double *x;
    ts_double *y;
    ts_double *z;

    long int neigh_idx1;
    long int neigh_idx2;
} ts_vtk_data;





ts_vesicle *vtk2vesicle(char *filename, ts_tape *tape);
ts_bool parse_vtk(char *filename, ts_vesicle *vesicle);
ts_bool vtk_index2vesicle(xmlNode *node, ts_vesicle *vesicle);



#endif
