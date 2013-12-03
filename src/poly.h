#ifndef _POLY_H
#define _POLY_H

#include"general.h"

ts_poly	*init_poly(ts_uint n, ts_vertex *grafted_vtx);

ts_poly_list *init_poly_list(ts_uint n_poly, ts_uint n_mono, ts_vertex_list *vlist);

ts_bool poly_free(ts_poly *poly);

ts_bool poly_list_free(ts_poly_list *poly_list);

#endif
