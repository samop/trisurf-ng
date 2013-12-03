#ifndef _POLY_H
#define _POLY_H

#include"general.h"

ts_poly	*init_poly(ts_uint n, ts_vertex *grafted_vtx);

ts_poly_list *init_poly_list(ts_uint n_poly, ts_uint n_mono, ts_vertex_list *vlist);

#endif
