#ifndef _H_SH
#define _H_SH
#include "general.h"

ts_double plgndr(ts_int l, ts_int m, ts_float x);
ts_double shY(ts_int l,ts_int m,ts_double theta,ts_double fi);


ts_bool precomputeShCoeff(ts_spharm *sph);
ts_spharm *sph_init(ts_vertex_list *vlist, ts_uint l);
ts_bool sph_free(ts_spharm *sph);
ts_double getR0(ts_vesicle *vesicle);
ts_bool preparationSh(ts_vesicle *vesicle, ts_double r0);
ts_bool calculateYlmi(ts_vesicle *vesicle);
#endif
