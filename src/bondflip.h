#ifdef _H_BONDFLIP
#define _H_BONDFLIP

ts_bool single_bondflip_timestep(ts_vesicle *vesicle, ts_bond *bond, ts_double *rn);

ts_bool ts_flip_bond(ts_vertex *k,ts_vertex *it,ts_vertex *km, ts_vertex *kp, ts_bond *bond);
#endif