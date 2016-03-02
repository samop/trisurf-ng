/* vim: set ts=4 sts=4 sw=4 noet : */
#ifndef _ENERGY_H
#define _ENERGY_H
ts_bool mean_curvature_and_energy(ts_vesicle *vesicle);
inline ts_bool energy_vertex(ts_vertex *vtx);
inline ts_bool bond_energy(ts_bond *bond,ts_poly *poly);
#endif
