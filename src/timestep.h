#ifndef _TIMESTEP_H
#define _TIMESTEP_H
ts_bool single_timestep(ts_vesicle *vesicle);
ts_bool run_simulation(ts_vesicle *vesicle, ts_uint mcsweeps, ts_uint inititer, ts_uint iterations);
#endif
