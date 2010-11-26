#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include "general.h"
#include "vertex.h"
#include "bond.h"
#include "initial_timestep.h"

ts_bool initial_distribution(ts_vesicle *vesicle){
    ts_fprintf(stderr,"Starting initial_distribution on vesicle with %u shells!...\n",vesicle->nshell);
	ts_bool retval;
	ts_vertex_list *vlist=&vesicle->vlist;
	ts_bond_list *blist=&vesicle->blist;
	ts_uint nshell=vesicle->nshell;
    ts_uint no_vertices=5*nshell*nshell+2;


    ts_fprintf(stderr,"Calling init_vertex_list...\n");
	retval = init_vertex_list(vlist,no_vertices,2);
	if(retval!=TS_SUCCESS) fatal("There was an error in calling init_vertex_list. Cannot continue.",1);

    retval = vertex_set_global_values(vesicle);
    ts_fprintf(stderr,"Calling init_bond_list...\n");
	retval = init_bond_list(blist);
     ts_fprintf(stderr,"initial_distribution finished!\n");
	return TS_SUCCESS;
} 

