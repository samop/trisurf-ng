/* vim: set ts=4 sts=4 sw=4 noet : */
#include<stdlib.h>
#include "general.h"
#include "cluster.h"
#include <math.h>




ts_cluster_list *init_cluster_list(){
	ts_cluster_list *cstlist=(ts_cluster_list *)malloc(sizeof(ts_cluster_list));
	cstlist->n=0;
	cstlist->cluster=NULL;
	return cstlist;
}

ts_cluster *new_cluster(ts_cluster_list *cstlist){
	
	cstlist->n++;
	cstlist->cluster=(ts_cluster **)realloc(cstlist->cluster,cstlist->n*sizeof(ts_cluster *));
	if(cstlist->cluster==NULL) fatal("Cannot reallocate memory for additional **ts_cluster.",100);
	cstlist->cluster[cstlist->n-1]=(ts_cluster *)calloc(1,sizeof(ts_cluster));
	if(cstlist->cluster[cstlist->n-1]==NULL) fatal("Cannot allocate memory for additional *ts_cluster.",100);
	return cstlist->cluster[cstlist->n-1];
}

ts_bool cluster_add_vertex(ts_cluster *cluster, ts_vertex *vtx){
	cluster->nvtx++;
	cluster->vtx=(ts_vertex **)realloc(cluster->vtx, cluster->nvtx*sizeof(ts_vertex *));
	cluster->vtx[cluster->nvtx-1]=vtx;
	vtx->cluster=cluster;
	return TS_SUCCESS;
}

ts_bool cluster_free(ts_cluster *cluster){
	if(cluster!=NULL){
		if(cluster->vtx!=NULL)
			free(cluster->vtx);
		free(cluster);
	}
	return TS_SUCCESS;
}

ts_bool cluster_list_free(ts_cluster_list *cstlist){
	ts_uint i;
	if(cstlist!=NULL){
		for(i=0;i<cstlist->n;i++){
			cluster_free(cstlist->cluster[i]);
		}
		free(cstlist);
	}
	return TS_SUCCESS;
}

