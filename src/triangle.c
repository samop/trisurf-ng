#include<stdlib.h>
#include<stdio.h>
#include "general.h"
#include "triangle.h"
#include<math.h>

ts_triangle_list *init_triangle_list(){
    ts_triangle_list *tlist=(ts_triangle_list *)malloc(sizeof(ts_triangle_list));
	tlist->n = 0;
	tlist->tria=NULL;
	return tlist;
}


ts_triangle *triangle_add(ts_triangle_list *tlist, ts_vertex *vtx1, ts_vertex *vtx2, ts_vertex *vtx3){
        if(vtx1==NULL || vtx2==NULL || vtx3==NULL){
            return NULL;
        }
        tlist->n++;
        tlist->tria=(ts_triangle **)realloc(tlist->tria,tlist->n*sizeof(ts_triangle *));
        if(tlist->tria==NULL) fatal("Cannot reallocate memory for additional ts_triangle.",5);

        tlist->tria[tlist->n-1]=(ts_triangle *)calloc(1,sizeof(ts_triangle));
        if(tlist->tria[tlist->n-1]==NULL) fatal("Cannot reallocate memory for additional ts_triangle.",5);
        tlist->tria[tlist->n-1]->data=(ts_triangle_data *)calloc(1,sizeof(ts_triangle_data));

        //NOW insert vertices!
        tlist->tria[tlist->n - 1]->idx=tlist->n-1;
        tlist->tria[tlist->n - 1]->data->vertex[0]=vtx1;
        tlist->tria[tlist->n - 1]->data->vertex[1]=vtx2;
        tlist->tria[tlist->n - 1]->data->vertex[2]=vtx3;
        return tlist->tria[tlist->n-1];
}


ts_bool triangle_add_neighbour(ts_triangle *tria, ts_triangle *ntria){
    if(tria==NULL || ntria==NULL) return TS_FAIL;
/*TODO: check if the neighbour already exists! Now there is no such check
 * because of the performance issue. */
	tria->data->neigh_no++;
	tria->data->neigh=realloc(tria->data->neigh,tria->data->neigh_no*sizeof(ts_triangle *));
	if(tria->data->neigh == NULL)
			fatal("Reallocation of memory failed during insertion of triangle neighbour in triangle_add_neighbour",3);
	tria->data->neigh[tria->data->neigh_no-1]=ntria;
   
/* we repeat the procedure for the neighbour */  
	ntria->data->neigh_no++;
	ntria->data->neigh=realloc(ntria->data->neigh,ntria->data->neigh_no*sizeof(ts_triangle *));
	if(ntria->data->neigh == NULL)
			fatal("Reallocation of memory failed during insertion of triangle neighbour in triangle_add_neighbour",3);
	ntria->data->neigh[ntria->data->neigh_no-1]=tria;
	return TS_SUCCESS;
}


ts_bool triangle_remove_neighbour(ts_triangle *tria, ts_triangle *ntria){
    ts_uint i,j=0; 
    if(tria==NULL || ntria==NULL) return TS_FAIL;

    for(i=0;i<tria->data->neigh_no;i++){
        if(tria->data->neigh[i]!=ntria){
            tria->data->neigh[j]=tria->data->neigh[i];
            j++;
        } 
    }
    if(j==i) {
        return TS_FAIL; 
        //fatal("In triangle_remove_neighbour: Specified neighbour does not exist for given triangle",3);
    }
    tria->data->neigh_no--;
    tria->data->neigh=(ts_triangle **)realloc(tria->data->neigh,tria->data->neigh_no*sizeof(ts_triangle *));
	if(tria->data->neigh == NULL){
		fatal("Reallocation of memory failed during removal of vertex neighbour in triangle_remove_neighbour",100);
	}
/* we repeat the procedure for neighbour */
    for(i=0;i<ntria->data->neigh_no;i++){
        if(ntria->data->neigh[i]!=tria){
            ntria->data->neigh[j]=ntria->data->neigh[i];
            j++;
        } 
    }
    if(j==i) {
        return TS_FAIL; 
        //fatal("In triangle_remove_neighbour: Specified neighbour does not exist for given triangle",3);
    }
    ntria->data->neigh_no--;
    ntria->data->neigh=(ts_triangle **)realloc(ntria->data->neigh,ntria->data->neigh_no*sizeof(ts_triangle *));
	if(ntria->data->neigh == NULL){
		fatal("Reallocation of memory failed during removal of vertex neighbour in triangle_remove_neighbour",100);
	}
    return TS_SUCCESS;
}

ts_bool triangle_normal_vector(ts_triangle *tria){
	ts_double x21,x31,y21,y31,z21,z31,xden;
	x21=tria->data->vertex[1]->data->x - tria->data->vertex[0]->data->x;
	x31=tria->data->vertex[2]->data->x - tria->data->vertex[0]->data->x;
	y21=tria->data->vertex[1]->data->y - tria->data->vertex[0]->data->y;
	y31=tria->data->vertex[2]->data->y - tria->data->vertex[0]->data->y;
	z21=tria->data->vertex[1]->data->z - tria->data->vertex[0]->data->z;
	z31=tria->data->vertex[2]->data->z - tria->data->vertex[0]->data->z;

	tria->data->xnorm=y21*z31 - z21*y31;
	tria->data->ynorm=z21*x31 - x21*z31;
	tria->data->znorm=x21*y31 - y21*x31;
	xden=tria->data->xnorm*tria->data->xnorm +
         tria->data->ynorm*tria->data->ynorm + 
         tria->data->znorm*tria->data->znorm;
#ifdef TS_DOUBLE_DOUBLE
	xden=sqrt(xden);
#endif
#ifdef TS_DOUBLE_FLOAT
	xden=sqrtf(xden);
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
	xden=sqrtl(xden);
#endif
	tria->data->xnorm=tria->data->xnorm/xden;
	tria->data->ynorm=tria->data->ynorm/xden;
	tria->data->znorm=tria->data->znorm/xden;	
	return TS_SUCCESS;
}


ts_bool triangle_data_free(ts_triangle_data *data){
    if(data->neigh!=NULL) free(data->neigh);
    free(data);
    return TS_SUCCESS;
}

ts_bool triangle_list_free(ts_triangle_list *tlist){
    ts_uint i;
    for(i=0;i<tlist->n;i++){
        triangle_data_free(tlist->tria[i]->data);
        free(tlist->tria[i]);
    }
    free(tlist->tria);
    free(tlist);  
    return TS_SUCCESS;
}

