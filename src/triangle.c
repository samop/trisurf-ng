#include<stdlib.h>
#include<stdio.h>
#include "general.h"
#include "triangle.h"
#include<math.h>

/** @brief Prepares the list for triangles.
  *
  * Create empty list for holding the information on triangles. Triangles are
  * added later on with triangle_add().
  * Returns pointer to the tlist datastructure it has created. This pointer must
  * be assigned to some variable or it will be lost.
  *
  *
  * Example of usage:
  *     ts_triangle_list *tlist;
  *		tlist=triangle_data_free();
  *
  *	    Initalized data structure for holding the information on triangles.
  *		
  */
ts_triangle_list *init_triangle_list(){
    ts_triangle_list *tlist=(ts_triangle_list *)malloc(sizeof(ts_triangle_list));
	tlist->n = 0;
	tlist->tria=NULL;
	return tlist;
}

/** @brief Add the triangle to the triangle list and create necessary data
 * structures.
  *
  * Add the triangle ts_triangle with ts_triangle_data to the ts_triangle_list.
  * The triangle list is resized, the ts_triangle is allocated and
  * ts_triangle_data is allocated and zeroed. The function receives 4 arguments:
  * ts_triangle_list *tlist as list of triangles and 3 ts_vertex *vtx as
  * vertices that are used to form a triangle. Returns a pointer to newly
  * created triangle. This pointer doesn't need assigning, since it is
  * referenced by triangle list.
  *
  * WARNING: Function can be accelerated a bit by removing the NULL checks.
  * However the time gained by removal doesn't justify the time spent by
  * debugging stupid NULL pointers.
  *
  * Example of usage:
  *		triangle_add(tlist, vlist->vtx[1], vlist->vtx[2], vlist->vtx[3]);
  *
  *	    Creates a triangle with given vertices and puts it into the list.
  *		
  */
ts_triangle *triangle_add(ts_triangle_list *tlist, ts_vertex *vtx1, ts_vertex *vtx2, ts_vertex *vtx3){
        if(vtx1==NULL || vtx2==NULL || vtx3==NULL){
            return NULL;
        }
        tlist->n++;
        tlist->tria=(ts_triangle **)realloc(tlist->tria,tlist->n*sizeof(ts_triangle *));
        if(tlist->tria==NULL) fatal("Cannot reallocate memory for additional ts_triangle.",5);

        tlist->tria[tlist->n-1]=(ts_triangle *)calloc(1,sizeof(ts_triangle));
        if(tlist->tria[tlist->n-1]==NULL) fatal("Cannot reallocate memory for additional ts_triangle.",5);
  //      tlist->tria[tlist->n-1]->data=(ts_triangle_data *)calloc(1,sizeof(ts_triangle_data));

        //NOW insert vertices!
        tlist->tria[tlist->n - 1]->idx=tlist->n-1;
        tlist->tria[tlist->n - 1]->vertex[0]=vtx1;
        tlist->tria[tlist->n - 1]->vertex[1]=vtx2;
        tlist->tria[tlist->n - 1]->vertex[2]=vtx3;
        return tlist->tria[tlist->n-1];
}

/** @brief Add the neigbour to triangles.
  *
  * Add the neigbour to the list of neighbouring triangles. The
  * neighbouring triangles are those, who share two vertices. Function resizes
  * the list and adds the pointer to neighbour. It receives two arguments of
  * ts_triangle type. It then adds second triangle to the list of first
  * triangle, but not the opposite. Upon
  * success it returns TS_SUCCESS, upon detecting NULL pointers 
  * returns TS_FAIL and it FATALY ends when the data structure
  * cannot be resized.
  *
  *
  * WARNING: Function can be accelerated a bit by removing the NULL checks.
  * However the time gained by removal doesn't justify the time spent by
  * debugging stupid NULL pointers.
  *
  * Example of usage:
  *		triangle_remove_neighbour(tlist->tria[3], tlist->tria[4]);
  *
  *	    Triangles 3 and 4 are not neighbours anymore.
  *		
  */

ts_bool triangle_add_neighbour(ts_triangle *tria, ts_triangle *ntria){
    if(tria==NULL || ntria==NULL) return TS_FAIL;
/*TODO: check if the neighbour already exists! Now there is no such check
 * because of the performance issue. */
	tria->neigh_no++;
	tria->neigh=realloc(tria->neigh,tria->neigh_no*sizeof(ts_triangle *));
	if(tria->neigh == NULL)
			fatal("Reallocation of memory failed during insertion of triangle neighbour in triangle_add_neighbour",3);
	tria->neigh[tria->neigh_no-1]=ntria;
  
 
/* we repeat the procedure for the neighbour */  
/*	ntria->data->neigh_no++;
	ntria->data->neigh=realloc(ntria->data->neigh,ntria->data->neigh_no*sizeof(ts_triangle *));
	if(ntria->data->neigh == NULL)
			fatal("Reallocation of memory failed during insertion of triangle neighbour in triangle_add_neighbour",3);
	ntria->data->neigh[ntria->data->neigh_no-1]=tria;
*/
	return TS_SUCCESS;
}

/** @brief Remove the neigbours from triangle.
  *
  * Removes the neigbour from the list of neighbouring triangles. The
  * neighbouring triangles are those, who share two vertices. Function resizes
  * the list and deletes the pointer to neighbour. It receives two arguments of
  * ts_triangle type. It then removes eachother form eachother's list. Upon
  * success it returns TS_SUCCESS, upon failure to find the triangle in the
  * neighbour list returns TS_FAIL and it FATALY ends when the datastructure
  * cannot be resized.
  *
  * WARNING: The function doesn't check whether the pointer is NULL or invalid. It is the
  * job of programmer to make sure the pointer is valid.
  *
  * WARNING: Function is slow. Do not use it often!
  *
  * Example of usage:
  *		triangle_remove_neighbour(tlist->tria[3], tlist->tria[4]);
  *
  *	    Triangles 3 and 4 are not neighbours anymore.
  *		
  */
ts_bool triangle_remove_neighbour(ts_triangle *tria, ts_triangle *ntria){
    ts_uint i,j=0; 
    if(tria==NULL || ntria==NULL) return TS_FAIL;

    for(i=0;i<tria->neigh_no;i++){
        if(tria->neigh[i]!=ntria){
            tria->neigh[j]=tria->neigh[i];
            j++;
        } 
    }
    if(j==i) {
        return TS_FAIL; 
        //fatal("In triangle_remove_neighbour: Specified neighbour does not exist for given triangle",3);
    }
    tria->neigh_no--;
//	fprintf(stderr,"*** tria_number=%d\n",tria->neigh_no);
    tria->neigh=(ts_triangle **)realloc(tria->neigh,tria->neigh_no*sizeof(ts_triangle *));
	if(tria->neigh == NULL){
		fatal("Reallocation of memory failed during removal of vertex neighbour in triangle_remove_neighbour",100);
	}
/* we repeat the procedure for neighbour */
	j=0;
    for(i=0;i<ntria->neigh_no;i++){
        if(ntria->neigh[i]!=tria){
            ntria->neigh[j]=ntria->neigh[i];
            j++;
        } 
    }
    if(j==i) {
        return TS_FAIL; 
        //fatal("In triangle_remove_neighbour: Specified neighbour does not exist for given triangle",3);
    }
    ntria->neigh_no--;
//	fprintf(stderr,"*** ntria_number=%d\n",ntria->neigh_no);
    ntria->neigh=(ts_triangle **)realloc(ntria->neigh,ntria->neigh_no*sizeof(ts_triangle *));
	if(ntria->neigh == NULL){
		fatal("Reallocation of memory failed during removal of vertex neighbour in triangle_remove_neighbour",100);
	}
    return TS_SUCCESS;
}


/** @brief Calculates normal vector of the triangle.

  *
  * Calculate normal vector of the triangle (xnorm, ynorm and znorm) and stores
  * information in underlying ts_triangle_data data_structure.
  *
  * Function receives one argument of type ts_triangle. It should be corectly
  * initialized with underlying data structure of type ts_triangle_data. the
  * result is stored in triangle->data->xnorm, triangle->data->ynorm,
  * triangle->data->znorm. Returns TS_SUCCESS on completion. 
  *
  * NOTE: Function uses math.h library. pow function implementation is selected
  * accordind to the setting in genreal.h
  *
  * Example of usage:
  *		triangle_normal_vector(tlist->tria[3]);
  *
  *	    Computes normals and stores information into tlist->tria[3]->xnorm,
  *	    tlist->tria[3]->ynorm, tlist->tria[3]->znorm.
  *		
  */
ts_bool triangle_normal_vector(ts_triangle *tria){
	ts_double x21,x31,y21,y31,z21,z31,xden;
	x21=tria->vertex[1]->x - tria->vertex[0]->x;
	x31=tria->vertex[2]->x - tria->vertex[0]->x;
	y21=tria->vertex[1]->y - tria->vertex[0]->y;
	y31=tria->vertex[2]->y - tria->vertex[0]->y;
	z21=tria->vertex[1]->z - tria->vertex[0]->z;
	z31=tria->vertex[2]->z - tria->vertex[0]->z;

	tria->xnorm=y21*z31 - z21*y31;
	tria->ynorm=z21*x31 - x21*z31;
	tria->znorm=x21*y31 - y21*x31;
	xden=tria->xnorm*tria->xnorm +
         tria->ynorm*tria->ynorm + 
         tria->znorm*tria->znorm;
#ifdef TS_DOUBLE_DOUBLE
	xden=sqrt(xden);
#endif
#ifdef TS_DOUBLE_FLOAT
	xden=sqrtf(xden);
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
	xden=sqrtl(xden);
#endif
	tria->xnorm=tria->xnorm/xden;
	tria->ynorm=tria->ynorm/xden;
	tria->znorm=tria->znorm/xden;	

/*  Here it is an excellent point to recalculate volume of the triangle and
 *  store it into datastructure. Volume is required at least by constant volume
 *  calculation of vertex move and bondflip and spherical harmonics. */
    tria->volume=(tria->vertex[0]->x+ tria->vertex[1]->x + tria->vertex[2]->x) * tria->xnorm + 
       (tria->vertex[0]->y+ tria->vertex[1]->y + tria->vertex[2]->y) * tria->ynorm + 
    (tria->vertex[0]->z+ tria->vertex[1]->z + tria->vertex[2]->z) * tria->znorm;
    tria->volume=-xden*tria->volume/18.0;
/*  Also, area can be calculated in each triangle */
    tria->area=xden/2;


	return TS_SUCCESS;
}






/** @brief Frees the memory allocated for data structure of triangle list
 * (ts_triangle_list)
  *
  * Function frees the memory of ts_triangle_list previously allocated. It
  * accepts one argument, the address of data structure. It destroys all
  * ts_triangle's in the list with underlying data (by calling
  * triangle_data_free()), and the list itself.
  *
  * Should be used eveytime the deletion of triangle list (created by
  * init_triangle_list() and altered by add_triangle() or remove_triangle()) is desired.
  *
  * WARNING: The function doesn't check whether the pointer is NULL or invalid. It is the
  * job of programmer to make sure the pointer is valid.
  *
  * WARNING: Careful when destroying triangle lists. There could be pointers to
  * that information remaining in structures like vertex_data. This pointers
  * will be rendered invalid by this operation and should not be used anymore.
  *
  * Example of usage:
  *		triangle_list_free(tlist);
  *
  *	    Clears all the information on triangles.
  *		
  */
ts_bool triangle_list_free(ts_triangle_list *tlist){
    ts_uint i;
    for(i=0;i<tlist->n;i++){
    	if(tlist->tria[i]->neigh!=NULL) free(tlist->tria[i]->neigh);
        free(tlist->tria[i]);
    }
    free(tlist->tria);
    free(tlist);  
    return TS_SUCCESS;
}

