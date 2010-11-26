#include<stdlib.h>
#include<stdio.h>
#include "general.h"
#include "triangle.h"
#include<math.h>

ts_bool init_triangle_list(ts_triangle_list *tlist){
	tlist->n = 0;
	tlist->triangle=NULL;
	return TS_SUCCESS;
}

ts_bool clear_triangle_values(ts_triangle *triang,ts_uint idx){
	triang->idx=idx;
	triang->neigh_no=0;
	triang->vertex[0]=NULL;
	triang->vertex[1]=NULL;
	triang->vertex[2]=NULL;
	triang->neigh=NULL;
    triang->xnorm=0;
    triang->ynorm=0;
    triang->znorm=0;
	return TS_SUCCESS;
}

ts_bool triangle_add(ts_triangle_list *tlist, ts_vertex *vtx1, ts_vertex *vtx2, ts_vertex *vtx3){
        tlist->n++;
        tlist->triangle=realloc(tlist->triangle,tlist->n*sizeof(ts_triangle));
        if(tlist->triangle==NULL) fatal("Cannot reallocate memory for additional ts_triangle.",5);
	clear_triangle_values(&tlist->triangle[tlist->n-1],tlist->n-1);
        //NOW insert vertices!
        tlist->triangle[tlist->n - 1].idx=tlist->n;
        tlist->triangle[tlist->n - 1].vertex[0]=vtx1;
        tlist->triangle[tlist->n - 1].vertex[1]=vtx2;
        tlist->triangle[tlist->n - 1].vertex[2]=vtx3;
        return TS_SUCCESS;
}


ts_bool triangle_add_neighbour(ts_triangle *tria, ts_triangle *ntria){
   // int i;
  //fprintf(stderr,"Sosedi so:\n");
  //  for(i=0;i<tria->neigh_no;i++)
  //      fprintf(stderr,"\t %d\n",tria->neigh[i]);
    

	tria->neigh_no++;
		//We need to reallocate space! The pointer *neight must be zero if not
		//having neighbours yet (if neigh_no was 0 at thime of calling
		tria->neigh=realloc(tria->neigh,tria->neigh_no*sizeof(ts_triangle *));
		if(tria->neigh == NULL){
			fatal("Reallocation of memory failed during insertion of triangle neighbour in triangle_add_neighbour",3);
		}
	tria->neigh[tria->neigh_no-1]=ntria;
   // fprintf(stderr,"dodajamo soseda %d!\n",ntria);

  //fprintf(stderr,"Sosedi so:\n");
  //  for(i=0;i<tria->neigh_no;i++)
  //      fprintf(stderr,"\t %d\n",tria->neigh[i]);
    


	return TS_SUCCESS;
}


ts_bool triangle_remove_neighbour(ts_triangle *tria, ts_triangle *ntria){
    ts_uint i,j=0;
/*    fprintf(stderr,"Sosedi so:\n");
    for(i=0;i<tria->neigh_no;i++)
        fprintf(stderr,"%d, ",tria->neigh[i]);
        fprintf(stderr,"\n");
  */  
    
    for(i=0;i<tria->neigh_no;i++){
        if(tria->neigh[i]!=ntria){
            tria->neigh[j]=tria->neigh[i];
            j++;

        } 
      //  else {
      //      fprintf(stderr,"was here once\n");

      //  }
    }

  //  fprintf(stderr,"Sosedi so:\n");
  //  for(i=0;i<tria->neigh_no;i++)
  //      fprintf(stderr,"\t %d\n",tria->neigh[i]);
    
//    fprintf(stderr,"iscemo soseda %d!\n",ntria);
    if(j==i) fatal("In triangle_remove_neighbour: Specified neighbour does not exist for given triangle",3);
         //   fprintf(stderr,"old nuber of neigh=%i\n",tria->neigh_no);
    tria->neigh_no--;
         //   fprintf(stderr,"new nuber of neigh=%i\n",tria->neigh_no);
    tria->neigh=realloc(tria->neigh,tria->neigh_no*sizeof(ts_triangle *));
	if(tria->neigh == NULL){
			fatal("Reallocation of memory failed during insertion of vertex neighbour in triangle_remove_neighbour",3);
		}
    return TS_SUCCESS;
}

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
	xden=tria->xnorm*tria->xnorm + tria->ynorm*tria->ynorm + tria->znorm*tria->znorm;
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
	return TS_SUCCESS;
}


ts_bool triangle_free(ts_triangle *triang){
    free(triang->neigh);
    return TS_SUCCESS;
}

ts_bool triangle_list_free(ts_triangle_list *tlist){
    int i;
    for(i=0;i<tlist->n;i++){
        triangle_free(&tlist->triangle[i]);
    }
    free(tlist->triangle);  
    return TS_SUCCESS;
}

