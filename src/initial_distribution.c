#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include "general.h"
#include "vertex.h"
#include "bond.h"
#include "vesicle.h"
#include "vertex.h"
#include "triangle.h"
#include "initial_distribution.h"
#include "energy.h"

ts_vesicle *initial_distribution_dipyramid(ts_uint nshell, ts_uint ncmax1, ts_uint ncmax2, ts_uint ncmax3, ts_double stepsize){
    ts_fprintf(stderr,"Starting initial_distribution on vesicle with %u shells!...\n",nshell);
	ts_bool retval;
    ts_uint no_vertices=5*nshell*nshell+2;
	
    ts_vesicle *vesicle=init_vesicle(no_vertices,ncmax1,ncmax2,ncmax3,stepsize);
    vesicle->nshell=nshell;
    retval = vtx_set_global_values(vesicle);
    retval = pentagonal_dipyramid_vertex_distribution(vesicle->vlist);
    retval = init_vertex_neighbours(vesicle->vlist);
    retval = init_sort_neighbours(vesicle->vlist);
    retval = init_vesicle_bonds(vesicle);
    retval = init_triangles(vesicle);
    retval = init_triangle_neighbours(vesicle);
    retval = init_common_vertex_triangle_neighbours(vesicle);
    retval = mean_curvature_and_energy(vesicle);
 ts_fprintf(stderr,"initial_distribution finished!\n");
	return vesicle;
} 


ts_bool pentagonal_dipyramid_vertex_distribution(ts_vertex_list *vlist){
	/* Some often used relations */
	const ts_double s1= sin(2.0*M_PI/5.0);
	const ts_double s2= sin(4.0*M_PI/5.0);
	const ts_double c1= cos(2.0*M_PI/5.0);
	const ts_double c2= cos(4.0*M_PI/5.0);

	/* Calculates projection lenght of an edge bond to pentagram plane */
	const ts_double xl0=A0/(2.0*sin(M_PI/5.0));
#ifdef TS_DOUBLE_DOUBLE
	const ts_double z0=sqrt(pow(A0,2)-pow(xl0,2));
#endif
#ifdef TS_DOUBLE_FLOAT
	const ts_double z0=sqrtf(powf(A0,2)-powf(xl0,2));
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
	const ts_double z0=sqrtl(powl(A0,2)-powl(xl0,2));
#endif
//	const z0=sqrt(A0*A0 -xl0*xl0); /* I could use pow function but if pow is used make a check on the float type. If float then powf, if long double use powl */

/*placeholder for the pointer to vertex datastructure list... DIRTY: actual pointer points towards invalid address, one position before actual beginning of the list... This is to solve the difference between 1 based indexing in original program in fortran and 0 based indexing in C. All algorithms remain unchanged because of this!*/
	ts_vertex **vtx=vlist->vtx -1 ; 


	ts_uint nshell=(ts_uint)( sqrt((ts_double)(vlist->n-2)/5));
//	printf("nshell=%u\n",nshell);
	ts_uint i,n0; // some for loop prereq
	ts_int j,k;
	ts_double dx,dy; // end loop prereq

	/* topmost vertex */
	vtx[1]->data->x=0.0;
	vtx[1]->data->y=0.0;
	vtx[1]->data->z=z0*(ts_double)nshell;
	
	/* starting from to in circular order on pentagrams */	
	for(i=1;i<=nshell;i++){
		n0=2+5*i*(i-1)/2; //-1 would be for the reason that C index starts from 0 
		vtx[n0]->data->x=0.0;
		vtx[n0]->data->y=(ts_double)i*xl0;
		vtx[n0+i]->data->x=vtx[n0]->data->y*s1;
		vtx[n0+i]->data->y=vtx[n0]->data->y*c1;
		vtx[n0+2*i]->data->x=vtx[n0]->data->y*s2;
		vtx[n0+2*i]->data->y=vtx[n0]->data->y*c2;
		vtx[n0+3*i]->data->x=-vtx[n0+2*i]->data->x;
		vtx[n0+3*i]->data->y=vtx[n0+2*i]->data->y;
		vtx[n0+4*i]->data->x=-vtx[n0+i]->data->x;
		vtx[n0+4*i]->data->y=vtx[n0+i]->data->y;
	}

	/* vertexes on the faces of the dipyramid */
	for(i=1;i<=nshell;i++){
		n0=2+5*i*(i-1)/2; // -1 would be because of C!
		for(j=1;j<=i-1;j++){
			dx=(vtx[n0]->data->x-vtx[n0+4*i]->data->x)/(ts_double)i;
			dy=(vtx[n0]->data->y-vtx[n0+4*i]->data->y)/(ts_double)i;
			vtx[n0+4*i+j]->data->x=(ts_double)j*dx+vtx[n0+4*i]->data->x;
			vtx[n0+4*i+j]->data->y=(ts_double)j*dy+vtx[n0+4*i]->data->y;
		}
		for(k=0;k<=3;k++){ // I would be worried about zero starting of for
			dx=(vtx[n0+(k+1)*i]->data->x - vtx[n0+k*i]->data->x)/(ts_double) i;
			dy=(vtx[n0+(k+1)*i]->data->y - vtx[n0+k*i]->data->y)/(ts_double) i;
			for(j=1; j<=i-1;j++){
				vtx[n0+k*i+j]->data->x= (ts_double)j*dx+vtx[n0+k*i]->data->x;
				vtx[n0+k*i+j]->data->y= (ts_double)j*dy+vtx[n0+k*i]->data->y;
			} 
		} 
	}

	for(i=1;i<=nshell;i++){
		n0= 2+ 5*i*(i-1)/2;
		for(j=0;j<=5*i-1;j++){
		vtx[n0+j]->data->z= z0*(ts_double)(nshell-i);   // I would be worried about zero starting of for
		}
	}

/* for botom part of dipyramide we calculate the positions of vertices */
	for(i=2+5*nshell*(nshell+1)/2;i<=vlist->n;i++){
		vtx[i]->data->x=vtx[vlist->n - i +1]->data->x;
		vtx[i]->data->y=vtx[vlist->n - i +1]->data->y;
		vtx[i]->data->z=-vtx[vlist->n - i +1]->data->z;
	}

	for(i=1;i<=vlist->n;i++){
		for(j=1;j<=vlist->n;j++){
			if(i!=j && vtx_distance_sq(vtx[i],vtx[j])<0.001){
				printf("Vertices %u and %u are the same!\n",i,j);
			}
		}
	}
	return TS_SUCCESS;
}



ts_bool init_vertex_neighbours(ts_vertex_list *vlist){
	ts_vertex **vtx=vlist->vtx -1; // take a look at dipyramid function for comment.
	const ts_double eps=0.001; //TODO: find out if you can use EPS from math.h
	ts_uint i,j;
	ts_double dist2; // Square of distance of neighbours
	/*this is not required if we zero all data in vertex structure at initialization */
	/*if we force zeroing at initialization this for loop can safely be deleted */
	//for(i=1;i<=vlist->n;i++){
	//	vtx[i].neigh_no=0;
	//}
	for(i=1;i<=vlist->n;i++){
		for(j=1;j<=vlist->n;j++){
			dist2=vtx_distance_sq(vtx[i],vtx[j]);
			if( (dist2>eps) && (dist2<(A0*A0+eps))){ 
	//if it is close enough, but not too much close (solves problem of comparing when i==j)
				vtx_add_neighbour(vtx[i],vtx[j]);
			}
		}
	//		printf ("vertex %u ima %u sosedov!\n",i,vtx[i]->data->neigh_no);
	}

	return TS_SUCCESS;
}

// TODO: with new datastructure can be rewritten.
ts_bool init_sort_neighbours(ts_vertex_list *vlist){
	ts_vertex **vtx=vlist->vtx -1; // take a look at dipyramid function for comment.
	ts_uint i,l,j,jj,jjj,k=0;   
    ts_double eps=0.001; // Take a look if EPS from math.h can be used

/*lets initialize memory for temporary vertex_list. Should we write a function instead */
    ts_vertex_list *tvlist=init_vertex_list(vlist->n);
    ts_vertex **tvtx=tvlist->vtx -1;  /* again to compensate for 0-indexing */

	ts_double dist2; // Square of distance of neighbours
    ts_double direct; // Something, dont know what, but could be normal of some kind
	for(i=1;i<=vlist->n;i++){
		k++; // WHY i IS NOT GOOD??
       	vtx_add_neighbour(tvtx[k], tvtx[vtx[i]->data->neigh[0]->idx+1]); //always add 1st
       	jjj=1;
       	jj=1;
       	for(l=2;l<=vtx[i]->data->neigh_no;l++){
           	for(j=2;j<=vtx[i]->data->neigh_no;j++){
               	dist2=vtx_distance_sq(vtx[i]->data->neigh[j-1],vtx[i]->data->neigh[jj-1]);
               	direct=vtx_direct(vtx[i],vtx[i]->data->neigh[j-1],vtx[i]->data->neigh[jj-1]);
               	if( (fabs(dist2-A0*A0)<=eps) && (direct>0.0) && (j!=jjj) ){
           			vtx_add_neighbour(tvtx[k],tvtx[vtx[i]->data->neigh[j-1]->idx+1]);
           			jjj=jj;
           			jj=j;
           			break;
           		}
       		}
       	}	
	}

    for(i=1;i<=vlist->n;i++){
        for(j=1;j<=vtx[i]->data->neigh_no;j++){
            if(vtx[i]->data->neigh_no!=tvtx[i]->data->neigh_no){ //doesn't work with nshell=1!
//                fprintf(stderr,"data1=%u data2=%u\n",vtx[i]->data->neigh_no,tvtx[i]->data->neigh_no);
                fatal("Number of neighbours not the same in init_sort_neighbours.",4);
            }
            //we must correct the pointers in original to point to their
            //neighbours according to indexes. Must be sure not to do it any
            //other way! Also, we need to repair the collection of bonds...
            vtx[i]->data->neigh[j-1]=vtx[tvtx[i]->data->neigh[j-1]->idx+1];
        }
    }

   // Must free memory for temporary vertex array to avoid memory leak! HERE! NOW!
//    free_vertex(tvlist.vertex,tvlist.n);
    vtx_list_free(tvlist);
	return TS_SUCCESS;
}


ts_bool init_vesicle_bonds(ts_vesicle *vesicle){
	ts_vertex_list *vlist=vesicle->vlist;
	ts_bond_list *blist=vesicle->blist;
	ts_vertex **vtx=vesicle->vlist->vtx - 1; // Because of 0 indexing
/* lets make correct clockwise ordering of in nearest neighbour list */
	ts_uint i,j,k;
	for(i=1;i<=vlist->n;i++){
		for(j=i+1;j<=vlist->n;j++){
			for(k=0;k<vtx[i]->data->neigh_no;k++){ // has changed 0 to < instead of 1 and <=
				if(vtx[i]->data->neigh[k]==vtx[j]){  //if addresses matches it is the same
					bond_add(blist,vtx[i],vtx[j]);
					break;
				}
			}
		}
	} 
/* Let's make a check if the number of bonds is correct */
    if((blist->n)!=3*(vlist->n-2)){
        ts_fprintf(stderr,"Number of bonds is %u should be %u!\n", blist->n, 3*(vlist->n-2));
        fatal("Number of bonds is not 3*(no_vertex-2).",4);
    }
	return TS_SUCCESS;
}



ts_bool init_triangles(ts_vesicle *vesicle){
	ts_uint i,j,jj,k;
	ts_vertex **vtx=vesicle->vlist->vtx -1; // difference between 0 indexing and 1 indexing
	ts_triangle_list *tlist=vesicle->tlist;
	ts_double dist, direct;
	ts_double eps=0.001; // can we use EPS from math.h?
	k=0;
	for(i=1;i<=vesicle->vlist->n;i++){
		for(j=1;j<=vtx[i]->data->neigh_no;j++){
			for(jj=1;jj<=vtx[i]->data->neigh_no;jj++){
		//		ts_fprintf(stderr,"%u: (%u,%u) neigh_no=%u ",i,j,jj,vtx[i].neigh_no);
        //      ts_fprintf(stderr,"%e, %e",vtx[i].neigh[j-1]->x,vtx[i].neigh[jj-1]->x);
				dist=vtx_distance_sq(vtx[i]->data->neigh[j-1],vtx[i]->data->neigh[jj-1]);
				direct=vtx_direct(vtx[i],vtx[i]->data->neigh[j-1],vtx[i]->data->neigh[jj-1]);				
				if(fabs(dist-A0*A0)<=eps && direct < 0.0 && vtx[i]->data->neigh[j-1]->idx+1 > i && vtx[i]->data->neigh[jj-1]->idx+1 >i){
					triangle_add(tlist,vtx[i],vtx[i]->data->neigh[j-1],vtx[i]->data->neigh[jj-1]);
				}	
			}	
		}
	}
/* We check if all triangles have 3 vertices and if the number of triangles
 * matches the theoretical value.
 */
	for(i=0;i<tlist->n;i++){
        k=0;
		for(j=0;j<3;j++){
			if(tlist->tria[i]->data->vertex[j]!=NULL)
            k++;
		}
            if(k!=3){
                fatal("Some triangles has less than 3 vertices..",4);
            }   
	} 
    if(tlist->n!=2*(vesicle->vlist->n -2)){
        ts_fprintf(stderr,"The number of triangles is %u but should be %u!\n",tlist->n,2*(vesicle->vlist->n -2));
        fatal("The number of triangles doesn't match 2*(no_vertex -2).",4);
    }
	return TS_SUCCESS;
}



ts_bool init_triangle_neighbours(ts_vesicle *vesicle){
	ts_uint i,j,nobo;
    ts_vertex *i1,*i2,*i3,*j1,*j2,*j3;
//	ts_vertex **vtx=vesicle->vlist->vtx -1; // difference between 0 indexing and 1 indexing
	ts_triangle_list *tlist=vesicle->tlist;
    ts_triangle **tria=tlist->tria -1;
    nobo=0;
    for(i=1;i<=tlist->n;i++){
        i1=tria[i]->data->vertex[0]; 
        i2=tria[i]->data->vertex[1]; 
        i3=tria[i]->data->vertex[2]; 
        for(j=1;j<=tlist->n;j++){
            if(j==i) continue;
            j1=tria[j]->data->vertex[0]; 
            j2=tria[j]->data->vertex[1]; 
            j3=tria[j]->data->vertex[2]; 
            if((i1==j1 && i3==j2) || (i1==j2 && i3==j3) || (i1==j3 && i3==j1)){
                    triangle_add_neighbour(tria[i],tria[j]);
                    nobo++;
            }
        }
    }
    for(i=1;i<=tlist->n;i++){
        i1=tria[i]->data->vertex[0]; 
        i2=tria[i]->data->vertex[1]; 
        i3=tria[i]->data->vertex[2]; 
        for(j=1;j<=tlist->n;j++){
            if(j==i) continue;
            j1=tria[j]->data->vertex[0]; 
            j2=tria[j]->data->vertex[1]; 
            j3=tria[j]->data->vertex[2]; 
            if((i1==j1 && i2==j3) || (i1==j3 && i2==j2) || (i1==j2 && i2==j1)){
                triangle_add_neighbour(tria[i],tria[j]);
                nobo++;
            }
        }
    }
    for(i=1;i<=tlist->n;i++){
        i1=tria[i]->data->vertex[0]; 
        i2=tria[i]->data->vertex[1]; 
        i3=tria[i]->data->vertex[2]; 
        for(j=1;j<=tlist->n;j++){
            if(j==i) continue;
            j1=tria[j]->data->vertex[0]; 
            j2=tria[j]->data->vertex[1]; 
            j3=tria[j]->data->vertex[2]; 
            if((i2==j1 && i3==j3) || (i2==j3 && i3==j2) || (i2==j2 && i3==j1)){
                triangle_add_neighbour(tria[i],tria[j]);
                nobo++;
            }
        }
    }
    if(nobo != vesicle->blist->n*2) {
            ts_fprintf(stderr,"Number of triangles= %u, number of bonds= %u\n",nobo/2, vesicle->blist->n);
            fatal("Number of triangle neighbour pairs differs from double the number of bonds!",4);
    }
    return TS_SUCCESS;
}


ts_bool init_common_vertex_triangle_neighbours(ts_vesicle *vesicle){
	ts_uint i,j,jp,k;
    ts_vertex *k1,*k2,*k3,*k4,*k5;
	ts_vertex **vtx=vesicle->vlist->vtx -1; // difference between 0 indexing and 1 indexing
	ts_triangle_list *tlist=vesicle->tlist;
    ts_triangle **tria=tlist->tria -1;

    for(i=1;i<=vesicle->vlist->n;i++){
        for(j=1;j<=vtx[i]->data->neigh_no;j++){
            k1=vtx[i]->data->neigh[j-1];
            jp=j+1;
            if(j == vtx[i]->data->neigh_no) jp=1;
            k2=vtx[i]->data->neigh[jp-1];
            for(k=1;k<=tlist->n;k++){		// VERY NON-OPTIMAL!!! too many loops (vlist.n * vtx.neigh * tlist.n )!
                k3=tria[k]->data->vertex[0];
                k4=tria[k]->data->vertex[1];
                k5=tria[k]->data->vertex[2];
//                ts_fprintf(stderr,"%u %u: k=(%u %u %u)\n",k1,k2,k3,k4,k5);
                if((vtx[i]==k3 && k1==k4 && k2==k5) ||
                (vtx[i]==k4 && k1==k5 && k2==k3) ||
                (vtx[i]==k5 && k1==k3 && k2==k4)){
          //          ts_fprintf(stderr, "Added to tristar! ");
                    vertex_add_tristar(vtx[i],tria[k]);
                }
            }
        }
/*        ts_fprintf(stderr,"TRISTAR for %u (%u):",i-1,vtx[i].tristar_no);
        for(j=0;j<vtx[i].tristar_no;j++){
            ts_fprintf(stderr," %u,",vtx[i].tristar[j]->idx);
        }
        ts_fprintf(stderr,"\n"); */
    }
    return TS_SUCCESS;
}


ts_bool init_normal_vectors(ts_triangle_list *tlist){
	/* Normals point INSIDE vesicle */
	ts_uint k;
	ts_triangle **tria=tlist->tria -1; //for 0 indexing
	for(k=1;k<=tlist->n;k++){
		triangle_normal_vector(tria[k]);	
	}
	return TS_SUCCESS;
}
