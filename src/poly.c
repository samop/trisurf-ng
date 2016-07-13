/* vim: set ts=4 sts=4 sw=4 noet : */
#include"general.h"
#include"poly.h"
#include<stdlib.h>
#include"vertex.h"
#include"bond.h"
#include<math.h>
#include"energy.h"
#include"cell.h"
ts_bool poly_assign_filament_xi(ts_vesicle *vesicle, ts_tape *tape){
	ts_uint i;

	for(i=0;i<vesicle->filament_list->n;i++){
 	vesicle->filament_list->poly[i]->k = tape->xi;
    	}
	
	return TS_SUCCESS;
}


ts_bool poly_assign_spring_const(ts_vesicle *vesicle){
	ts_uint i;

	for(i=0;i<vesicle->poly_list->n;i++){
 	vesicle->poly_list->poly[i]->k = vesicle->spring_constant;
    	}
	
	return TS_SUCCESS;
}

ts_poly	*init_poly(ts_uint n, ts_vertex *grafted_vtx){
	ts_poly	*poly=(ts_poly *)calloc(1,sizeof(ts_poly));
	poly->vlist = init_vertex_list(n);
	poly->blist = init_bond_list();
	if (grafted_vtx!=NULL){
		poly->grafted_vtx = grafted_vtx;
		grafted_vtx->grafted_poly = poly;
	}

	ts_uint i;
	for(i=0;i<n-1;i++){
		vtx_add_cneighbour(poly->blist, poly->vlist->vtx[i], poly->vlist->vtx[i+1]);
		vtx_add_neighbour(poly->vlist->vtx[i+1], poly->vlist->vtx[i]);
	}

	for(i=0;i<poly->blist->n;i++){
	poly->blist->bond[i]->bond_length=sqrt(vtx_distance_sq(poly->blist->bond[i]->vtx1,poly->blist->bond[i]->vtx2));
	bond_energy(poly->blist->bond[i],poly);
	}

	return poly;
}


ts_poly_list *init_poly_list(ts_uint n_poly, ts_uint n_mono, ts_vertex_list *vlist, ts_vesicle *vesicle){
	ts_poly_list *poly_list=(ts_poly_list *)calloc(1,sizeof(ts_poly_list));
	poly_list->poly	= (ts_poly **)calloc(n_poly,sizeof(ts_poly *));
	ts_uint i=0,j=0,k; //idx;
	ts_uint gvtxi;
	ts_double xnorm,ynorm,znorm,normlength;
	ts_double dphi,dh;

	// Grafting polymers:
	if (vlist!=NULL){
		if (n_poly > vlist->n){fatal("Number of polymers larger than numbero f vertices on a vesicle.",310);}
	
		while(i<n_poly){
			gvtxi = rand() % vlist->n;
			if (vlist->vtx[gvtxi]->grafted_poly == NULL){
			poly_list->poly[i] = init_poly(n_mono, vlist->vtx[gvtxi]);
			i++;
			}
		}
	}
	else
	{
		for(i=0;i<n_poly;i++){
			poly_list->poly[i] = init_poly(n_mono, NULL);
		}
	}

	poly_list->n = n_poly;

	if (vlist!=NULL){
	/* Make straight grafted poylmers normal to membrane (polymer brush). Dist. between poly vertices put to 1*/
		ts_int intpoly=vesicle->tape->internal_poly;
		ts_int cellidx;
		ts_double posX,posY,posZ,prevPosX,prevPosY,prevPosZ, phi,theta;
		ts_bool retval;
		ts_int l;
		for (i=0;i<poly_list->n;i++){
	
			xnorm=0.0;
			ynorm=0.0;
			znorm=0.0;
			for (j=0;j<poly_list->poly[i]->grafted_vtx->tristar_no;j++){
				xnorm-=poly_list->poly[i]->grafted_vtx->tristar[j]->xnorm;
				ynorm-=poly_list->poly[i]->grafted_vtx->tristar[j]->ynorm;
				znorm-=poly_list->poly[i]->grafted_vtx->tristar[j]->znorm;	
			}
			normlength=sqrt(xnorm*xnorm+ynorm*ynorm+znorm*znorm);
			if(intpoly && i%2){
				normlength=-normlength;
			}
			xnorm=xnorm/normlength;
			ynorm=ynorm/normlength;
			znorm=znorm/normlength;

			//prepare starting position for building the polymeres
			prevPosX=poly_list->poly[i]->grafted_vtx->x;
			prevPosY=poly_list->poly[i]->grafted_vtx->y;
			prevPosZ=poly_list->poly[i]->grafted_vtx->z;
			for (j=0;j<poly_list->poly[i]->vlist->n;j++){
				//trying to go towards normal
				posX=prevPosX+(ts_double)xnorm;
				posY=prevPosY+(ts_double)ynorm;
				posZ=prevPosZ+(ts_double)znorm;
				k=0;
				l=0;
				while(1){
					poly_list->poly[i]->vlist->vtx[j]->x = posX;
					poly_list->poly[i]->vlist->vtx[j]->y = posY;
					poly_list->poly[i]->vlist->vtx[j]->z = posZ;
					cellidx=vertex_self_avoidance(vesicle, poly_list->poly[i]->vlist->vtx[j]);
					//distance from neighbors
					//retval=TS_SUCCESS;
					/*for(k=0;k<poly_list->poly[i]->vlist->vtx[j]->neigh_no;k++){
						dist=vtx_distance_sq(poly_list->poly[i]->vlist->vtx[j],poly_list->poly[i]->vlist->vtx[j]->neigh[k]);
						if(dist<1.0 || dist>vesicle->dmax){
							retval=TS_FAIL;
							printf("dist_fail! %e\n", dist);
						}
					}*/
					//if(retval!=TS_FAIL){
						//check occupation number
						retval=cell_occupation_number_and_internal_proximity(vesicle->clist,cellidx,poly_list->poly[i]->vlist->vtx[j]);
					//}
					if(retval==TS_SUCCESS){
						retval=cell_add_vertex(vesicle->clist->cell[cellidx],poly_list->poly[i]->vlist->vtx[j]);
						break;
					}
					else{
				//		printf("%d %d Cannot put the vesicle here. Finding another position\n",i,j);
						theta=drand48()*M_PI-M_PI/2;
						phi=drand48()*2*M_PI;
						posX=prevPosX+sin(phi)*cos(theta);
						posY=prevPosY+sin(phi)*sin(theta);
						posZ=prevPosZ+cos(theta);
						//randomly change the normal.
					}
					k++;
					if(k>1000){
						//lets choose another grafting vertex;
						while(1){
							gvtxi = rand() % vesicle->vlist->n;
							if (vesicle->vlist->vtx[gvtxi]->grafted_poly == NULL){
								ts_fprintf(stdout,"Found new potential grafting vertex %d for poly %d\n",gvtxi,i);
								poly_list->poly[i]->grafted_vtx->grafted_poly=NULL;
								poly_list->poly[i]->grafted_vtx = vesicle->vlist->vtx[gvtxi];
								vesicle->vlist->vtx[gvtxi]->grafted_poly = poly_list->poly[i];
								l++;
								k=0;




			xnorm=0.0;
			ynorm=0.0;
			znorm=0.0;
			int o;
			for (o=0;o<poly_list->poly[i]->grafted_vtx->tristar_no;o++){
				xnorm-=poly_list->poly[i]->grafted_vtx->tristar[o]->xnorm;
				ynorm-=poly_list->poly[i]->grafted_vtx->tristar[o]->ynorm;
				znorm-=poly_list->poly[i]->grafted_vtx->tristar[o]->znorm;	
			}
			normlength=sqrt(xnorm*xnorm+ynorm*ynorm+znorm*znorm);
			if(intpoly && i%2){
				normlength=-normlength;
			}
			xnorm=xnorm/normlength;
			ynorm=ynorm/normlength;
			znorm=znorm/normlength;

			//prepare starting position for building the polymeres
			posX=poly_list->poly[i]->grafted_vtx->x+xnorm;
			posY=poly_list->poly[i]->grafted_vtx->y+ynorm;
			posZ=poly_list->poly[i]->grafted_vtx->z+znorm;
			





								break;
							}
						}
						if(l>1000)
							fatal("Cannot make internal polymeres. No space inside the vesicle?",1001);
					}
				}
				prevPosX=posX;
				prevPosY=posY;
				prevPosZ=posZ;
			}
		}
	}
	else
	{
	/* Make filaments inside the vesicle. Helix with radius... Dist. between poly vertices put to 1*/
	ts_double a,R,H,tantheta,h,r,phi,A0=1.2;

		a = A0*(ts_double)vesicle->nshell;
		R = A0*((ts_double)vesicle->nshell)/(2.0*sin(M_PI/5.0));
		H = sqrt(a*a - R*R);
		tantheta = sqrt(R*R - a*a/4.0)/H;
		
		h = -H + sqrt(vesicle->clist->dmin_interspecies)*1.5;
		r = (H-fabs(h))*tantheta - sqrt(vesicle->clist->dmin_interspecies)*1.5;
		dphi = 2.0*asin(1.0/2.0/r)*1.001;
		dh = dphi/2.0/M_PI*1.001;
		phi=0.0;
		for(i=0;i<poly_list->n;i++){
			for (j=0;j<poly_list->poly[i]->vlist->n;j++){
				h = h + dh;
				r = (H-fabs(h))*tantheta - sqrt(vesicle->clist->dmin_interspecies)*1.5;
				dphi = 2.0*asin(1.0/2.0/r)*1.001;
				dh = dphi/2.0/M_PI*1.001;
				phi+=dphi;
				//ji = j + i*poly_list->poly[i]->vlist->n;
				poly_list->poly[i]->vlist->vtx[j]->x = r*cos(phi);
				poly_list->poly[i]->vlist->vtx[j]->y = r*sin(phi);
				poly_list->poly[i]->vlist->vtx[j]->z = h;// ji*dh - (dh*poly_list->n*poly_list->poly[i]->vlist->n/2.0);
			}
		}
	}

		//index correction for polymeres. Important, since each vtx has to have unique id
/*	idx=vlist->n;
	for(i=0;i<n_poly;i++){
		for(j=0;j<n_mono;j++,idx++){

			poly_list->poly[i]->vlist->vtx[j]->idx=idx;

		}
	}
*/

	return poly_list;
}


ts_bool poly_free(ts_poly *poly){

	if (poly->grafted_vtx!=NULL){
		poly->grafted_vtx->grafted_poly=NULL;
	}
	vtx_list_free(poly->vlist);
	bond_list_free(poly->blist);
	free(poly);

	return TS_SUCCESS;
}

ts_bool poly_list_free(ts_poly_list *poly_list){
	ts_uint i;

	for(i=0;i<poly_list->n;i++){
		poly_free(poly_list->poly[i]);
	}
	free(poly_list->poly);
	free(poly_list);
	
	return TS_SUCCESS;
}


ts_poly *remove_poly_with_index(ts_poly_list *poly_list, ts_uint idx){
	ts_uint i;
	ts_poly *removed_poly=poly_list->poly[idx];

	poly_list->n--; //decrease the total number of polymeres
	for(i=idx;i<poly_list->n;i++){ //move the rest of the polymeres up.
		poly_list->poly[i]=poly_list->poly[i+1];
//		poly_list->poly[idx]->idx=idx;
	}
	
	return removed_poly;
}


ts_bool remove_random_polymeres(ts_poly_list *poly_list, ts_uint number){

	ts_uint i, idx;
	ts_poly *poly;

	ts_poly **new_poly_array;
	if(number>poly_list->n) fatal("The number of polymeres to be removed from the list is greater than the number of total polymeres in the list",999);
	for(i=number;i>0;i--){
		idx=rand() % poly_list->n;
		poly=remove_poly_with_index(poly_list, idx);
		poly_free(poly);
	}
	printf("Addr before %ld\n", (long)poly_list->poly);
	new_poly_array=(ts_poly **)calloc(poly_list->n,sizeof(ts_poly *));
	for(i=0;i<poly_list->n;i++){
		new_poly_array[i]=poly_list->poly[i];
	}
	free(poly_list->poly);
	poly_list->poly=new_poly_array;
	printf("Addr after %ld\n", (long)poly_list->poly);
	return TS_SUCCESS;
}

