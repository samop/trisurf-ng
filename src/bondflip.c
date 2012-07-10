#include<stdlib.h>
#include<math.h>
#include "general.h"
#include "vertex.h"
#include "bond.h"
#include "triangle.h"
#include "vesicle.h"
#include "energy.h"
#include "timestep.h"
#include "cell.h"
#include "bondflip.h"
//#include "io.h"
#include<stdio.h>

ts_bool single_bondflip_timestep(ts_vesicle *vesicle, ts_bond *bond, ts_double *rn){
/*c  Vertex and triangle (lm and lp) indexing for bond flip:
c      +----- k-------+              +----- k ------+
c      |lm1 / | \ lp1 |              |lm1 /   \ lp1 |
c      |  /   |   \   |              |  /       \   |
c      |/     |     \ |     FLIP     |/    lm     \ |
c     km  lm  | lp   kp    --->      km ---------- kp  
c      |\     |     / |              |\    lp     / |  
c      |  \   |   /   |              |  \       /   |
c      |lm2 \ | / lp2 |              |lm2 \   / lp2 |
c      +------it------+              +----- it -----+
c
*/
    ts_vertex *it=bond->vtx1;
    ts_vertex *k=bond->vtx2;
    ts_uint nei,neip,neim;
    ts_uint i; //j;
    ts_double oldenergy, delta_energy;
 //   ts_triangle *lm=NULL,*lp=NULL, *lp1=NULL, *lp2=NULL, *lm1=NULL, *lm2=NULL;

    ts_vertex *kp,*km;

    if(it->neigh_no< 3) return TS_FAIL;
    if(k->neigh_no< 3) return TS_FAIL;
    if(k==NULL || it==NULL){
        fatal("In bondflip, number of neighbours of k or it is less than 3!",999);
    }

    nei=0;
    for(i=0;i<it->neigh_no;i++){ // Finds the nn of it, that is k 
        if(it->neigh[i]==k){
            nei=i;
            break;
        }
    }
    neip=nei+1;  // I don't like it.. Smells like I must have it in correct order
    neim=nei-1;
    if(neip>=it->neigh_no) neip=0;
    if((ts_int)neim<0) neim=it->neigh_no-1; /* casting is essential... If not
there the neim is never <0 !!! */
  //  fprintf(stderr,"The numbers are: %u %u\n",neip, neim);
    km=it->neigh[neim];  // We located km and kp
    kp=it->neigh[neip];

    if(km==NULL || kp==NULL){
        fatal("In bondflip, cannot determine km and kp!",999);
    }

  //  fprintf(stderr,"I WAS HERE! after the 4 vertices are known!\n");

/* test if the membrane is wrapped too much, so that kp is nearest neighbour of
 * km. If it is true, then don't flip! */
    for(i=0;i<km->neigh_no;i++){
        if(km->neigh[i] == kp) return TS_FAIL;
    }
 //   fprintf(stderr,"Membrane didn't wrap too much.. Continue.\n");
/* if bond would be too long, return... */
    if(vtx_distance_sq(km,kp) > vesicle->dmax ) return TS_FAIL;
 //   fprintf(stderr,"Bond will not be too long.. Continue.\n");

/* we make a bond flip. this is different than in original fortran */
// 0. step. Get memory prior the flip
  oldenergy=0;
  oldenergy+=k->xk* k->energy;
  oldenergy+=kp->xk* kp->energy;
  oldenergy+=km->xk* km->energy;
  oldenergy+=it->xk* it->energy;
//  for(i=0;i<k->neigh_no;i++) oldenergy+=k->neigh[i]->xk*k->neigh[i]->energy;
//  for(i=0;i<kp->neigh_no;i++) oldenergy+=kp->neigh[i]->xk*kp->neigh[i]->energy;
//  for(i=0;i<km->neigh_no;i++) oldenergy+=km->neigh[i]->xk*km->neigh[i]->energy;
//  for(i=0;i<it->neigh_no;i++) oldenergy+=it->neigh[i]->xk*it->neigh[i]->energy;
/*
fprintf(stderr,"*** Naslov k=%ld\n",(long)k);
fprintf(stderr,"*** Naslov it=%ld\n",(long)it);
fprintf(stderr,"*** Naslov km=%ld\n",(long)km);
fprintf(stderr,"*** Naslov kp=%ld\n",(long)kp);

for(i=0;i<k->neigh_no;i++)
    fprintf(stderr,"k sosed=%ld\n",(long)k->neigh[i]);
for(i=0;i<it->neigh_no;i++)
    fprintf(stderr,"it sosed=%ld\n",(long)it->neigh[i]);

for(i=0;i<km->neigh_no;i++)
    fprintf(stderr,"km sosed=%ld\n",(long)km->neigh[i]);
for(i=0;i<kp->neigh_no;i++)
    fprintf(stderr,"kp sosed=%ld\n",(long)kp->neigh[i]);


*/
//    fprintf(stderr,"I WAS HERE! Before bondflip!\n");
    ts_flip_bond(k,it,km,kp, bond);
//    fprintf(stderr,"I WAS HERE! Bondflip successful!\n");

/* Calculating the new energy */
  delta_energy=0;
  for(i=0;i<k->neigh_no;i++) energy_vertex(k->neigh[i]);
  for(i=0;i<kp->neigh_no;i++) energy_vertex(kp->neigh[i]);
  for(i=0;i<km->neigh_no;i++) energy_vertex(km->neigh[i]);
  for(i=0;i<it->neigh_no;i++) energy_vertex(it->neigh[i]);
  delta_energy+=k->xk* k->energy;
  delta_energy+=kp->xk* kp->energy;
  delta_energy+=km->xk* km->energy;
  delta_energy+=it->xk* it->energy;
//  for(i=0;i<k->neigh_no;i++) delta_energy+=k->neigh[i]->xk*k->neigh[i]->energy;
//  for(i=0;i<kp->neigh_no;i++) delta_energy+=kp->neigh[i]->xk*kp->neigh[i]->energy;
//  for(i=0;i<km->neigh_no;i++) delta_energy+=km->neigh[i]->xk*km->neigh[i]->energy;
//  for(i=0;i<it->neigh_no;i++) delta_energy+=it->neigh[i]->xk*it->neigh[i]->energy;
  delta_energy-=oldenergy;
 // fprintf(stderr,"I WAS HERE! Got energy!\n");
/* MONTE CARLO */
    if(delta_energy>=0){
#ifdef TS_DOUBLE_DOUBLE
        if(exp(-delta_energy)< drand48() )
#endif
#ifdef TS_DOUBLE_FLOAT
        if(expf(-delta_energy)< (ts_float)drand48())
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
        if(expl(-delta_energy)< (ts_ldouble)drand48())
#endif
        {
            //not accepted, reverting changes
       //     fprintf(stderr,"Failed to move, due to MC\n");

//            ts_flip_bond(km,kp,it,k, bond);
            ts_flip_bond(kp,km,k,it, bond);
                

/*
fprintf(stderr,"*** Naslov k=%d\n",k);
fprintf(stderr,"*** Naslov it=%d\n",it);
fprintf(stderr,"*** Naslov km=%d\n",km);
fprintf(stderr,"*** Naslov kp=%d\n",kp);
for(i=0;i<k->neigh_no;i++)
    fprintf(stderr,"k sosed=%d\n",k->neigh[i]);
for(i=0;i<it->neigh_no;i++)
    fprintf(stderr,"it sosed=%d\n",it->neigh[i]);


for(i=0;i<km->neigh_no;i++)
    fprintf(stderr,"km sosed=%d\n",km->neigh[i]);
for(i=0;i<kp->neigh_no;i++)
    fprintf(stderr,"kp sosed=%d\n",kp->neigh[i]);
*/



        //    fprintf(stderr,"Reverted condition!\n");
            return TS_FAIL;
        }
    }
        //    fprintf(stderr,"Success\n");


/* IF BONDFLIP ACCEPTED, THEN RETURN SUCCESS! */
    return TS_SUCCESS;
}


ts_bool ts_flip_bond(ts_vertex *k,ts_vertex *it,ts_vertex *km, ts_vertex *kp,
ts_bond *bond){

    ts_triangle *lm=NULL,*lp=NULL, *lp1=NULL, *lm2=NULL;
    ts_uint i,j; //lmidx, lpidx;
if(k==NULL || it==NULL || km==NULL || kp==NULL){
    fatal("ts_flip_bond: You called me with invalid pointers to vertices",999);
}
// 1. step. We find lm and lp from k->tristar !
    for(i=0;i<it->tristar_no;i++){
        for(j=0;j<k->tristar_no;j++){
            if((it->tristar[i] == k->tristar[j])){ //ce gre za skupen trikotnik
                if((it->tristar[i]->vertex[0] == km || it->tristar[i]->vertex[1]
== km || it->tristar[i]->vertex[2]== km )){
                lm=it->tristar[i];
         //       lmidx=i;
                }
                else
                {
                lp=it->tristar[i];
         //       lpidx=i;
                }

            }
        }
    }
if(lm==NULL || lp==NULL) fatal("ts_flip_bond: Cannot find triangles lm and lp!",999);

//we look for important triangles lp1 and lm2.

 for(i=0;i<k->tristar_no;i++){
        for(j=0;j<kp->tristar_no;j++){
                if((k->tristar[i] == kp->tristar[j]) && k->tristar[i]!=lp){ //ce gre za skupen trikotnik
                    lp1=k->tristar[i];
            }
        }
}

 for(i=0;i<it->tristar_no;i++){
        for(j=0;j<km->tristar_no;j++){
            if((it->tristar[i] == km->tristar[j]) && it->tristar[i]!=lm){ //ce gre za skupen trikotnik
                    lm2=it->tristar[i];
            } 
        }
    }
/*
// DEBUG TESTING!
fprintf(stderr,"*** Naslov k=%d\n",k);
fprintf(stderr,"*** Naslov it=%d\n",it);
fprintf(stderr,"*** Naslov km=%d\n",km);
fprintf(stderr,"*** Naslov kp=%d\n",kp);

for(i=0;i<k->neigh_no;i++)
    fprintf(stderr,"k sosed=%d\n",k->neigh[i]);
for(i=0;i<it->neigh_no;i++)
    fprintf(stderr,"it sosed=%d\n",it->neigh[i]);


// END DEBUG TESTING!
*/
if(lm2==NULL || lp1==NULL) fatal("ts_flip_bond: Cannot find triangles lm2 and lp1!",999);

/*
//DEBUG TESTING
fprintf(stderr,"1. step: lm, lm2, lp1 and lp found!\n");
fprintf(stderr,"--- Naslov lm=%ld",(long)lm);


fprintf(stderr,"   vtxs(%ld, %ld, %ld)\n",(long)lm->vertex[0],(long)lm->vertex[1], (long)lm->vertex[2]);
fprintf(stderr,"--- Naslov lp=%ld",(long)lp);
fprintf(stderr,"   vtxs(%ld, %ld, %ld)\n",(long)lp->vertex[0],(long)lp->vertex[1], (long)lp->vertex[2]);
fprintf(stderr,"--- Naslov lm2=%ld",(long)lm2);
fprintf(stderr,"   vtxs(%ld, %ld, %ld)\n",(long)lm2->vertex[0],(long)lm2->vertex[1], (long)lm2->vertex[2]);
fprintf(stderr,"--- Naslov lp1=%ld",(long)lp1);
fprintf(stderr,"   vtxs(%ld, %ld, %ld)\n",(long)lp1->vertex[0],(long)lp1->vertex[1], (long)lp1->vertex[2]);

for(i=0;i<lm->neigh_no;i++)
    fprintf(stderr,"lm sosed=%ld\n",(long)lm->neigh[i]);
for(i=0;i<lp->neigh_no;i++)
    fprintf(stderr,"lp sosed=%ld\n",(long)lp->neigh[i]);
// END DEBUG TESTING
*/
/*
// DEBUG TESTING!

for(i=0;i<3;i++){

    if(lp1->neigh[i]==lp) fprintf(stderr,"Nasel sem par lp1->lp\n");
    if(lp->neigh[i]==lp1) fprintf(stderr,"Nasel sem par lp->lp1\n");
    if(lm2->neigh[i]==lm) fprintf(stderr,"Nasel sem par lm2->lm\n");
    if(lm->neigh[i]==lm2) fprintf(stderr,"Nasel sem par lm->lm2\n");
}
// END DEBUG TESTING!
*/


// 2. step. We change the triangle vertices... (actual bond flip)
    for(i=0;i<3;i++) if(lm->vertex[i]== it) lm->vertex[i]= kp;
    for(i=0;i<3;i++) if(lp->vertex[i]== k) lp->vertex[i]= km;
//fprintf(stderr,"2. step: actual bondflip made\n");
// 2a. step. If any changes in triangle calculations must be done, do it here!
//   * normals are recalculated here
    triangle_normal_vector(lp);
    triangle_normal_vector(lm);
//fprintf(stderr,"2a. step: triangle normals recalculated\n");
// 3. step. Correct neighbours in vertex_list


            vtx_remove_neighbour(k,it);
//            vtx_remove_neighbour(it,k);
//fprintf(stderr,"3. step (PROGRESS): removed k and it neighbours\n");
    
            //Tukaj pa nastopi tezava... Kam dodati soseda?
            vtx_insert_neighbour(km,kp,k);
            vtx_insert_neighbour(kp,km,it);
//            vertex_add_neighbour(km,kp); //pazi na vrstni red.
//            vertex_add_neighbour(kp,km);
//fprintf(stderr,"3. step: vertex neighbours corrected\n");

// 3a. step. If any changes to ts_vertex, do it here!
//   bond_length calculatons not required for it is done in energy.c

// 4. step. Correct bond_list (don't know why I still have it!)
            bond->vtx1=km;
            bond->vtx2=kp;
//fprintf(stderr,"4. step: bondlist corrected\n");


// 5. step. Correct neighbouring triangles 
   
    triangle_remove_neighbour(lp,lp1);
  //  fprintf(stderr,".\n");
    triangle_remove_neighbour(lp1,lp);
  //  fprintf(stderr,".\n");
    triangle_remove_neighbour(lm,lm2);
  //  fprintf(stderr,".\n");
    triangle_remove_neighbour(lm2,lm);
   
    triangle_add_neighbour(lm,lp1);    
    triangle_add_neighbour(lp1,lm);
    triangle_add_neighbour(lp,lm2);  //Vrstni red?!
    triangle_add_neighbour(lm2,lp);

//fprintf(stderr,"5. step: triangle neigbours corrected\n");


// 6. step. Correct tristar for vertices km, kp, k and it
            vertex_add_tristar(km,lp);  // Preveri vrstni red!
            vertex_add_tristar(kp,lm);
            vtx_remove_tristar(it,lm);
            vtx_remove_tristar(k,lp);
//fprintf(stderr,"6. step: tristar corrected\n");

/*
//DEBUG TESTING
fprintf(stderr,"--- Naslov lm=%d",lm);


fprintf(stderr,"   vtxs(%d, %d, %d)\n",lm->vertex[0],lm->vertex[1], lm->vertex[2]);
fprintf(stderr,"--- Naslov lp=%d",lp);
fprintf(stderr,"   vtxs(%d, %d, %d)\n",lp->vertex[0],lp->vertex[1], lp->vertex[2]);
fprintf(stderr,"--- Naslov lm2=%d",lm2);
fprintf(stderr,"   vtxs(%d, %d, %d)\n",lm2->vertex[0],lm2->vertex[1], lm2->vertex[2]);
fprintf(stderr,"--- Naslov lp1=%d",lp1);
fprintf(stderr,"   vtxs(%d, %d, %d)\n",lp1->vertex[0],lp1->vertex[1], lp1->vertex[2]);

for(i=0;i<lm->neigh_no;i++)
    fprintf(stderr,"lm sosed=%d\n",lm->neigh[i]);
for(i=0;i<lp->neigh_no;i++)
    fprintf(stderr,"lp sosed=%d\n",lp->neigh[i]);
// END DEBUG TESTING
*/
  energy_vertex(k);
  energy_vertex(kp);
  energy_vertex(km);
  energy_vertex(it);


// END modifications to data structure!


    return TS_SUCCESS;
}
