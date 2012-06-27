#include "general.h"
#include "sh.h"
#include "vesicle.h"
#include "initial_distribution.h"
#include <math.h>

int main(){
ts_vesicle *vesicle=initial_distribution_dipyramid(17,60,60,60,0.15);
//parsetape(vesicle,&i);

//these four must come from parsetype!
vesicle->dmax=1.67*1.67;
vesicle->stepsize=0.15;
vesicle->clist->max_occupancy=8;
vesicle->bending_rigidity=25.0;
fprintf(stderr,"xk=%f\n",vesicle->bending_rigidity);

vesicle->sphHarmonics=sph_init(vesicle->vlist, 10);
int i,j;
ts_double area;
for(i=1;i<=vesicle->sphHarmonics->l;i++){
    for(j=1;j<=2*i+1;j++){
    fprintf(stderr,"co(%d,%d)=%e\n",i,j,vesicle->sphHarmonics->co[i][j]);
    }
}
ts_double r0;
vesicle_volume(vesicle);
fprintf(stderr,"Volume=%e\n",vesicle->volume);
r0=getR0(vesicle);
fprintf(stderr,"r0=%e\n",r0);
area=0;
for(i=0;i<vesicle->tlist->n;i++){
	area+=vesicle->tlist->tria[i]->area;
}
fprintf(stderr,"area_dipyramid=%e\n",area);

preparationSh(vesicle,r0);
calculateYlmi(vesicle);
calculateUlm(vesicle);

for(i=0;i<vesicle->sphHarmonics->l;i++){
    for(j=0;j<2*i+1;j++){
    fprintf(stderr,"ulm(%d,%d)=%e\n",i,j+1,vesicle->sphHarmonics->ulm[i][j]);
    }
}



sph_free(vesicle->sphHarmonics);
vesicle_free(vesicle);
return 0;
}
