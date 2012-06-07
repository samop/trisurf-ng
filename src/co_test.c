#include "general.h"
#include "sh.h"
#include "vesicle.h"
#include "initial_distribution.h"

int main(){
ts_vesicle *vesicle=initial_distribution_dipyramid(4,60,60,60,0.15);
//parsetape(vesicle,&i);

//these four must come from parsetype!
vesicle->dmax=1.67*1.67;
vesicle->stepsize=0.15;
vesicle->clist->max_occupancy=8;
vesicle->bending_rigidity=25.0;
fprintf(stderr,"xk=%f\n",vesicle->bending_rigidity);

vesicle->sphHarmonics=sph_init(vesicle->vlist, 10);
int i,j;
for(i=0;i<vesicle->sphHarmonics->l;i++){
    for(j=0;j<2*i+1;j++){
    fprintf(stderr,"co(%d,%d)=%f\n",i,j,vesicle->sphHarmonics->co[i][j]);
    }
}


sph_free(vesicle->sphHarmonics);
//vesicle_free(vesicle);
return 0;
}
