
#include <string.h>
#include "confuse.h"

int main(void)
{

   long int nshell=17,ncxmax=60, ncymax=60, nczmax=60;
    double xk0=25.0, dmax=1.67,stepsize=0.15;
    cfg_opt_t opts[] = {
        CFG_SIMPLE_INT("nshell", &nshell),
        CFG_SIMPLE_FLOAT("dmax", &dmax),
        CFG_SIMPLE_FLOAT("xk0",&xk0),
        CFG_SIMPLE_FLOAT("stepsize",&stepsize),
        CFG_SIMPLE_INT("nxmax", &ncxmax),
        CFG_SIMPLE_INT("nymax", &ncymax),
        CFG_SIMPLE_INT("nzmax", &nczmax),
        CFG_END()
    };
    cfg_t *cfg;    
    int retval;
    cfg = cfg_init(opts, 0);
    retval=cfg_parse(cfg, "tape_new");

      printf("nshell: %i\n", nshell);
    printf("dmax: %f\n", dmax);
    printf("xk0: %f\n", xk0);
    printf("stepsize: %f\n", stepsize);
    printf("nxmax: %i\n", ncxmax);
    printf("nymax: %i\n", ncymax);
    printf("nzmax: %i\n", nczmax);
    return 0;
}

