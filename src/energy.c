#include<stdlib.h>
#include "general.h"
#include "energy.h"
#include "vertex.h"
#include<math.h>
#include<stdio.h>
ts_bool mean_curvature_and_energy(ts_vesicle *vesicle){

    ts_uint i, jj, jjp;
    
    ts_vertex_list *vlist=&vesicle->vlist;
    ts_vertex *vtx=vlist->vertex;

    for(i=0;i<vlist->n;i++){
        //should call with zero index!!!
        energy_vertex(&vtx[i]);
    }

    return TS_SUCCESS;
}


inline ts_bool energy_vertex(ts_vertex *vtx){
//    ts_vertex *vtx=&vlist->vertex[n]-1; // Caution! 0 Indexed value!
//    ts_triangle *tristar=vtx->tristar-1;
    ts_uint jj;
    ts_uint jjp,jjm;
    ts_vertex *j,*jp, *jm;
    ts_triangle *jt;
    ts_double s=0,xh=0,yh=0,zh=0,txn=0,tyn=0,tzn=0;
    ts_double x1,x2,x3,ctp,ctm,tot,xlen;
    ts_double h,ht;
    for(jj=1; jj<=vtx->neigh_no;jj++){
        jjp=jj+1;
        if(jjp>vtx->neigh_no) jjp=1;
        jjm=jj-1;
        if(jjm<1) jjm=vtx->neigh_no;
        j=vtx->neigh[jj-1];
        jp=vtx->neigh[jjp-1];
        jm=vtx->neigh[jjm-1];
        jt=vtx->tristar[jj-1];
        x1=vertex_distance_sq(vtx,jp); //shouldn't be zero!
        x2=vertex_distance_sq(j,jp); // shouldn't be zero!
        x3=(j->x-jp->x)*(vtx->x-jp->x)+
           (j->y-jp->y)*(vtx->y-jp->y)+
           (j->z-jp->z)*(vtx->z-jp->z);
        
#ifdef TS_DOUBLE_DOUBLE
        ctp=x3/sqrt(x1*x2-x3*x3);
#endif
#ifdef TS_DOUBLE_FLOAT
        ctp=x3/sqrtf(x1*x2-x3*x3);
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
        ctp=x3/sqrtl(x1*x2-x3*x3);
#endif
        x1=vertex_distance_sq(vtx,jm);
        x2=vertex_distance_sq(j,jm);
        x3=(j->x-jm->x)*(vtx->x-jm->x)+
           (j->y-jm->y)*(vtx->y-jm->y)+
           (j->z-jm->z)*(vtx->z-jm->z);
#ifdef TS_DOUBLE_DOUBLE
        ctm=x3/sqrt(x1*x2-x3*x3);
#endif
#ifdef TS_DOUBLE_FLOAT
        ctm=x3/sqrtf(x1*x2-x3*x3);
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
        ctm=x3/sqrtl(x1*x2-x3*x3);
#endif
        tot=ctp+ctm;
        tot=0.5*tot;
        xlen=vertex_distance_sq(j,vtx);
#ifdef  TS_DOUBLE_DOUBLE 
        vtx->bond_length[jj-1]=sqrt(xlen); 
#endif
#ifdef  TS_DOUBLE_FLOAT
        vtx->bond_length[jj-1]=sqrtf(xlen); 
#endif
#ifdef  TS_DOUBLE_LONGDOUBLE 
        vtx->bond_length[jj-1]=sqrtl(xlen); 
#endif

        vtx->bond_length_dual[jj-1]=tot*vtx->bond_length[jj-1];

        s+=tot*xlen;
        xh+=tot*(j->x - vtx->x);
        yh+=tot*(j->y - vtx->y);
        zh+=tot*(j->z - vtx->z);
        txn+=jt->xnorm;
        tyn+=jt->ynorm;
        tzn+=jt->znorm;
    }
    
    h=xh*xh+yh*yh+zh*zh;
    ht=txn*xh+tyn*yh + tzn*zh;
    s=s/4.0; 
#ifdef TS_DOUBLE_DOUBLE
    if(ht>=0.0) {
        vtx->curvature=sqrt(h);
    } else {
        vtx->curvature=-sqrt(h);
    }
#endif
#ifdef TS_DOUBLE_FLOAT
    if(ht>=0.0) {
        vtx->curvature=sqrtf(h);
    } else {
        vtx->curvature=-sqrtf(h);
    }
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
    if(ht>=0.0) {
        vtx->curvature=sqrtl(h);
    } else {
        vtx->curvature=-sqrtl(h);
    }
#endif
    vtx->energy=0.5*s*(vtx->curvature/s-vtx->c)*(vtx->curvature/s-vtx->c);

    return TS_SUCCESS;
}
