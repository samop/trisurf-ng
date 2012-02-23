#include<stdlib.h>
#include "general.h"
#include "energy.h"
#include "vertex.h"
#include<math.h>
#include<stdio.h>
ts_bool mean_curvature_and_energy(ts_vesicle *vesicle){

    ts_uint i;
    
    ts_vertex_list *vlist=vesicle->vlist;
    ts_vertex **vtx=vlist->vtx;

    for(i=0;i<vlist->n;i++){
        energy_vertex(vtx[i]);
        
    }

    return TS_SUCCESS;
}


inline ts_bool energy_vertex(ts_vertex *vtx){
//    ts_vertex *vtx=&vlist->vertex[n]-1; // Caution! 0 Indexed value!
//    ts_triangle *tristar=vtx->tristar-1;
    ts_vertex_data *data=vtx->data;
    ts_uint jj;
    ts_uint jjp,jjm;
    ts_vertex *j,*jp, *jm;
    ts_triangle *jt;
    ts_double s=0,xh=0,yh=0,zh=0,txn=0,tyn=0,tzn=0;
    ts_double x1,x2,x3,ctp,ctm,tot,xlen;
    ts_double h,ht;
    for(jj=1; jj<=data->neigh_no;jj++){
        jjp=jj+1;
        if(jjp>data->neigh_no) jjp=1;
        jjm=jj-1;
        if(jjm<1) jjm=data->neigh_no;
        j=data->neigh[jj-1];
        jp=data->neigh[jjp-1];
        jm=data->neigh[jjm-1];
//        printf("tristar_no=%u, neigh_no=%u, jj=%u\n",data->tristar_no,data->neigh_no,jj);
        jt=data->tristar[jj-1];
        x1=vtx_distance_sq(vtx,jp); //shouldn't be zero!
        x2=vtx_distance_sq(j,jp); // shouldn't be zero!
        x3=(j->data->x-jp->data->x)*(data->x-jp->data->x)+
           (j->data->y-jp->data->y)*(data->y-jp->data->y)+
           (j->data->z-jp->data->z)*(data->z-jp->data->z);
        
#ifdef TS_DOUBLE_DOUBLE
        ctp=x3/sqrt(x1*x2-x3*x3);
#endif
#ifdef TS_DOUBLE_FLOAT
        ctp=x3/sqrtf(x1*x2-x3*x3);
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
        ctp=x3/sqrtl(x1*x2-x3*x3);
#endif
        x1=vtx_distance_sq(vtx,jm);
        x2=vtx_distance_sq(j,jm);
        x3=(j->data->x-jm->data->x)*(data->x-jm->data->x)+
           (j->data->y-jm->data->y)*(data->y-jm->data->y)+
           (j->data->z-jm->data->z)*(data->z-jm->data->z);
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
        xlen=vtx_distance_sq(j,vtx);
#ifdef  TS_DOUBLE_DOUBLE 
        data->bond[jj-1]->bond_length=sqrt(xlen); 
#endif
#ifdef  TS_DOUBLE_FLOAT
        data->bond[jj-1]->bond_length=sqrtf(xlen); 
#endif
#ifdef  TS_DOUBLE_LONGDOUBLE 
        data->bond[jj-1]->bond_length=sqrtl(xlen); 
#endif

        data->bond[jj-1]->bond_length_dual=tot*data->bond[jj-1]->bond_length;

        s+=tot*xlen;
        xh+=tot*(j->data->x - data->x);
        yh+=tot*(j->data->y - data->y);
        zh+=tot*(j->data->z - data->z);
        txn+=jt->xnorm;
        tyn+=jt->ynorm;
        tzn+=jt->znorm;
    }
    
    h=xh*xh+yh*yh+zh*zh;
    ht=txn*xh+tyn*yh + tzn*zh;
    s=s/4.0; 
#ifdef TS_DOUBLE_DOUBLE
    if(ht>=0.0) {
        data->curvature=sqrt(h);
    } else {
        data->curvature=-sqrt(h);
    }
#endif
#ifdef TS_DOUBLE_FLOAT
    if(ht>=0.0) {
        data->curvature=sqrtf(h);
    } else {
        data->curvature=-sqrtf(h);
    }
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
    if(ht>=0.0) {
        data->curvature=sqrtl(h);
    } else {
        data->curvature=-sqrtl(h);
    }
#endif
// What is vtx->data->c?????????????? Here it is 0!
// c is forced curvature energy for each vertex. Should be set to zero for
// norman circumstances.
    data->energy=0.5*s*(data->curvature/s-data->c)*(data->curvature/s-data->c);

    return TS_SUCCESS;
}
