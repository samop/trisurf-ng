#include<math.h>
#include<stdlib.h>
#include "general.h"
#include "sh.h"

/* Gives you legendre polynomials. Taken from NR, p. 254 */
ts_double plgndr(ts_int l, ts_int m, ts_float x){
	ts_double fact, pll, pmm, pmmp1, somx2;
	ts_int i,ll;

#ifdef TS_DOUBLE_DOUBLE
	if(m<0 || m>l || fabs(x)>1.0)
		fatal("Bad arguments in routine plgndr",1);
#endif
#ifdef TS_DOUBLE_FLOAT
	if(m<0 || m>l || fabsf(x)>1.0)
		fatal("Bad arguments in routine plgndr",1);
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
	if(m<0 || m>l || fabsl(x)>1.0)
		fatal("Bad arguments in routine plgndr",1);
#endif
	pmm=1.0;
	if (m>0) {
#ifdef TS_DOUBLE_DOUBLE
		somx2=sqrt((1.0-x)*(1.0+x));
#endif
#ifdef TS_DOUBLE_FLOAT
		somx2=sqrtf((1.0-x)*(1.0+x));
#endif
#ifdef TS_DOUBLE_LONGDOUBLE
		somx2=sqrtl((1.0-x)*(1.0+x));
#endif
		fact=1.0;
		for (i=1; i<=m;i++){
			pmm *= -fact*somx2;
			fact +=2.0;
		}
	}

	if (l == m) return pmm;
	else {
		pmmp1=x*(2*m+1)*pmm;
		if(l==(m+1)) return(pmmp1);
		else {
			pll=0; /* so it can not be uninitialized */
			for(ll=m+2;ll<=l;ll++){
				pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
				pmm=pmmp1;
				pmmp1=pll;
			}
			return(pll);
		}
	}
}


/*Computes Y(l,m,theta,fi) (Miha's definition that is different from common definition for  factor srqt(1/(2*pi)) */
ts_double shY(ts_int l,ts_int m,ts_double theta,ts_double fi){
	ts_double fac1, fac2, K;
	int i;

	if(l<0 || m>l || m<-l)
		fatal("Error using shY function!",1);

	fac1=1.0;
	for(i=1; i<=l-abs(m);i++){
		fac1 *= i;
	}
	fac2=1.0;
	for(i=1; i<=l+abs(m);i++){
		fac2 *= i;
	}

	if(m==0){
		K=sqrt(1.0/(2.0*M_PI));
	}
	else if (m>0) {
		K=sqrt(1.0/(M_PI))*cos(m*fi);
	} 
	else {
		//K=pow(-1.0,abs(m))*sqrt(1.0/(2.0*M_PI))*cos(m*fi);
		if(abs(m)%2==0)
		K=sqrt(1.0/(M_PI))*cos(m*fi);
		else
		K=-sqrt(1.0/(M_PI))*cos(m*fi);
	}
	
	return K*sqrt((2.0*l+1.0)/2.0*fac1/fac2)*plgndr(l,abs(m),cos(theta));	
}
