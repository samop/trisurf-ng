#include<stdio.h>
#include<stdlib.h>
#include "general.h"
#include<stdarg.h>

ts_uint ts_fprintf(FILE *fd, char *fmt, ...){
if(quiet) return TS_SUCCESS;
    va_list ap;
    va_start(ap,fmt);
vfprintf(fd, fmt, ap); /* Call vfprintf */
va_end(ap); /* Cleanup the va_list */
return TS_SUCCESS;
}

void err(char *text){
	ts_fprintf(stderr,"Err: %s\n", text);
}

void fatal(char *text, ts_int errcode){
	ts_fprintf(stderr,"Fatal: %s. TERMINATED!\n", text);
	exit(errcode);
}


