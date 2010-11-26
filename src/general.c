#include<stdio.h>
#include<stdlib.h>
#include "general.h"
#include<stdarg.h>

ts_uint ts_fprintf(FILE *fd, char *fmt, va_list ap){
if(quiet) return TS_SUCCESS;
fprintf(fd, fmt, ap); /* Call vprintf */
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


