#include<stdio.h>
#include<stdlib.h>
#include "general.h"
#include<stdarg.h>

#include <sys/time.h>
#include <unistd.h>
#include <time.h>


ts_uint ts_fprintf(FILE *fd, char *fmt, ...){
if(quiet) return TS_SUCCESS;
	va_list ap;
	va_start(ap,fmt);
	char tmbuf[255];
	struct timeval now;
  	gettimeofday(&now, 0);
	strftime(tmbuf, sizeof tmbuf, "%Y-%m-%d %H:%M:%S", localtime(&now.tv_sec));
fprintf(fd, "[%s] ",tmbuf); 
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


