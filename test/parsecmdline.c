#include<stdio.h>
#include<stdlib.h>
#include<string.h>
int main(int argv, char *argc[]){


	char *commands, *backup, *saveptr, *saveopptr, *command, *operator[2], *operand;
	int i,j;
	if(argv!=2){
		fprintf(stderr, "Error. Usage: parsecmdline cmd1=1,cmd2=2,...\n");
		exit(1);
	}
	commands=(char *)malloc(10000*sizeof(char));
    backup=commands;
	strcpy(commands,argc[1]);
	


	for(i=0; ;i++, commands=NULL){
		//breaks comma separated list of commands into specific commands.
		command=strtok_r(commands,",",&saveptr);	
		if(command==NULL) break;
		fprintf(stdout,"Command %d: %s\n",i,command);	
		//extracts name of command and value of command into operator[2] array.
		for(j=0; j<2;j++,command=NULL){
			operator[j]=strtok_r(command,"=",&saveopptr);
			if(operator[j]==NULL) break;
			fprintf(stdout," ---> Operator %d: %s\n",j,operator[j]);		
		}
		//1. check: must have 2 operators.
		if(j!=2) fprintf(stderr,"Error. Command no. %d is not formatted properly.\n",i);
		//2. check: must be named properly.
		//3. check: must be of right format (integer, double, string, ...)


	}

	free(backup);
    return 0;
}
