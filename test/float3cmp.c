#include<stdio.h>

/** Compares positions of the vertices in two files to given accuracy
  *
  * @author Samo Penic
  * @date to be done
  *
  * Usage: float3cmp vertexfile1 vertexfile2 eps
  *
  * Compares files with 3 floats in a line for similarities. Result of
  * comparison is reported on stdout.
*/
int main(int argv, char *argc[]){

	if(argv!=4){
		fprintf(stderr,"Error. Usage: float3cmp vertexfile1 vertexfile2 eps\n");
		exit(1);
	}
	
	fprintf("Comparing %s and %s for equality. Criterion: abs(float1-float2)< %f\n", argc[1], argc[2],argc[3]);

return 0;
} 
