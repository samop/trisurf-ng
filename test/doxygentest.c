/** @file doxygentest.c
    @author Samo Penic
    @date 4.3.2010
*/

#include<stdio.h>


/** Defines the number of n/a items */
#define N_A 10

/** @brief Structure that holds some info on three members 
  *
  * This detailed description describes the structure more precisely.
  * One member is also extra documented!
  */
typedef struct{
    int i;  
    char *a; /**< If you think you don't understand this, add a note to member
*/
    long *ptr;
} test;

//typedef t_test test;

/** @brief This is a dummy function to test doxygen documentation 
 *
 *  Some further explanation of the testfunction comes in the next lines
 *  when the brief has been satisfied.
 *      @param param1 is an integer
 *      @param param2 is a float
 *      @return function returns nothing
 *
 *  TODO: Still not fully functional algorithm
 */
void testfunction(int param1, float param2){


    int i; /** This will add an extra line in the documentation block */


}

void main(){
    test i;
}
