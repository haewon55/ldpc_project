#include <stdio.h>
#include <stdlib.h>


int main(){
	int **a;
	a = malloc(4*sizeof(int*));
	int i,j;
	for (i = 0; i < 4; i++){
		a[i] = malloc(3*sizeof(int));
	}
	for (i=0; i< 4; i++){
		for(j=0; j<3; j++){
			a[i][j] = i + j;
		}
	}
	for (i=0; i< 4; i++){
		for(j=0; j<3; j++){
			printf("%d  ", a[i][j]);	
		}
		printf("\n");
	}
return 0;
}
