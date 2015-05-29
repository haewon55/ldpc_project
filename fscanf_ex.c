#include <stdio.h>
#include <stdlib.h>

int main(void *arg){
	FILE *fp;
	fp = fopen("./results/1024_768_10iter_10000try_0.0500BER_CORE3_result", "r");
	int i, ans;
	for (i = 0; i < 100; i++){
		fscanf(fp, "%d\t", &ans);
		printf("%d\n", ans);
	}
		return EXIT_SUCCESS;
}
