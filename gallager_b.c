#include<stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#define MAX_RANDOM     LONG_MAX   // Maximum value of random() 
#define NODES             16384   // Maximum of number of code/check nodes
#define J                    17   // Maximum number of checks per code bit
#define K                    17   // Maximum number of code bits per check

int max_size_M;
int max_size_N;
int size_M[NODES];                // Size of row index sets
int size_N[NODES];                // Size of column index sets
int set_M[NODES][J];              // Column sets for code nodes
int set_N[NODES][K];              // Row sets for check nodes

int n, k;	// size of parity check matrix
int iter;	// number of iterations
char filename[100];	// file name for parity check matrix
char **parity_matrix;	// parity check matrix

// ---------------
// NODE STRUCTURES
// ---------------


struct node {
	int deg;
	int *neighbors;
	char *recv_msg;
	char *send_msg;
};

int main(int argc, char *argv[]){
	int opt;
	FILE *fp;
	while((opt=getopt(argc, argv, "n:k:f:i:")) != -1){
		switch(opt){
			case 'n':
				n = atoi(optarg);
				break;
			case 'k':
				k = atoi(optarg);
				break;
			case 'i':
				iter = atoi(optarg);
				break;
			case 'f':
				sprintf(filename, "./%s", optarg);
				fp = fopen(filename, "r");
				if (fp == NULL){
					fprintf(stderr, "Error opening file:%s\n", strerror(errno));
					exit(EXIT_FAILURE);
				}
				break;
			default: /* '?' */
				fprintf(stderr, "Usage: %s [-n n] [-k k] [-i iterations] [-f filename]\n", argv[0]);
				exit(EXIT_FAILURE);
		}
	}

	// Reading parity matrix 
	parity_matrix = malloc(n * sizeof(char*));
	int i = 0, j = 0;
	for (i=0; i < n; i++){
		parity_matrix[i] = malloc(k*sizeof(char));
	}
 
 
	
	for (i=0; i < n; i++){
		fscanf(fp, "%s", parity_matrix[i]);
	}

	// Test printf for parity matrix reading
	for (i=0; i< n; i++){
		for(j=0; j<k; j++){
			if (parity_matrix[i][j]=='1'){
				parity_matrix[i][j] = 1;
			} else{
				parity_matrix[i][j] = 0;
			}

			printf("%d	", parity_matrix[i][j]);
		}
		printf("\n");
	}

	struct node var_nodes[n];
	struct node check_nodes[k];

	/* Set variable nodes */
	for (i = 0; i < n; i++){
		var_nodes[i].deg = 0;
		/* This step can be skipped if it's regular code */
		for(j = 0; j < k ; j++){
			if (parity_matrix[i][j] == 1){
				var_nodes[i].deg++;
			}
		}

		var_nodes[i].neighbors = malloc(var_nodes[i].deg * sizeof(int));
		var_nodes[i].send_msg = malloc(var_nodes[i].deg * sizeof(char));
		var_nodes[i].recv_msg = malloc(var_nodes[i].deg * sizeof(char));
		
		int c = 0;
		for(j = 0; j <k; j++){
			if (parity_matrix[i][j] ==1){
				var_nodes[i].send_mst[c] = signal[i];
				var_nodes[i].neighbors[c++] = j;
			}
		}		
	}
	
	/* Set check nodes */
	for (i = 0; i < n; i++){
		check_nodes[i].deg = 0;
		/* This step can be skipped if it's regular code */
		for(j = 0; j < k ; j++){
			if (parity_matrix[i][j] == 1){
				check_nodes[i].deg++;
			}
		}

		check_nodes[i].neighbors = malloc(check_nodes[i].deg * sizeof(int));
		check_nodes[i].send_msg = malloc(check_nodes[i].deg * sizeof(char));
		check_nodes[i].recv_msg = malloc(check_nodes[i].deg * sizeof(char));
		
		int c = 0;
		for(j = 0; j <k; j++){
			if (parity_matrix[i][j] ==1){
				check_nodes[i].neighbors[c++] = j;
			}
		}		
	}


	/* Start message passing algorithm */
	int l;
	for (l = 0; l < iter; l++){


	}



	return 0;
}






