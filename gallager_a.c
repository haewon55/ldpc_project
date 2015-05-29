#include<stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <time.h>

#define ITER_CHUNK 1000


int n, k;	// size of parity check matrix
int iter;	// number of iterations per one signal
int sig_iter; // number of signals to try
char filename[100];	// file name for parity check matrix
char **parity_matrix;	// parity check matrix
float ber = 0.01;

// ---------------
// NODE STRUCTURES
// ---------------


struct node {
	int deg;
	int *neighbors;
	int *neighbor_ports;
	int port_temp;
	char *recv_msg;
	char *send_msg;
};

int main(int argc, char *argv[]){
	int i, j, l, c;
	int opt;
	FILE *fp;
	while((opt=getopt(argc, argv, "n:k:f:i:e:t:")) != -1){
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
			case 'e':
				ber = atof(optarg);
				break;
			case 't':
				sig_iter = atoi(optarg);
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
				fprintf(stderr, "Usage: %s [-n n] [-k k] [-i #iterations] [-e bit-error rate] [-t #signals] [-f filename]\n", argv[0]);
				exit(EXIT_FAILURE);
		}
	}

	/* Print out informations */
	printf("Testing (%d, %d) code\n", n, k);
	printf("Bit error rate : %.5f\n", ber);
	printf("# iterations in message passing algorithm: %d", iter);
	printf("# of signals to try: %d\n", sig_iter);
	printf("Reading parity matrix from the file: %s\n", filename);

	/* Reading parity matrix from the file */ 
	parity_matrix = malloc(k * sizeof(char*));
	for (i=0; i < k; i++){
		parity_matrix[i] = malloc(n*sizeof(char));
	}
	
	for (i=0; i < k; i++){
		fscanf(fp, "%s", parity_matrix[i]);
	}

	for (i=0; i< k; i++){
		for(j=0; j<n; j++){
			if (parity_matrix[i][j]=='1'){
				parity_matrix[i][j] = 1;
			} else{
				parity_matrix[i][j] = 0;
			}
		}
	}
	fclose(fp);

	/* initialize var/check nodes */
	struct node var_nodes[n];
	struct node check_nodes[k];


	/* Initialize degrees.
	   This step can be skipped if it's regular code */
	for (i = 0; i < n; i++){
		var_nodes[i].deg = 0;
	}
	for (j = 0; j < k; j++){
		check_nodes[j].deg = 0;
	}

	for (i = 0; i < k; i++){
		for(j = 0; j < n ; j++){
			if (parity_matrix[i][j] == 1){
				var_nodes[j].deg++;
				check_nodes[i].deg++;
			}
		}
	}

	for (i = 0; i < n; i++){
		var_nodes[i].neighbors = malloc(var_nodes[i].deg * sizeof(int));
		var_nodes[i].neighbor_ports = malloc(var_nodes[i].deg * sizeof(int));
		var_nodes[i].send_msg = malloc(var_nodes[i].deg * sizeof(char));
		var_nodes[i].recv_msg = malloc(var_nodes[i].deg * sizeof(char));
		var_nodes[i].port_temp = 0;
	}

	for (i = 0; i < k; i++){
		check_nodes[i].neighbors = malloc(check_nodes[i].deg * sizeof(int));
		check_nodes[i].neighbor_ports = malloc(check_nodes[i].deg * sizeof(int));
		check_nodes[i].send_msg = malloc(check_nodes[i].deg * sizeof(char));
		check_nodes[i].recv_msg = malloc(check_nodes[i].deg * sizeof(char));
		check_nodes[i].port_temp = 0;
	}

	for(i = 0; i < k; i++){
		for(j = 0; j < n; j++){
			if(parity_matrix[i][j] == 1){
				var_nodes[j].neighbors[var_nodes[j].port_temp] = i;
				check_nodes[i].neighbors[check_nodes[i].port_temp] = j;
				var_nodes[j].neighbor_ports[var_nodes[j].port_temp] = check_nodes[i].port_temp;
				check_nodes[i].neighbor_ports[check_nodes[i].port_temp] = var_nodes[j].port_temp;
				var_nodes[j].port_temp++;
				check_nodes[i].port_temp++;
			}
		}
	}

/*
	// Test nodes initialization
	for (i = 0; i < k; i++){
		printf("variable node %d deg: %d\n", i, check_nodes[i].deg);
		for (j = 0; j< check_nodes[i].deg; j++){
			printf("%d  ", check_nodes[i].neighbors[j]);
		}
		printf("\n");
	}
*/


	int t;
	int parity;
	/* TODO: Need improvement.
	   Don't know how to deal with this threshold. */
	int b = 9;
	int iter_chunks = 0;

	FILE *write_fp;
	char write_filename[100];
	sprintf(write_filename, "%d_%d_%diter_%dtry_%.4fBER_result", n, k,iter,sig_iter, ber); 
	write_fp = fopen(write_filename, "w");
	
	char signal[n];
	char belief[n];
	int error_per_iter;
	time_t init_time;
	float r;
	srand((unsigned)time(&init_time));

	clock_t start, end;
    double cpu_time_used;
	      
	start = clock();
	
	for (t = 0; t < sig_iter; t++){
	if ((t+1) % ITER_CHUNK == 0){
		printf("%d Done.\n", ++iter_chunks * ITER_CHUNK);
	}	

	/* Signal initialization */
	error_per_iter = 0;
	for (i = 0; i < n; i++){
		r = (double)rand() / (double)RAND_MAX;
		if (r < ber){
			signal[i] = 1;
			error_per_iter++;
		} else {
			signal[i] = 0;
		}	
		belief[i] = signal[i];
		/* Setting first variable messages */
		for(c = 0; c < var_nodes[i].deg; c++){
			var_nodes[i].send_msg[c] = signal[i];
		}
	}
	
	// fprintf(write_fp, "%d\n", error_per_iter);


	/* Start message passing algorithm */
	for (l = 0; l < iter; l++){
		/* check node operations
			- receive message from variable nodes
			- perform parity checking
			- generate send messages
	   */
		error_per_iter = 0;
		for(j = 0; j < k; j++){
			parity = 0;
			for (c = 0; c < check_nodes[j].deg; c++){
				check_nodes[j].recv_msg[c] = var_nodes[check_nodes[j].neighbors[c]].send_msg[check_nodes[j].neighbor_ports[c]];
				parity += check_nodes[j].recv_msg[c]; 
			}
			for (c = 0; c < check_nodes[j].deg; c++){
				check_nodes[j].send_msg[c] = (parity + check_nodes[j].recv_msg[c]) % 2;
			}
		}

		/* variable node operations
		    - receive message from check nodes
			- generate send messages
			- beleif update
		*/
		for (i = 0; i < n; i++){
			parity = 0;
			for(c = 0; c < var_nodes[i].deg; c++){
				var_nodes[i].recv_msg[c] = check_nodes[var_nodes[i].neighbors[c]].send_msg[var_nodes[i].neighbor_ports[c]];
				parity += var_nodes[i].recv_msg[c];
			}
			for(c = 0; c < var_nodes[i].deg; c++){
				if ((parity - var_nodes[i].recv_msg[c]) >=  b){
					var_nodes[i].send_msg[c] = 1;
				} else if((parity - var_nodes[i].recv_msg[c]) < (var_nodes[i].deg - b) ){
					var_nodes[i].send_msg[c] = 0;
				} else{
					var_nodes[i].send_msg[c] = signal[i];
				}
			}

			if (parity >= (var_nodes[i].deg/2)){
				belief[i] = 1;
				error_per_iter++;
			} else {
				belief[i] = 0;
			}
		}
		fprintf(write_fp, "%d\t", error_per_iter);
	}
	fprintf(write_fp, "\n");
	}
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("%d completed\n", sig_iter);
	printf("%.2f s elapsed\n", cpu_time_used);

	fclose(write_fp);
	return 0;
}




