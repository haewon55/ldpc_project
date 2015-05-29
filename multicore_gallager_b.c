#define _GNU_SOURCE
#include <sched.h>
#include<stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <pthread.h>
#include <pthread.h>

#define ITER_CHUNK 1000
#define NUM_CORES 8

int n, k;	// size of parity check matrix
int iter;	// number of iterations per one signal
int sig_iter; // number of signals to try
char filename[100];	// file name for parity check matrix
char **parity_matrix;	// parity check matrix
float ber = 0.01;
int num_cores = NUM_CORES;

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


struct thread_data{
	int sig_iter_per_core;
	int cpu;
};



int core_affinitize(int core){
	cpu_set_t *cmask;
	int n;
	int ret;
	/*For NUMA setting*/
	struct bitmask *bmask;
	FILE *fp;
	int phy_id;
	n = sysconf(_SC_NPROCESSORS_ONLN);
	if (core < 0 || core >= (int) n) {
		errno = -EINVAL;
		return -1;
	}
	if((cmask = CPU_ALLOC(n)) == NULL)
		perror("CPU ALLOC");

	CPU_ZERO_S(n, cmask);
	CPU_SET_S(core, n, cmask);
	ret = sched_setaffinity(0, n, cmask);
	CPU_FREE(cmask);
	printf("Core %d affinitzied\n", core);
	return ret;
}


void *thr_func(void *arg){

	struct thread_data *thr_data = (struct thread_data *)arg;
	int sig_iter = thr_data->sig_iter_per_core;
	int cpu = thr_data->cpu;

	core_affinitize(cpu);

	/* initialize var/check nodes */
	int i,j,l,t,c=0;
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
	
	int parity;
	/* TODO: Need improvement.
	   Don't know how to deal with this threshold. */
	int b = 10;
	int iter_chunks = 0;
	FILE *write_fp;
	char write_filename[200];
	sprintf(write_filename, "./results/%s_%diter_%dsig_%.4fBER_CORE%d", filename,iter,sig_iter, ber, cpu); 
	write_fp = fopen(write_filename, "w");
	
	char signal[n];
	//char belief[n];
	int error_per_iter;
	float r;

	for (t = 0; t < sig_iter; t++){
	if ((t+1) % ITER_CHUNK == 0){
		printf("CORE %d: %d Done.\n", cpu, ++iter_chunks * ITER_CHUNK);
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
		//belief[i] = signal[i];
		/* Setting first variable messages */
		for(c = 0; c < var_nodes[i].deg; c++){
			var_nodes[i].send_msg[c] = signal[i];
		}
	}
	
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
				} else if((parity - var_nodes[i].recv_msg[c]) <= (var_nodes[i].deg - b) ){
					var_nodes[i].send_msg[c] = 0;
				} else{
					var_nodes[i].send_msg[c] = signal[i];
				}
			}

			if (parity >= (var_nodes[i].deg/2)){
				//belief[i] = 1;
				error_per_iter++;
			} else {
				//belief[i] = 0;
			}
		}
		fprintf(write_fp, "%d\t", error_per_iter);
	}
	fprintf(write_fp, "\n");
	}
	printf("%d completed\n", sig_iter);

	fclose(write_fp);
	return 0;

}



int main(int argc, char *argv[]){
	int i, j;
	int opt;
	FILE *fp;
	
	while((opt=getopt(argc, argv, "n:k:f:i:e:t:c:")) != -1){
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
			case 'c':
				num_cores = atoi(optarg);
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

	
	time_t init_time;
	srand((unsigned)time(&init_time));
		
	/* Multi thread version */
	pthread_t thr[num_cores];
	int ci, tc_result;
	struct thread_data thread_datas[num_cores];

	for(ci = 0; ci < num_cores; ci++){
		thread_datas[ci].sig_iter_per_core = sig_iter/num_cores;
		if (ci == (num_cores-1)){
			thread_datas[ci].sig_iter_per_core += sig_iter%num_cores;
		}
		thread_datas[ci].cpu = ci;
		if((tc_result = pthread_create(&thr[ci], NULL, thr_func, &thread_datas[ci]))){
			fprintf(stderr, "pthred create error: %d\n", tc_result);
			return EXIT_FAILURE;
		}
	}

	for(ci = 0; ci < num_cores; ci++){
		pthread_join(thr[ci], NULL);
	}

	printf("Summing up the result..\n");
	float ber_at_iter[iter];
	char sum_filename[200];
	for (i = 0; i < iter; i++){
			ber_at_iter[i] = 0.0;
	}
	int sig_iter_per_core, val;
	for(ci = 0; ci < num_cores; ci++){
		sig_iter_per_core = thread_datas[ci].sig_iter_per_core;
		sprintf(sum_filename, "./results/%s_%diter_%dsig_%.4fBER_CORE%d", filename,iter,sig_iter_per_core, ber, ci); 
		fp = fopen(sum_filename, "r");
		for(j = 0; j < sig_iter_per_core; j++){
			for(i = 0; i < iter; i++){
				fscanf(fp, "%d\t", &val);
				ber_at_iter[i] += (double)val/(double)n;
			}
		}
	}
	for(i = 0; i < iter; i++){
		ber_at_iter[i] = (double)ber_at_iter[i]/(double)sig_iter;
	}

	printf("BER at iteration");
	for (i = 0; i < iter; i++){
		printf("  %d\t\t", i+1);
	}
	printf("\n                ");
	for(i = 0; i< iter; i++){
		printf("%.4f\t\t", ber_at_iter[i]);
	}
	printf("\n");
	return EXIT_SUCCESS;

	}




