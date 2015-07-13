#define _GNU_SOURCE
#include <sched.h>
#include <stdio.h>
#include <stdbool.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include <pthread.h>
#include "flexible_gallager.h"

#define ITER_CHUNK 1000
#define NUM_CORES sysconf(_SC_NPROCESSORS_ONLN)

int n, m;	// size of parity check matrix
int iter;	// number of iterations per one signal
int max_iter; // maximum number of bp iterations to run
float target_ber; // target ber to achieve after bp
char filepath[256];
char filename[100];	// file name for parity check matrix
char **parity_matrix;	// parity check matrix
float ber = 0.01;
int num_cores = 1;
bool is_flexible = false;
bool is_sequential = false;
/** Used only for flexible codes */
int depth;	// depth of the flexible code
int level;	// level of depth you want to use for decoding
int **deg_pairs;	// degree pairs for each depth


void set_up_nodes(*var_nodes, *check_nodes){
	int i,j; // integers for for loops

	for (i = 0; i < _m; i++){
		for(j = 0; j < n ; j++){
			if (parity_matrix[i][j] == 1){
				var_nodes[j].deg++;
				check_nodes[i].deg++;
			}
		}
	}

	/* Initialize each var/check node 
	 * Memory alloccation according to each degree
	 */
	for (i = 0; i < n; i++){
		var_nodes[i].neighbors = malloc(var_nodes[i].deg * sizeof(int));
		var_nodes[i].neighbor_ports = malloc(var_nodes[i].deg * sizeof(int));
		var_nodes[i].send_msg = malloc(var_nodes[i].deg * sizeof(char));
		var_nodes[i].recv_msg = malloc(var_nodes[i].deg * sizeof(char));
		var_nodes[i].port_temp = 0;
	}

	for (i = 0; i < _m; i++){
		check_nodes[i].neighbors = malloc(check_nodes[i].deg * sizeof(int));
		check_nodes[i].neighbor_ports = malloc(check_nodes[i].deg * sizeof(int));
		check_nodes[i].send_msg = malloc(check_nodes[i].deg * sizeof(char));
		check_nodes[i].recv_msg = malloc(check_nodes[i].deg * sizeof(char));
		check_nodes[i].port_temp = 0;
	}

	/* Initialize connections between nodes
	 * Construct port-to-port connections
	 */
	for(i = 0; i < _m; i++){
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
}


int*  calc_threshold(*var_nodes, *check_nodes){
	double Al = 0;
	for (d = 0; d < level; d++){
		Al += (double)m_l[d] / (double)m_l_tot * pow((double)(1-2*ber), deg_pairs[d][1] - 1);
	}
	int b[n];
	for (i = 0; i < n; i++){
		b[i] = ceil((log((1-ber)/ber)/log((1+Al)/(1-Al)) + var_nodes[i].deg -1)/2);
		b[i] = b[i] > var_nodes[i].deg? var_nodes[i].deg : b[i];
	}

	return b;
}

void beleif_proapagation(*var_nodes, *check_nodes){
	int parity;
	char signal[n];
	int error_per_iter;
	float r;

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
		for(j = 0; j < _m; j++){
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
				if ((parity - var_nodes[i].recv_msg[c]) >=  b[i]){
					var_nodes[i].send_msg[c] = 1;
				} else if((parity - var_nodes[i].recv_msg[c]) <= (var_nodes[i].deg - b[i]) ){
					var_nodes[i].send_msg[c] = 0;
				} else{
					var_nodes[i].send_msg[c] = signal[i];
				}
			}

			if (parity > (var_nodes[i].deg/2)){
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


void *flexible_thr_func(void *arg){

	/* initialize var/check nodes */
	int i,j,l,t,d,c=0;
	int m_l[depth], m_l_tot = 0;
	struct node var_nodes[n];
	struct node check_nodes[m];
	int _m = m;

	/* Initialize degrees.
	   This steddp can be skipped if it's a regular code */

	for(d = 0; d < depth; d++){
		m_l[d] = n * deg_pairs[d][0]/deg_pairs[d][1];
	}

	for(d = 0; d < level; d++){
		m_l_tot += m_l[d];
	}

	if (~is_sequential) _m = m_l_tot;


	for (i = 0; i < n; i++){
		var_nodes[i].deg = 0;
	}
	for (j = 0; j < _m; j++){
		check_nodes[j].deg = 0;
	}

	for (i = 0; i < _m; i++){
		for(j = 0; j < n ; j++){
			if (parity_matrix[i][j] == 1){
				var_nodes[j].deg++;
				check_nodes[i].deg++;
			}
		}
	}


	/* Initialize each var/check node 
	 * Memory alloccation according to each degree
	 */
	for (i = 0; i < n; i++){
		var_nodes[i].neighbors = malloc(var_nodes[i].deg * sizeof(int));
		var_nodes[i].neighbor_ports = malloc(var_nodes[i].deg * sizeof(int));
		var_nodes[i].send_msg = malloc(var_nodes[i].deg * sizeof(char));
		var_nodes[i].recv_msg = malloc(var_nodes[i].deg * sizeof(char));
		var_nodes[i].port_temp = 0;
	}

	for (i = 0; i < _m; i++){
		check_nodes[i].neighbors = malloc(check_nodes[i].deg * sizeof(int));
		check_nodes[i].neighbor_ports = malloc(check_nodes[i].deg * sizeof(int));
		check_nodes[i].send_msg = malloc(check_nodes[i].deg * sizeof(char));
		check_nodes[i].recv_msg = malloc(check_nodes[i].deg * sizeof(char));
		check_nodes[i].port_temp = 0;
	}

	/* Initialize connections between nodes
	 * Construct port-to-port connections
	 */
	for(i = 0; i < _m; i++){
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
	double Al = 0;
	for (d = 0; d < level; d++){
		Al += (double)m_l[d] / (double)m_l_tot * pow((double)(1-2*ber), deg_pairs[d][1] - 1);
	}
	int b[n];
	for (i = 0; i < n; i++){
		b[i] = ceil((log((1-ber)/ber)/log((1+Al)/(1-Al)) + var_nodes[i].deg -1)/2);
		b[i] = b[i] > var_nodes[i].deg? var_nodes[i].deg : b[i];
	}

	printf("%d\n", b[0]);
	int iter_chunks = 0;
	FILE *write_fp;
	char write_filename[200];
	if (is_flexible) sprintf(write_filename, "./results_details/%s_%diter_%dsig_%.4fBER_%dlevel_CORE%d", filename,iter,sig_iter, ber, level, cpu);
	else sprintf(write_filename, "./results_details/%s_%diter_%dsig_%.4fBER_CORE%d", filename,iter,sig_iter, ber, cpu); 
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
		for(j = 0; j < _m; j++){
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
				if ((parity - var_nodes[i].recv_msg[c]) >=  b[i]){
					var_nodes[i].send_msg[c] = 1;
				} else if((parity - var_nodes[i].recv_msg[c]) <= (var_nodes[i].deg - b[i]) ){
					var_nodes[i].send_msg[c] = 0;
				} else{
					var_nodes[i].send_msg[c] = signal[i];
				}
			}

			if (parity > (var_nodes[i].deg/2)){
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
	bool n_flag = false, m_flag=false, i_flag=false, e_flag=false, M_flag = false, T_flag=false, f_flag=false, p_flag = false, l_flag=false, F_flag = false;
	while((opt=getopt(argc, argv, "n:m:f:i:e:t:c:l:Fh")) != -1){
		switch(opt){
			case 'n':
				n = atoi(optarg);
				n_flag = true;
				break;
			case 'm':
				m = atoi(optarg);
				m_flag = true;
				break;
			case 'i':
				iter = atoi(optarg);
				i_flag = true;
				break;
			case 'e':
				ber = atof(optarg);
				e_flag = true;
				break;
			case 'M':
				max_iter = atoi(optarg);
				M_flag = true;
				break;
			case 'T':
				target_ber = atof(optarg);
				T_flag = true;
				break;
			case 'p':
				sprintf(filepath, "%s", optarg);
				p_flag = true;
				break;
			case 'f':
				sprintf(filename, "%s", optarg);
				f_flag = true;
				break;
			case 'c':
				num_cores = (NUM_CORES > atoi(optarg))? atoi(optarg):NUM_CORES; 
				break;
			case 'F':
				is_flexible = true;
				F_flag = true;
				break;
			case 'l':
				level = atoi(optarg);
				l_flag = true;
				break;
			case 'h':
				fprintf(stderr, "Usage: %s [-n n] [-m m] [-i #iterations] [-e bit-error rate] [-t #signals] [-f filename] [-c #cores] [-F] [-l #levels] \n", argv[0]);
				exit(EXIT_FAILURE);
			default: /* '?' */
				fprintf(stderr, "Usage: %s [-n n] [-m m] [-i #iterations] [-e bit-error rate] [-t #signals] [-f filename] [-c #cores] [-F] [-l #levels] \n", argv[0]);
				exit(EXIT_FAILURE);
		}
	}

	/* Check if all necessary arguments were give */
	if (!(n_flag && m_flag && i_flag && e_flag && T_flag && f_flag)){
		fprintf(stderr, "Usage: %s [-n n] [-m m] [-i #iterations] [-e bit-error rate] [-t #signals] [-f filename] [-c #cores] [-F] [-l #levels] \n", argv[0]);
		exit(EXIT_FAILURE);
	}

	if(F_flag && !l_flag){
		fprintf(stderr, "You should specify how many levels you want to use with flexible(-F) option.");
		exit(EXIT_FAILURE);
	}

	/* Open input files */
	if(T_flag){
		strcat(filepath,"/");
		strcat(filepath, filename);	
	} else{
		strcpy(filepath, "./");
		strcat(filepath, filename);
	}
	/* Opening file for parity check matrix */
	fp = fopen(filepath, "r");
	if (fp == NULL){
		fprintf(stderr, "Error opening file:%s\n", strerror(errno));
		exit(EXIT_FAILURE);
	}
	/* Opening file for code structure */
	char code_str_filepath[256];
	strcpy(code_str_filepath, filepath);
	strcat(filepath, "_code_str");
	code_str_fp = fopen(code_str_filepath, "r");
	if (code_str_fp == NULL) {
		fprintf(stderr, "Error opening file:%s\n", strerror(errno));
		exit(EXIT_FAILURE);
	}
	
	/* Print out informations */
	printf("Testing (%d, %d) ", n, m);
	if (is_flexible) printf("flexible ");	
	printf("code using %d cores\n", num_cores);
	printf("Bit error rate : %.5f\n", ber);
	printf("# iterations in message passing algorithm: %d\n", iter);
	printf("# of signals to try: %d\n", sig_iter);
	printf("Reading parity matrix from the file: %s\n", filename);


	/* Read Code Structure */
	/* 1. Read the depth */
	fscanf(code_str_fp, "%d", &depth);
	/* 2. Read min/max degree */
	fscanf(code_str_fp, "%d %d", &min_deg, &max_deg);
	/* 3. Read length, degree distribution for each level */
	int m[depth];
	float rho[depth][max_deg-min_deg];
	for (i=0; i < depth; i++){
		fscanf(code_str_fp, "%d", &m[i]);
		for(j = min_deg; j <= max_deg; j++){
			fscanf(code_str_fp, "%.4f" &&rho[i][j]);
		}	
	}
	/* 4. Decide the code length to run */
	m = 0;
	for (i = 0; i < level; i++){
		m += m[i];
	}


	/* Reading a parity matrix */ 
	parity_matrix = malloc(m * sizeof(char*));
	for (i=0; i < m; i++){
		parity_matrix[i] = malloc(n*sizeof(char));
	}
	for (i=0; i < m; i++){
		fscanf(fp, "%s", parity_matrix[i]);
	}
	for (i=0; i< m; i++){
		for(j=0; j<n; j++){
			if (parity_matrix[i][j]=='1'){
				parity_matrix[i][j] = 1;
			} else{
				parity_matrix[i][j] = 0;
			}
		}
	}
	fclose(fp);

	/* Initialze the node structures & connections */


	/* Initialize the seed for random number generation */
	time_t init_time;
	srand((unsigned)time(&init_time));
	
	/* Calculating threshold */
	double Al = 0;
	for (d = 0; d < level; d++){	
		Al += (double)m_l[d] / (double)m_l_tot * pow((double)(1-2*ber), deg_pairs[d][1] - 1);
	}
	
	int b[n];
	for (i = 0; i < n; i++){
		b[i] = ceil((log((1-ber)/ber)/log((1+Al)/(1-Al)) + var_nodes[i].deg -1)/2);
		b[i] = b[i] > var_nodes[i].deg? var_nodes[i].deg : b[i];
	}
	

	struct node var_nodes[n];
	struct node check_nodes[m];
	
	set_up_nodes(var_nodes, check_nodes);

	int num_iter = 1/target_ber * 100;
	int t;

	for (t = 0; t < num_iter; t++){
		error_per_iter = beleif_propagation(var_nodes, check_nodes);

	}

	/*printf("Summing up the result..\n");
	float ber_at_iter[iter];
	char sum_filename[200];
	for (i = 0; i < iter; i++){
			ber_at_iter[i] = 0.0;
	}
	int sig_iter_per_core, val;
	for(ci = 0; ci < num_cores; ci++){
		sig_iter_per_core = thread_datas[ci].sig_iter_per_core;
		if (is_flexible) sprintf(sum_filename, "./results_details/%s_%diter_%dsig_%.4fBER_%dlevel_CORE%d", filename,iter,sig_iter_per_core, ber, level, ci);
	   	else sprintf(sum_filename, "./results_details/%s_%diter_%dsig_%.4fBER_CORE%d", filename,iter,sig_iter_per_core, ber, ci); 
		fp = fopen(sum_filename, "r");
		for(j = 0; j < sig_iter_per_core; j++){
			for(i = 0; i < iter; i++){
				fscanf(fp, "%d\t", &val);
				ber_at_iter[i] += (double)val/(double)n;
			}
		}
		fclose(fp);
	}
	for(i = 0; i < iter; i++){
		ber_at_iter[i] = (double)ber_at_iter[i]/(double)sig_iter;
	}

	if (is_flexible) sprintf(sum_filename, "./results/%s_%diter_%dsig_%.4fBER_%dlevel", filename,iter,sig_iter, ber, level);
	else sprintf(sum_filename, "./results/%s_%diter_%dsig_%.4fBER", filename,iter,sig_iter, ber); 

	fp = fopen(sum_filename, "w");
	printf("BER at iteration");
	for (i = 0; i < iter; i++){
		printf("  %d\t\t", i+1);
	}
	printf("\n                ");
	fprintf(fp, "%d\n", iter);
	for(i = 0; i< iter; i++){
		fprintf(fp, "%.10f\n", ber_at_iter[i]); 
		printf("%.4f\t\t", ber_at_iter[i]);
	}
	printf("\n");
	fclose(fp);*/
	return EXIT_SUCCESS;

	}




