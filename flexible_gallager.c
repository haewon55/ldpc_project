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
int max_iter; // maximum number of bp iterations to run
float target_ber; // target ber to achieve after bp
int min_deg, max_deg;
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


void set_up_nodes(struct node* var_nodes, struct node* check_nodes, int _m, float* rho){
	int i,j; // integers for for loops

	for (i = 0; i < _m; i++){
		check_nodes[i].deg = 0;
	}
	for (j = 0; j < n; j++){
		var_nodes[j].deg = 0;
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
	int b; // integer to store temporary value of b
	for (i = 0; i < n; i++){
		var_nodes[i].neighbors = malloc(var_nodes[i].deg * sizeof(int));
		var_nodes[i].neighbor_ports = malloc(var_nodes[i].deg * sizeof(int));
		var_nodes[i].send_msg = malloc(var_nodes[i].deg * sizeof(char));
		var_nodes[i].recv_msg = malloc(var_nodes[i].deg * sizeof(char));
		var_nodes[i].port_temp = 0;
		/* Calculate threshold b*/
		double Al = 0;
		for (j = 0; j < _m; j++){
			Al += pow((double)(1-2*ber), check_nodes[j].deg)/(double)_m;
		}
		b  = ceil((log((1-ber)/ber)/log((1+Al)/(1-Al)) + var_nodes[i].deg -1)/2);
		var_nodes[i].b =  b > var_nodes[i].deg? var_nodes[i].deg : b;
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


int* belief_propagation(struct node* var_nodes, struct node* check_nodes, int _m){
	int i, j, l; // integers for for loops
	int parity;
	char signal[n];
	float r;
	int *error_per_iter = malloc(max_iter * sizeof(int));
	for (i = 0; i < max_iter; i++){
		error_per_iter[i] = 0;
	}

	/* Signal initialization */
	for (i = 0; i < n; i++){
		r = (double)rand() / (double)RAND_MAX;
		if (r < ber){
			signal[i] = 1;
		} else {
			signal[i] = 0;
		}

		/* Setting first variable messages */
		for(j = 0; j < var_nodes[i].deg; j++){
			var_nodes[i].send_msg[j] = signal[i];
		}
	}

	/* Start message passing algorithm */
	for (l = 0; l < max_iter; l++){
		/* check node operations
			- receive message from variable nodes
			- perform parity checking
			- generate send messages
	   */
		for(i = 0; i < _m; i++){
			parity = 0;
			for (j = 0; j < check_nodes[j].deg; j++){
				check_nodes[i].recv_msg[j] = var_nodes[check_nodes[i].neighbors[j]].send_msg[check_nodes[i].neighbor_ports[j]];
				parity += check_nodes[i].recv_msg[j]; 
			}
			for (j = 0; j < check_nodes[i].deg; j++){
				check_nodes[i].send_msg[j] = (parity + check_nodes[i].recv_msg[j]) % 2;
			}
		}

		/* variable node operations
		    - receive message from check nodes
			- generate send messages
			- belief update
		*/
		for (i = 0; i < n; i++){
			parity = 0;
			for(j = 0; j < var_nodes[i].deg; j++){
				var_nodes[i].recv_msg[j] = check_nodes[var_nodes[i].neighbors[j]].send_msg[var_nodes[i].neighbor_ports[j]];
				parity += var_nodes[i].recv_msg[j];
			}
			for(j = 0; j < var_nodes[i].deg; j++){
				if ((parity - var_nodes[i].recv_msg[j]) >=  var_nodes[i].b){
					var_nodes[i].send_msg[j] = 1;
				} else if((parity - var_nodes[i].recv_msg[j]) <= (var_nodes[i].deg - var_nodes[i].b) ){
					var_nodes[i].send_msg[j] = 0;
				} else{
					var_nodes[i].send_msg[j] = signal[i];
				}
			}

			if (parity > (var_nodes[i].deg/2)){
				error_per_iter[l]++;
			}
		}
	}

	return error_per_iter;
}


int main(int argc, char *argv[]){
	int i, j;
	int opt;
	FILE *fp, *code_str_fp;
	bool n_flag = false, m_flag=false, e_flag=false, M_flag = false, T_flag=false, f_flag=false, p_flag = false, l_flag=false;
	while((opt=getopt(argc, argv, "n:m:e:M:T:p:f:l:h")) != -1){
		switch(opt){
			case 'n':
				n = atoi(optarg);
				n_flag = true;
				break;
			case 'm':
				m = atoi(optarg);
				m_flag = true;
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
			case 'l':
				level = atoi(optarg);
				l_flag = true;
				break;
			case 'h':
				fprintf(stderr,  argv[0]);
				exit(EXIT_FAILURE);
			default: /* '?' */
				fprintf(stderr, usage, argv[0]);
				exit(EXIT_FAILURE);
		}
	}

	/* Check if all necessary arguments were give */
	if (!(n_flag && m_flag && e_flag && M_flag && T_flag && f_flag)){
		fprintf(stderr, usage, argv[0]);
		exit(EXIT_FAILURE);
	}

	/* Open input files */
	if(p_flag){
		strcat(filepath,"/");
		strcat(filepath, filename);	
	} else{
		strcpy(filepath, "./");
		strcat(filepath, filename);
	}
	/* Opening file for parity check matrix */
	fp = fopen(filepath, "r");
	if (fp == NULL){
		fprintf(stderr, usage, "Error opening file:%s\n", strerror(errno));
		exit(EXIT_FAILURE);
	}
	/* Opening file for code structure */
	char code_str_filepath[256];
	strcpy(code_str_filepath, filepath);
	strcat(code_str_filepath, "_code_str");
	code_str_fp = fopen(code_str_filepath, "r");
	if (code_str_fp == NULL) {
		fprintf(stderr, "Error opening file:%s\n", code_str_filepath);
		exit(EXIT_FAILURE);
	}
	
	/* Print out informations */
	printf("Testing (%d, %d) ", n, m);
	printf("code using %d cores\n", num_cores);
	printf("Bit error rate : %.5f\n", ber);
	printf("Target BER to achieve: %.4f\n", target_ber);
	printf("Reading parity matrix from the file: %s\n", filename);


	/* Read Code Structure */
	/* 1. Read the depth */
	fscanf(code_str_fp, "%d", &depth);
	is_flexible = depth > 1;
	/* 2. Read min/max degree */
	fscanf(code_str_fp, "%d %d", &min_deg, &max_deg);
	/* 3. Read length, degree distribution for each level */
	int m[depth];
	float rho[depth][max_deg-min_deg+1];
	float temp = 0.0;
	for (i=0; i < depth; i++){
		fscanf(code_str_fp, "%d", &m[i]);
		for(j = min_deg; j <= max_deg; j++){
			fscanf(code_str_fp, "%.4f", &temp);
			rho[i][j-min_deg] = temp;
		}	
	}
	/* 4. Adjust the decoding length/rho */
	int m_adapt = 0;
	float rho_adapt[max_deg-min_deg+1];
	for(i = 0; i <= (max_deg-min_deg); i++){
		rho_adapt[i] = 0;
	}
	for (i = 0; i < level; i++){
		m_adapt += m[i];
		for(j = 0; j < (max_deg-min_deg); j++){
			rho_adapt[j] += rho[i][j];
		}
	}
	fclose(code_str_fp);

	/* Print out the read code structure */
	if(is_flexible){
		printf("Flexible Code with depth %d\n", depth);
		printf("Using level upto %d\n", level);
		printf("Decoding length = %d\n", m_adapt);
	}
	/* Read a parity check matrix */ 
	parity_matrix = malloc(m_adapt * sizeof(char*));
	for (i=0; i < m_adapt; i++){
		parity_matrix[i] = malloc(n*sizeof(char));
	}
	for (i=0; i < m_adapt; i++){
		fscanf(fp, "%s", parity_matrix[i]);
	}
	for (i=0; i< m_adapt; i++){
		for(j=0; j<n; j++){
			if (parity_matrix[i][j]=='1'){
				parity_matrix[i][j] = 1;
			} else{
				parity_matrix[i][j] = 0;
			}
		}
	}
	fclose(fp);

	printf("Done Reading Files\n");

	/* Initialze the node structures & connections */
	struct node var_nodes[n];
	struct node check_nodes[m_adapt];
	
	set_up_nodes(var_nodes, check_nodes, m_adapt, rho_adapt);

	/* Initialize the seed for random number generation */
	time_t init_time;
	srand((unsigned)time(&init_time));

	/* Start running belief propagation */
	int num_sig = 100/(target_ber *  n);
	int t;
	unsigned int error_per_iter[max_iter];
	int* result;
	for (i = 0; i < max_iter; i++){
		error_per_iter[i] = 0;
	}
	for (t = 0; t < num_sig; t++){
		if (t== 0){
			printf("Start running belief propagation\n");
			printf("Total %d signals will be tried\n", num_sig);
		} else if(t % ITER_CHUNK== 0){
			printf("%d Completed.\n", t);
		}
		result = belief_propagation(var_nodes, check_nodes, m_adapt);
		for (i = 0; i < max_iter; i++){
			error_per_iter[i] += result[i];
		}
	}

	/* Print out the result */
	printf("Summing up the result..\n");
	float ber_at_iter[max_iter];
	for (i = 0; i < max_iter; i++){
		ber_at_iter[i] = (float)error_per_iter[i]/ (float)(num_sig*n);
	}
	
	char sum_filename[256];
	if (is_flexible) sprintf(sum_filename, "%s_%dIter_%.4fTarget_%.4fBER_%dLevel", filename,max_iter,target_ber, ber, level);
	else sprintf(sum_filename, "%s_%dIter_%.4fTarget_%.4fBE", filename,max_iter,target_ber, ber); 

	fp = fopen(sum_filename, "w");
	printf("BER at iteration\n");
	for (i = 0; i < max_iter; i++){
		printf("%5d\t", i+1);
	}
	printf("\n");
	for(i = 0; i< max_iter; i++){
		fprintf(fp, "%.10f\n", ber_at_iter[i]); 
		printf("%.4f\t", ber_at_iter[i]);
	}
	printf("\n");
	fclose(fp);
	return EXIT_SUCCESS;

}


