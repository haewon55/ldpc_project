#define _GNU_SOURCE

#define ITER_CHUNK 1000
#define NUM_CORES sysconf(_SC_NPROCESSORS_ONLN)


/* Node structure */
struct node {
	int deg;
	int *neighbors;
	int *neighbor_ports;
	int port_temp;
	char *recv_msg;
	char *send_msg;
	int b; /* threshold b. only used for var nodes */
};


char* usage = "Usage: %s [-n n] [-m m] [-M max #iterations] [-e bit-error rate] [-T target ber] [-p filepath] [-f filename] [-l #levels] \n";
