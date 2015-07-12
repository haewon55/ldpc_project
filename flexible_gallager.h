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
};


