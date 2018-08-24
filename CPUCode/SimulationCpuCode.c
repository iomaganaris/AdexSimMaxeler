#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Maxfiles.h"
#include "MaxSLiCInterface.h"
typedef unsigned long long timestamp;
//Time stamp function for measuring DFE runtime
static timestamp getTimestamp(){

	struct timeval now;
	gettimeofday(&now, NULL);
	return now.tv_usec + (timestamp)now.tv_sec * 1000000;

}
void print_neurons(double* x, int M){
	printf("vt\n");
	for(int i = 0; i < M*6; i+=6){
		printf("%.8e\n",x[i]);
	}
	printf("vm\n");
	for(int i = 1; i < M*6; i+=6){
		printf("%.8e\n",x[i]);
	}
	printf("I\n");
	for(int i = 2; i < M*6; i+=6){
		printf("%.8e\n",x[i]);
	}
	printf("x\n");
	for(int i = 3; i < M*6; i+=6){
		printf("%.8e\n",x[i]);
	}
	printf("Spike\n");
	for(int i = 4; i < M*6; i+=6){
		printf("%lf\n",x[i]);
	}
}
void print_synapses(double* syn, int N_S, int N_T){
	/*printf("conn\n");
	for(int i = 0; i < N_S*12; i+=12){
		for(int j = 0; j < N_T; j++){
			printf("%d, ", syn[j*12*N_S+i]);
			//if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		printf("\n");
	}*/
	/*printf("test\n");
	for(int i =0; i < N_S; i++){
				for(int j = 3; j < N_T*12; j+=12){
					printf("%f \n", syn[i*N_T*12+j]);
					//if((i*N_T+j+1)%4 == 0) printf("\n");
				}
				//printf("\n");
			}*/
	/*for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++){
			//printf("conn= %d, w= %lf, FFp= %lf, FBp= %lf, FBn= %lf, R= %lf, u= %lf, U= %lf, A= %lf, lastup= %lf, target_I= %lf\n",syn[i][j].conn,syn[i][j].w,syn[i][j].FFp,syn[i][j].FBp,syn[i][j].FBn,syn[i][j].R,syn[i][j].u,syn[i][j].U,syn[i][j].A,syn[i][j].lastupdate,syn[i][j].target_I);
			printf("%d ", syn[i][j].conn);
			//if((i*N_T+j)%4 == 0) printf("\n");
		}
		printf("\n");
	}*/
	/*printf("\nw\n");
	for(int i = 1; i < N_S*12; i+=12){
		for(int j = 0; j < N_T; j++){
			printf("%.8e, ", syn[j*12*N_S+i]);
			//if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		printf("\n");
	}
	printf("\nFFp\n");
	for(int i = 2; i < N_S*12; i+=12){
		for(int j = 0; j < N_T; j++){
			printf("%.8e, ", syn[j*12*N_S+i]);
			//if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		printf("\n");
	}
	printf("\nFBp\n");
	for(int i = 3; i < N_S*12; i+=12){
		for(int j = 0; j < N_T; j++){
			printf("%.8e, ", syn[j*12*N_S+i]);
			//if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		printf("\n");
	}
	printf("\nFBn\n");
	for(int i = 4; i < N_S*12; i+=12){
			for(int j = 0; j < N_T; j++){
				printf("%.8e, ", syn[j*12*N_S+i]);
				//if((i*N_T+j+1)%4 == 0) printf("\n");
			}
			printf("\n");
		}
	printf("\nR\n");
	for(int i = 5; i < N_S*12; i+=12){
			for(int j = 0; j < N_T; j++){
				printf("%.8e, ", syn[j*12*N_S+i]);
				//if((i*N_T+j+1)%4 == 0) printf("\n");
			}
			printf("\n");
		}*/
	printf("\nu\n");
	for(int i = 6; i < N_S*12; i+=12){
			for(int j = 0; j < N_T; j++){
				printf("%.8e, ", syn[j*12*N_S+i]);
				//if((i*N_T+j+1)%4 == 0) printf("\n");
			}
			printf("\n");
		}
	/*printf("\nU\n");
	for(int i = 7; i < N_S*12; i+=12){
			for(int j = 0; j < N_T; j++){
				printf("%.8e, ", syn[j*12*N_S+i]);
				//if((i*N_T+j+1)%4 == 0) printf("\n");
			}
			printf("\n");
		}*/
	printf("\nA\n");
	for(int i = 8; i < N_S*12; i+=12){
			for(int j = 0; j < N_T; j++){
				printf("%.8e, ", syn[j*12*N_S+i]);
				//if((i*N_T+j+1)%4 == 0) printf("\n");
			}
			printf("\n");
		}
	/*printf("\nlastupdate\n");
	for(int i = 9; i < N_S*12; i+=12){
			for(int j = 0; j < N_T; j++){
				printf("%f, ", syn[j*12*N_S+i]);
				//if((i*N_T+j+1)%4 == 0) printf("\n");
			}
			printf("\n");
		}
	printf("\ntarget_I\n");
	for(int i = 10; i < N_S*12; i+=12){
			for(int j = 0; j < N_T; j++){
				printf("%.8e, ", syn[j*12*N_S+i]);
				//if((i*N_T+j+1)%4 == 0) printf("\n");
			}
			printf("\n");
		}
	printf("\n");*/
}

int main(void)
{
	timestamp t0,t1;
	const int size = 384;
	uint32_t N = 384;
	uint32_t M = 384;
	uint32_t steps = 3;

	uint32_t N_S = 0;//100;
	uint32_t N_Group_S = M;
	uint32_t N_Group_T = M;

	double input1_pos = 25;
	double input2_pos = 75;
	double rad = 5;

	int size_d = size * 6;
	int sizeBytes_d_var = size_d * sizeof(double);
	//int sizeBytes_d_var = sizeBytes_d * 6;
	int sizeBytes_i = size * sizeof(int);
	double *x = malloc(sizeBytes_d_var);
	//double *y = malloc(sizeBytes_d);
	//double *s = malloc(sizeBytes_d);
	double scalar = 3.0;
	/*double *vm = malloc(sizeBytes_d);
	double *vt = malloc(sizeBytes_d);
	double *I = malloc(sizeBytes_d);
	double *x = malloc(sizeBytes_d);
	double *Spike = malloc(sizeBytes_i);*/


	// TODO Generate input data
	/*for(int i = 0; i<size; ++i) {
		//x[i] = i/10.0;
		y[i] = i/20.0;
	}*/


	int adex_param_size = 16;	// word alligned size of adex parameters array
	int adex_param_array_size = adex_param_size*sizeof(double);	// adex params array bytes
	double *adex_params = malloc(adex_param_array_size);
	//double *out = malloc(adex_param_array_size);
	adex_params[0] = 281*1e-12; // C
	adex_params[1] = 30*1e-9; // gL
	adex_params[2] = (281*1e-12) / (30*1e-9); // taum
	adex_params[3] = -70.6*1e-3; // EL
	adex_params[4] = 2*1e-3; // DeltaT
	adex_params[5] = -50.4*1e-3; // vti
	adex_params[6] = -45*1e-3; // vtrest
	adex_params[7] = 18*1e-3; // VTmax
	adex_params[8] = 50*1e-3; // tauvt
	adex_params[9] = 144*1e-3; // tauw
	adex_params[10] = 4*1e-9; // c
	adex_params[11] = 0.0805*1e-9; // b
	adex_params[12] = -70.6*1e-3; // Vr
	adex_params[13] = 1*1e-3; // defaultclock_dt

	int stdp_param_size = 32;
	int stdp_param_array_size = stdp_param_size*sizeof(double);
	double *stdp_params = malloc(stdp_param_array_size);
	stdp_params[0] = 0 * 1e-3;  //mV Ee
	stdp_params[1] = 2 * 1e-3;	//ms taue
	stdp_params[2] = 50; //Hz Fon
	stdp_params[3] = 3; //Hz Foff
	stdp_params[4] = 100 * 1e-10; // s
	stdp_params[5] = 2.0; 	//Amax
	stdp_params[6] = 0; 	//Amin
	stdp_params[7] = 0.1; 	//Ainit
	stdp_params[8] = 1.0;	//Umax
	stdp_params[9] = 0;		//Umin
	stdp_params[10] = 0.1;	//Uinit
	stdp_params[11] = 0;	//dFBn
	stdp_params[12] = 0;	//dFBp
	stdp_params[13] = 0;	//dFFp
	stdp_params[14] = 50 * 1e-3;	//ms tau_u
	stdp_params[15] = 200 * 1e-3;	//ms tau_r
	stdp_params[16] = 0.1771;		//AFBn
	stdp_params[17] = 0.0327 * 1e3 * 1e-3;	//ms tau_FBn
	stdp_params[18] = 0.1548;				//AFBp
	stdp_params[19] = 0.2302 * 1e3 * 1e-3;	//ms tau_FBp
	stdp_params[20] = 0.0618;				//AFFp
	stdp_params[21] = 0.0666 * 1e3 * 1e-3;	//ms tauFFp
	stdp_params[22] =  0.15;				//etaU
	stdp_params[23] = 0.15;					//etaA

	// Initialization of Neuron Variables
	for(uint32_t i = 0; i<M*6; i+=6) {
		//vt[i] = adex_params[6];
		//vm[i] = adex_params[3];
		//I[i] = 0;
		//x[i] = 0;
		//Spike[i] = 0;
		x[i] = adex_params[6];
		x[i+1] = adex_params[3];
		x[i+2] = 0;//i;
		x[i+3] = 0;
		x[i+4] = 1;
	}

	// Initialization of Synapses Variables
	//Synapse *syn[N_S+N_Group_S];
	int syn_size = (N_S+N_Group_S) * N_Group_T * 12;
	int syn_sizeBytes = syn_size * sizeof(double);
	double* syn = malloc(syn_sizeBytes);
	//double syn[384*384*12*sizeof(double)];
	//for(int i = 0; i < N_S+N_Group_S; i++) syn[i] = (double*)malloc(sizeof(double) * (N_Group_T));
	int con = 1;
	//Nx(M*12)
	/*for(int i = 0; i < N_Group_S+N_S; i++)
		for(int j = 0; j < N_Group_T*12; j+=12)
			syn[i*N_Group_T*12+j] = 1;
	// Initialization of Synapses for Neurons
	int init_const = 0;
	for(int i = 0; i < N_Group_S; i++){
		for(int j = 0; j < N_Group_T*12; j+=12){
			if (syn[i*N_Group_T*12+j]) {
				//syn[i][j].conn = 1;	// all connected
				//Connectivity, initialization now happens only in connected synapses.
				syn[i*N_Group_T*12+j+3] = 0;
				syn[i*N_Group_T*12+j+4] = 0;
				syn[i*N_Group_T*12+j+5] = 1;
				syn[i*N_Group_T*12+j+7] = exp(-(((pow((init_const+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(stdp_params[8]-stdp_params[9])+stdp_params[9];	// takes time
				syn[i*N_Group_T*12+j+8] = exp(-(((pow((init_const+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(stdp_params[5]-stdp_params[6])+stdp_params[6];	// takes time
				init_const++;
				//syn[i][j].U = exp(-(((pow(((i)+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(Umax-Umin)+Umin;	// takes time
				//syn[i][j].A = exp(-(((pow(((i)+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(Amax-Amin)+Amin;	// takes time
			}
		}
	}
	// Initialization of Synapses for external input (bottom rows)
	for(int i = N_Group_S; i < N_Group_S+N_S; i++){
		for(int j = 0; j < N_Group_T*12; j+=12){
			if (syn[i*N_Group_T*12+j]) {
				//Connectivity, initialization now happens only in connected synapses.
				//syn[i][j].conn = 1;	// all connected
				syn[i*N_Group_T*12+j+3] = init_const;//0;
				syn[i*N_Group_T*12+j+4] = 0;
				syn[i*N_Group_T*12+j+5] = 1;
				//syn[i][j].U = exp(-(((pow((((i-N_Group_S)*N_Group_T/M+j)+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(Umax-Umin)+Umin;	// takes time
				//syn[i][j].A = exp(-(((pow((((i-N_Group_S)*N_Group_T/M+j)+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(Amax-Amin)+Amin;	// takes time
				syn[i*N_Group_T*12+j+7] = exp(-(((pow(((init_const)+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(stdp_params[8]-stdp_params[9])+stdp_params[9];	// takes time
				syn[i*N_Group_T*12+j+8] = exp(-(((pow(((init_const)+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(stdp_params[5]-stdp_params[6])+stdp_params[6];	// takes time
				init_const++;
			}
		}
	}*/

	//Mx(N*12)
	for(int i = 0; i < (N_Group_S+N_S)*12; i+=12)
		for(int j = 0; j < N_Group_T; j++)
			syn[i+j*(N_Group_S+N_S)*12] = 1;
	// Initialization of Synapses for Neurons
	int init_const = 0;
	for(int i = 0; i < N_Group_S*12; i+=12){
		for(int j = 0; j < N_Group_T; j++){
			if (syn[i+j*(N_Group_S+N_S)*12]) {
				//syn[i][j].conn = 1;	// all connected
				//Connectivity, initialization now happens only in connected synapses.
				syn[i+j*(N_Group_S+N_S)*12+2] = init_const;//0;
				syn[i+j*(N_Group_S+N_S)*12+3] = init_const;//0;
				syn[i+j*(N_Group_S+N_S)*12+4] = init_const;//0;
				syn[i+j*(N_Group_S+N_S)*12+5] = init_const;//1;
				syn[i+j*(N_Group_S+N_S)*12+7] = init_const;//exp(-(((pow((init_const+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(stdp_params[8]-stdp_params[9])+stdp_params[9];	// takes time
				syn[i+j*(N_Group_S+N_S)*12+8] = init_const;//exp(-(((pow((init_const+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(stdp_params[5]-stdp_params[6])+stdp_params[6];	// takes time
				init_const++;
				//syn[i][j].U = exp(-(((pow(((i)+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(Umax-Umin)+Umin;	// takes time
				//syn[i][j].A = exp(-(((pow(((i)+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(Amax-Amin)+Amin;	// takes time
			}
		}
	}
	// Initialization of Synapses for external input (bottom rows)
	for(int i = N_Group_S*12; i < (N_Group_S+N_S)*12; i+=12){
		for(int j = 0; j < N_Group_T; j++){
			if (syn[i+j*(N_Group_S+N_S)*12]) {
				//Connectivity, initialization now happens only in connected synapses.
				//syn[i][j].conn = 1;	// all connected
				syn[i+j*(N_Group_S+N_S)*12+2] = init_const;//0;
				syn[i+j*(N_Group_S+N_S)*12+3] = init_const;//0;//init_const;//0;
				syn[i+j*(N_Group_S+N_S)*12+4] = init_const;//0;
				syn[i+j*(N_Group_S+N_S)*12+5] = init_const;//1;
				//syn[i][j].U = exp(-(((pow((((i-N_Group_S)*N_Group_T/M+j)+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(Umax-Umin)+Umin;	// takes time
				//syn[i][j].A = exp(-(((pow((((i-N_Group_S)*N_Group_T/M+j)+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(Amax-Amin)+Amin;	// takes time
				syn[i+j*(N_Group_S+N_S)*12+7] = init_const;//exp(-(((pow(((init_const)+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(stdp_params[8]-stdp_params[9])+stdp_params[9];	// takes time
				syn[i+j*(N_Group_S+N_S)*12+8] = init_const;//exp(-(((pow(((init_const)+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(stdp_params[5]-stdp_params[6])+stdp_params[6];	// takes time
				init_const++;
			}
		}
	}
	// Initialization of Input Neurons Spikes
	int InputSpikes_size = N_S * steps;
	int InputSpikes_size_bytes = InputSpikes_size * sizeof(int);
	int* InputSpikes = malloc(InputSpikes_size_bytes);
	for(int i = 0; i < steps; i++){
		for(int j = 0; j < N_S; j++){
			InputSpikes[i*N_S+j] = 1;
		}
	}

	//print_synapses(syn, N_S+N_Group_S, N_Group_T);
	////print_synapses(syn, N_Group_T, N_S+N_Group_S);

	printf("Writing to LMem.\n");
	Simulation_writeLMem(size_d, 0, x);
	Simulation_writeLMem(syn_size, size_d, syn);
	Simulation_writeLMem(InputSpikes_size, size_d + syn_size, InputSpikes);
	//UpdateSynapses_post_writeLMem(InputSpikes_size/sizeof(int), sizeBytes_d_var/sizeof(double) + syn_size/sizeof(double), InputSpikes);

	/*UpdateNeurons_writeLMem(0, sizeBytes_d, vm);
	UpdateNeurons_writeLMem(sizeBytes_d, sizeBytes_d, vt);
	UpdateNeurons_writeLMem(sizeBytes_d * 2, sizeBytes_d, I);
	UpdateNeurons_writeLMem(sizeBytes_d * 3, sizeBytes_d, x);
	UpdateNeurons_writeLMem(sizeBytes_d * 4, sizeBytes_i, Spike);*/

	/*printf("Adex Parameters\n");
	for(int i = 0; i<14; i++){
		printf("adex_params[%d] = %.8e\n",i,adex_params[i]);
	}*/

	printf("Running on DFE.\n");	
	fflush(stdout);

	long BurstLengthInBytes = 384;
	//UpdateNeurons(M, N, adex_param_size, stdp_param_size, steps, adex_params, stdp_params);
	//UpdateNeurons(scalar, M, N, M, adex_param_size, stdp_param_size, steps, adex_params, stdp_params, y ,s);
	t0 = getTimestamp();
	//printf("loopLength = %d",UpdateSynapses_post_get_UpdateSynapses_postKernel_loopLength());
	Simulation(N_Group_S, N_Group_T, N_S, adex_param_size, stdp_param_size, steps, BurstLengthInBytes, adex_params, stdp_params);
	t1 = getTimestamp();
	printf("DFE Runtime = %9.7f seconds\n", (t1-t0)/1000000.0 );
	fflush(stdout);
	Simulation_readLMem(size_d, 0, x);
	Simulation_readLMem(syn_size, size_d, syn);
	/*UpdateNeurons_readLMem(0, sizeBytes_d, vm);
	UpdateNeurons_readLMem(sizeBytes_d, sizeBytes_d, vt);
	UpdateNeurons_readLMem(sizeBytes_d * 2, sizeBytes_d, I);
	UpdateNeurons_readLMem(sizeBytes_d * 3, sizeBytes_d, x);
	UpdateNeurons_readLMem(sizeBytes_d * 4, sizeBytes_i, Spike);*/

	print_neurons(x, N_Group_T);
	//print_synapses(syn, N_S+N_Group_S, N_Group_T);

	// TODO Use result data
	/*for(int i=0; i<size; ++i)
		if (s[i] != x[i] + y[i] + scalar)
			return 1;*/

	/*printf("vt\n");
	for(int i = 0; i < M*6; i+=6){
		printf("%.8e\n",x[i]);
	}
	printf("vm\n");
	for(int i = 1; i < M*6; i+=6){
		printf("%.8e\n",x[i]);
	}*/
	/*printf("I\n");
	for(int i = 2; i < M*6; i+=6){
		printf("%.8e\n",x[i]);
	}*/
	//print_synapses(syn, N_Group_T, N_S+N_Group_S);
	/*printf("x\n");
	for(int i = 3; i < M*6; i+=6){
		printf("%.8e\n",x[i]);
	}
	printf("Spike\n");
	for(int i = 4; i < M*6; i+=6){
		printf("%.8e\n",x[i]);
	}*/

	printf("Done.\n");
	return 0;
}
