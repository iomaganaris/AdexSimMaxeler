#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Maxfiles.h"
#include "MaxSLiCInterface.h"

#define NxM
//#define MxM
//#define NxMxM	//not working probably
#define allconnected

double defaultclock_dt = 1*1e-3;	//ms

//double taum_S = 10 * 1e-3; //ms
double Ee = 0 * 1e-3; //mV
double taue = 2 * 1e-3; //ms
double Fon = 50; //Hz
double Foff = 3; //Hz

#ifdef NxM
	double s = 100 * 1e-10;//100*1e-10;
#endif
#ifdef MxM
	double s = 500000;	//for testing
#endif
double Amax = 2.0;
double Amin = 0;
double Ainit = 0.1;
double Umax = 1.0;
double Umin = 0;
double Uinit = 0.1;

double dFBn = 0;
double dFBp = 0;
double dFFp = 0;

//#Short-term plasticity params
double tau_u = 50 * 1e-3;	//ms
double tau_r = 200 * 1e-3;	//ms

//#prepostSTDP params: AFBn tau_FBn AFBp tau_FBp AFFp tau_FFp
//double params[6] = {0.1771,    0.0327,    0.1548,    0.2302,    0.0618,    0.0666};
double AFBn = 0.1771;
double tau_FBn = 0.0327 * 1e3 * 1e-3;	//ms
double AFBp = 0.1548;
double tau_FBp = 0.2302 * 1e3 * 1e-3;	//ms
double AFFp = 0.0618;
double tau_FFp = 0.0666 * 1e3 * 1e-3;	//ms
//#etaU = 0.35
double etaU = 0.15;
double etaA = 0.15;
//#etaA = 0.35

//# Adex Parameters
double C = 281*1e-12;	//pF
double gL = 30*1e-9; //nS
double taum = 281*1e-12 / 30*1e-9;	// C/gL	// double initilization of taum(?)
double EL = -70.6*1e-3;	//mV
double DeltaT = 2*1e-3;	//mV
double vti = -50.4*1e-3;	//mV
//#vtrest = vti + 5 * DeltaT
double vtrest = -45*1e-3;	//mV
double VTmax = 18*1e-3;	//mV
double tauvt = 50*1e-3;	//ms

double tauw = 144*1e-3;	//ms
double c = 4*1e-9;	//ns
double b = 0.0805*1e-9;	//nA
double Vr = -70.6*1e-3;	//mV

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
void print_synapses(double* syn, int N_S, int N_S_pad, int N_Group_S, int N_Group_S_pad, int N_Group_T, int N_Group_T_pad){
	printf("conn\n");
	for(int i = 0; i < N_Group_S*12; i+=12){
		for(int j = 0; j < N_Group_T; j++){
			printf("%d, ", (int)syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
			//if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		printf("\n");
	}
	for(int i = N_Group_S_pad*12; i < (N_Group_S_pad+N_S)*12; i+=12){
		for(int j = 0; j < N_Group_T; j++){
			printf("%f, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
			//if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		printf("\n");
	}
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
	printf("\nw\n");
	for(int i = 1; i < N_Group_S*12; i+=12){
		for(int j = 0; j < N_Group_T; j++){
			printf("%.8e, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
			//if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		printf("\n");
	}
	for(int i = N_Group_S_pad*12+1; i < (N_Group_S_pad+N_S)*12; i+=12){
		for(int j = 0; j < N_Group_T; j++){
			printf("%.8e, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
			//if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		printf("\n");
	}
	printf("\nFFp\n");
	for(int i = 2; i < N_Group_S*12; i+=12){
		for(int j = 0; j < N_Group_T; j++){
			printf("%.8e, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
			//if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		printf("\n");
	}
	for(int i = N_Group_S_pad*12+2; i < (N_Group_S_pad+N_S)*12; i+=12){
		for(int j = 0; j < N_Group_T; j++){
			printf("%.8e, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
			//if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		printf("\n");
	}
	printf("\nFBp\n");
	for(int i = 3; i < N_Group_S*12; i+=12){
		for(int j = 0; j < N_Group_T; j++){
			printf("%.8e, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
			//if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		printf("\n");
	}
	for(int i = N_Group_S_pad*12+3; i < (N_Group_S_pad+N_S)*12; i+=12){
		for(int j = 0; j < N_Group_T; j++){
			printf("%.8e, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
			//if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		printf("\n");
	}
	printf("\nFBn\n");
	for(int i = 4; i < N_Group_S*12; i+=12){
		for(int j = 0; j < N_Group_T; j++){
			printf("%.8e, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
			//if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		printf("\n");
	}
	for(int i = N_Group_S_pad*12+4; i < (N_Group_S_pad+N_S)*12; i+=12){
		for(int j = 0; j < N_Group_T; j++){
			printf("%.8e, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
			//if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		printf("\n");
	}
	printf("\nR\n");
	for(int i = 5; i < N_Group_S*12; i+=12){
		for(int j = 0; j < N_Group_T; j++){
			printf("%.8e, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
			//if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		printf("\n");
	}
	for(int i = N_Group_S_pad*12+5; i < (N_Group_S_pad+N_S)*12; i+=12){
		for(int j = 0; j < N_Group_T; j++){
			printf("%.8e, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
			//if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		printf("\n");
	}
	printf("\nu\n");
	for(int i = 6; i < N_Group_S*12; i+=12){
		for(int j = 0; j < N_Group_T; j++){
			printf("%.8e, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
			//if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		printf("\n");
	}
	for(int i = N_Group_S_pad*12+6; i < (N_Group_S_pad+N_S)*12; i+=12){
		for(int j = 0; j < N_Group_T; j++){
			printf("%.8e, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
			//if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		printf("\n");
	}
	printf("\nU\n");
	for(int i = 7; i < N_Group_S*12; i+=12){
		for(int j = 0; j < N_Group_T; j++){
			printf("%.8e, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
			//if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		printf("\n");
	}
	for(int i = N_Group_S_pad*12+7; i < (N_Group_S_pad+N_S)*12; i+=12){
		for(int j = 0; j < N_Group_T; j++){
			printf("%.8e, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
			//if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		printf("\n");
	}
	printf("\nA\n");
	for(int i = 8; i < N_Group_S*12; i+=12){
		for(int j = 0; j < N_Group_T; j++){
			printf("%.8e, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
			//if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		printf("\n");
	}
	for(int i = N_Group_S_pad*12+8; i < (N_Group_S_pad+N_S)*12; i+=12){
		for(int j = 0; j < N_Group_T; j++){
			printf("%.8e, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
			//if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		printf("\n");
	}
	printf("\nlastupdate\n");
	for(int i = 9; i < N_Group_S*12; i+=12){
		for(int j = 0; j < N_Group_T; j++){
			printf("%f, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
			//if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		printf("\n");
	}
	for(int i = N_Group_S_pad*12+9; i < (N_Group_S_pad+N_S)*12; i+=12){
		for(int j = 0; j < N_Group_T; j++){
			printf("%f, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
			//if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		printf("\n");
	}
	printf("\ntarget_I\n");
	for(int i = 10; i < N_Group_S*12; i+=12){
		for(int j = 0; j < N_Group_T; j++){
			printf("%.8e, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
			//if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		printf("\n");
	}
	for(int i = N_Group_S_pad*12+10; i < (N_Group_S_pad+N_S)*12; i+=12){
		for(int j = 0; j < N_Group_T; j++){
			printf("%.8e, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
			//if((i*N_T+j+1)%4 == 0) printf("\n");
		}
		printf("\n");
	}
	printf("\n");
	fflush(stdout);
}

int main(void)
{
	timestamp t0,t1;
	int BurstLengthBytes = 384;
		FILE *f = fopen("Spikes.txt", "w");
		FILE *g = fopen("Neurons_I.txt", "w");
		FILE *h = fopen("array_A.txt", "w");
	//neofytou connectivity testing
		FILE *in = fopen("connections.txt", "r");
		int nruns = 1;

		for(int nrun = 0; nrun < nruns; nrun++){
		double realtime = 0;
		double stime = 0.004; //second
		double stime2 = 50; //second

		double resolution_export = 10 * 1e-3; //every x ms


		int N = 100;
		#ifdef MxM
			N = 0;
		#endif
		int M = 2;

		uint32_t N_S;//100;
		uint32_t N_Group_S;
		uint32_t N_Group_T;

		#ifdef NxM
			N_S = N;//100;
			N_Group_S = 0;
			N_Group_T = M;	//for the simulation we have, normally is N
			#endif

		#ifdef MxM
			N_S = 0;//100;
			N_Group_S = M;
			N_Group_T = M;	//for the simulation we have, normally is N
		#endif


		uint32_t N_S_pad = N_S, N_Group_S_pad = N_Group_S, N_Group_T_pad = N_Group_T;
		if(N_S%BurstLengthBytes != 0) N_S_pad = (N_S/BurstLengthBytes)*BurstLengthBytes+BurstLengthBytes;
		if(N_Group_S%BurstLengthBytes != 0) N_Group_S_pad = (N_Group_S/BurstLengthBytes)*BurstLengthBytes+BurstLengthBytes;
		if(N_Group_T%BurstLengthBytes != 0) N_Group_T_pad = (N_Group_T/BurstLengthBytes)*BurstLengthBytes+BurstLengthBytes;

		printf("N_S_pad = %d, N_Group_S_pad = %d, N_Group_T_pad = %d\n", N_S_pad, N_Group_S_pad, N_Group_T_pad);

		double input1_pos = 25;
		double input2_pos = 75;
		double rad = 5;

		//Define input 1
		double *F_input1;
		F_input1 = (double*)malloc(sizeof(double)*N_S);
		for(int i = 0; i < N_S; i++){
			F_input1[i] = Foff;			//maybe is not needed
			F_input1[i] = exp(-(pow((i+1)-input1_pos,2)/(2.0*pow(rad,2))))*(Fon-Foff)+Foff; //Define gaussian input
		}

		//Define input 2
		double *F_input2;
		F_input2 = (double*)malloc(sizeof(double)*N_S);
		for(int i = 0; i < N_S; i++){
			F_input2[i] = Foff;			//maybe is not needed
			F_input2[i] = exp(-(pow((i+1)-input2_pos,2)/(2.0*pow(rad,2))))*(Fon-Foff)+Foff; //Define gaussian input
		}

		/*Poisson *input;
		input = (Poisson*)malloc(sizeof(Poisson)*N_S);
		for(int i = 0; i<N_S; i++){				// Initialization of Poisson Neurons
			input[i].GaussArray = F_input1[i];
			input[i].Spike = 0;
		}*/

		// ******** ADEX DE Parameters Initialization *********** //
		int adex_param_size = 16;	// word alligned size of adex parameters array
		int adex_param_array_size = adex_param_size*sizeof(double);	// adex params array bytes
		double *adex_params = malloc(adex_param_array_size);
		//double *out = malloc(adex_param_array_size);
		adex_params[0] = C;
		adex_params[1] = gL;
		adex_params[2] = taum;
		adex_params[3] = EL;
		adex_params[4] = DeltaT;
		adex_params[5] = vti;
		adex_params[6] = vtrest;
		adex_params[7] = VTmax;
		adex_params[8] = tauvt;
		adex_params[9] = tauw;
		adex_params[10] = c;
		adex_params[11] = b;
		adex_params[12] = Vr;
		adex_params[13] = defaultclock_dt;

		// *********** STDP DE Parameters Initialization ************ //
		int stdp_param_size = 32;
		int stdp_param_array_size = stdp_param_size*sizeof(double);
		double *stdp_params = malloc(stdp_param_array_size);
		stdp_params[0] = Ee;
		stdp_params[1] = taue;
		stdp_params[2] = Fon;
		stdp_params[3] = Foff;
		stdp_params[4] = s;
		stdp_params[5] = Amax;
		stdp_params[6] = Amin;
		stdp_params[7] = Ainit;
		stdp_params[8] = Umax;
		stdp_params[9] = Umin;
		stdp_params[10] = Uinit;
		stdp_params[11] = dFBn;
		stdp_params[12] = dFBp;
		stdp_params[13] = dFFp;
		stdp_params[14] = tau_u;
		stdp_params[15] = tau_r;
		stdp_params[16] = AFBn;
		stdp_params[17] = tau_FBn;
		stdp_params[18] = AFBp;
		stdp_params[19] = tau_FBp;
		stdp_params[20] = AFFp;
		stdp_params[21] = tau_FFp;
		stdp_params[22] = etaU;
		stdp_params[23] = etaA;

		// ********* Neuron Initialization ******** //
		int size_d = N_Group_T_pad * 6;
		int sizeBytes_d_var = size_d * sizeof(double);
		double *x = malloc(sizeBytes_d_var);
		for(uint32_t i = 0; i<M*6; i+=6) {
			//vt[i] = adex_params[6];
			//vm[i] = adex_params[3];
			//I[i] = 0;
			//x[i] = 0;
			//Spike[i] = 0;
			x[i] = vtrest;
			//x[i+1] = adex_params[3];
			#ifdef NxM
				x[i+1] = EL;
			#endif
			#ifdef MxM
				x[i+1] = vtrest + 0.005;//EL;
			#endif
			x[i+2] = 0;//i;
			x[i+3] = 0;
			x[i+4] = 0;	// Spike
		}

		// ******** Synapses Initialization ******** //
		int syn_size = (N_S_pad+N_Group_S_pad) * N_Group_T_pad * 12;
		int syn_sizeBytes = syn_size * sizeof(double);
		double* syn = malloc(syn_sizeBytes);
		#ifndef allconnected
			int con = 0;

			for(int i = 0; i < N_Group_S; i++){
				for(int j = 0; j < N_Group_T; j++){
					fscanf (in, "%d", &con);
					if (con == 1) syn[i+j*(N_Group_S+N_S)*12] = 1;
					else syn[i+j*(N_Group_S+N_S)*12] = 0;
				}
			}

			for(int i = N_Group_S; i < N_Group_S+N_S; i++){
				for(int j = 0; j < N_Group_T; j++){
					fscanf (in, "%d", &con);
					if (con == 1) syn[i+j*(N_Group_S+N_S)*12] = 1;
					else syn[i+j*(N_Group_S+N_S)*12] = 0;
				}
			}
		#else
			for(int j = 0; j < N_Group_T_pad; j++)
				for(int i = 0; i < N_Group_S_pad*12; i+=12)
					if((i < N_Group_S*12) & (j < N_Group_T)) syn[i+j*(N_Group_S_pad+N_S_pad)*12] = 1;
					else syn[i+j*(N_Group_S_pad+N_S_pad)*12] = 0;
			for(int j = 0; j < N_Group_T_pad; j++)
				for(int i = N_Group_S_pad*12; i < (N_Group_S_pad+N_S_pad)*12; i+=12)
					if((i-N_Group_S_pad*12 < N_S*12) & (j < N_Group_T)) syn[i+j*(N_Group_S_pad+N_S_pad)*12] = 1;
					else syn[i+j*(N_Group_S_pad+N_S_pad)*12] = 0;
		#endif
		int init_const = 0;
		for(int i = 0; i < N_Group_S_pad*12; i+=12){
			for(int j = 0; j < N_Group_T_pad; j++){
				if (syn[i+j*(N_Group_S_pad+N_S_pad)*12]) {
					//syn[i][j].conn = 1;	// all connected
					//Connectivity, initialization now happens only in connected synapses.
					//syn[i+j*(N_Group_S+N_S)*12+2] = init_const;//0;
					syn[i+j*(N_Group_S_pad+N_S_pad)*12+3] = init_const;//0;
					syn[i+j*(N_Group_S_pad+N_S_pad)*12+4] = init_const;//0;
					syn[i+j*(N_Group_S_pad+N_S_pad)*12+5] = init_const;//1;
					#ifdef MxM
					syn[i+j*(N_Group_S_pad+N_S_pad)*12+6] = 1;	// for testing
					#endif
					syn[i+j*(N_Group_S_pad+N_S_pad)*12+7] = init_const;//exp(-(((pow((init_const+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(stdp_params[8]-stdp_params[9])+stdp_params[9];	// takes time
					syn[i+j*(N_Group_S_pad+N_S_pad)*12+8] = init_const;//exp(-(((pow((init_const+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(stdp_params[5]-stdp_params[6])+stdp_params[6];	// takes time
					init_const++;
					//syn[i][j].U = exp(-(((pow(((i)+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(Umax-Umin)+Umin;	// takes time
					//syn[i][j].A = exp(-(((pow(((i)+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(Amax-Amin)+Amin;	// takes time
				}
			}
		}
		// Initialization of Synapses for external input (bottom rows)
		for(int i = N_Group_S_pad*12; i < (N_Group_S_pad+N_S_pad)*12; i+=12){
			for(int j = 0; j < N_Group_T_pad; j++){
				if (syn[i+j*(N_Group_S_pad+N_S_pad)*12]) {
					//Connectivity, initialization now happens only in connected synapses.
					//syn[i][j].conn = 1;	// all connected
					//syn[i+j*(N_Group_S+N_S)*12+2] = init_const;//0;
					syn[i+j*(N_Group_S_pad+N_S_pad)*12+3] = init_const;//0;//init_const;//0;
					syn[i+j*(N_Group_S_pad+N_S_pad)*12+4] = init_const;//0;
					syn[i+j*(N_Group_S_pad+N_S_pad)*12+5] = init_const;//1;
					//syn[i][j].U = exp(-(((pow((((i-N_Group_S)*N_Group_T/M+j)+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(Umax-Umin)+Umin;	// takes time
					//syn[i][j].A = exp(-(((pow((((i-N_Group_S)*N_Group_T/M+j)+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(Amax-Amin)+Amin;	// takes time
					syn[i+j*(N_Group_S_pad+N_S_pad)*12+7] = init_const;//exp(-(((pow(((init_const)+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(stdp_params[8]-stdp_params[9])+stdp_params[9];	// takes time
					syn[i+j*(N_Group_S_pad+N_S_pad)*12+8] = init_const;//exp(-(((pow(((init_const)+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(stdp_params[5]-stdp_params[6])+stdp_params[6];	// takes time
					init_const++;
				}
			}
		}

		uint32_t timesteps = stime/defaultclock_dt;
		printf("timesteps=%d\n",timesteps);

		// ********** Input Neurons Spikes Initialization ********** //
		int SpikeArray_size = N_S_pad * timesteps * sizeof(int);
		int* SpikeArray = malloc(SpikeArray_size);
		for(int i = 0; i < timesteps; i++){
			for(int j = 0; j < N_S_pad; j++){
				SpikeArray[i*N_S_pad+j] = i%2;
			}
		}

		srand(time(NULL));

		double init_time = 0, load_time = 0, unload_time = 0, memory_reads = 0, memory_writes = 0, DFE_time = 0;
		timestamp t0, t1;


		#ifdef NxM
			//if(t == 0 ) SpikeArray[89+N_Group_T] = 1;
			//else SpikeArray[89+N_Group_T] = 0;

			//if(t == 2 || t == 25) SpikeArray[73+N_Group_T] = 1;
			//else SpikeArray[73+N_Group_T] = 0;

			//if(t*defaultclock_dt == 0.001) SpikeArray[25+N_Group_S] = 1;
			//else SpikeArray[25+N_Group_S] = 0;

			/*if(t == 4) SpikeArray[23] = 1;
			else SpikeArray[23] = 0;*/
			//SpikeArray[4*N_S_pad+23] = 1;

			/*if(t == 5){
				SpikeArray[20] = 1;
				SpikeArray[45] = 1;

			}
			else {
				SpikeArray[20] = 0;
				SpikeArray[45] = 0;

			}*/
			//SpikeArray[5*N_S_pad+20] = 1;
			//SpikeArray[5*N_S_pad+45] = 1;
			//if(t == 6 ) SpikeArray[81+N_Group_T] = 1;
			//else SpikeArray[81+N_Group_T] = 0;

			/*if(t == 7 ) SpikeArray[7] = 1;
			else SpikeArray[7] = 0;*/
			//SpikeArray[7*N_S_pad+7] = 1;
			/*if(t == 7 || t == 11 || t == 13 || t == 20 || t == 26 || t == 28 || t == 31) SpikeArray[27] = 1;
			else SpikeArray[27] = 0;*/
			//SpikeArray[7*N_S_pad+27] = 1;
			//SpikeArray[11*N_S_pad+27] = 1;
			//SpikeArray[13*N_S_pad+27] = 1;
			//SpikeArray[20*N_S_pad+27] = 1;
			//SpikeArray[26*N_S_pad+27] = 1;
			//SpikeArray[28*N_S_pad+27] = 1;
			//SpikeArray[31*N_S_pad+27] = 1;

			/*if(t == 9 ) SpikeArray[25] = 1;
			else SpikeArray[25] = 0;*/
			//SpikeArray[9*N_S_pad+25] = 1;

			/*if(t == 11 ) SpikeArray[11] = 1;
			else SpikeArray[11] = 0;*/
			//SpikeArray[11*N_S_pad+11] = 1;

			/*if(t == 15 ) SpikeArray[29] = 1;
			else SpikeArray[29] = 0;*/
			//SpikeArray[15*N_S_pad+29] = 1;

			/*if(t == 19 ) SpikeArray[32] = 1;
			else SpikeArray[32] = 0;*/
			//SpikeArray[29*N_S_pad+32] = 1;

			/*if(t == 25){
				SpikeArray[41] = 1;
				//SpikeArray[73+N_Group_T] = 1;

			}
			else {
				SpikeArray[41] = 0;
				//SpikeArray[73+N_Group_T] = 0;

			}*/
			//SpikeArray[25*N_S_pad+41] = 1;

			/*if(t == 26 ) SpikeArray[14] = 1;
			else SpikeArray[14] = 0;*/
			//SpikeArray[26*N_S_pad+14] = 1;

			/*if(t == 28 || t == 30 ) SpikeArray[22] = 1;
			else SpikeArray[22] = 0;*/
			//SpikeArray[28*N_S_pad+22] = 1;
			//SpikeArray[30*N_S_pad+22] = 1;

			/*if(t == 28 ) SpikeArray[31] = 1;
			else SpikeArray[31] = 0;*/
			//SpikeArray[28*N_S_pad+31] = 1;

			/*if(t == 29 ) SpikeArray[24] = 1;
			else SpikeArray[24] = 0;*/
			//SpikeArray[29*N_S_pad+24] = 1;
		#endif

		printf("Writing to LMem.\n");
		t0 = getTimestamp();
		Simulation_writeLMem(size_d, 0, x);
		printf("Neurons Write Done.\n");
		fflush(stdout);
		Simulation_writeLMem(syn_size, size_d, syn);
		printf("Synapses Write Done.\n");
		fflush(stdout);
		Simulation_writeLMem(SpikeArray_size, size_d + syn_size, SpikeArray);
		printf("Spikes Write Done.\n");
		fflush(stdout);
		t1 = getTimestamp();
		memory_writes = (t1-t0)/1000000.0;
		printf("Running on DFE.\n");
		fflush(stdout);

		long BurstLengthInBytes = 384;
		//UpdateNeurons(M, N, adex_param_size, stdp_param_size, steps, adex_params, stdp_params);
		//UpdateNeurons(scalar, M, N, M, adex_param_size, stdp_param_size, steps, adex_params, stdp_params, y ,s);
		t0 = getTimestamp();
		//printf("loopLength = %d",UpdateSynapses_post_get_UpdateSynapses_postKernel_loopLength());
		Simulation(N_Group_S_pad, N_Group_T_pad, N_S_pad, adex_param_size, stdp_param_size, timesteps, BurstLengthInBytes, adex_params, stdp_params);
		t1 = getTimestamp();
		DFE_time = (t1-t0)/1000000.0;
		printf("DFE Runtime = %9.7f seconds\n", DFE_time);
		fflush(stdout);

		t0 = getTimestamp();
		Simulation_readLMem(size_d, 0, x);
		Simulation_readLMem(syn_size, size_d, syn);
		t1 = getTimestamp();
		memory_reads = (t1-t0)/1000000.0;

		printf("memory_reads = %lf, memory_writes = %lf, DFE_time = %lf\n", memory_reads, memory_writes, DFE_time);
		fflush(stdout);
		print_neurons(x, N_Group_T);
		print_synapses(syn,N_S, N_S_pad, N_Group_S, N_Group_S_pad, N_Group_T, N_Group_T_pad);
		fflush(stdout);
	}

	fclose(f);
	fclose(g);
	fclose(h);
	fclose(in);
	return 0;
}
