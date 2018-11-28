#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include "Maxfiles.h"
#include "MaxSLiCInterface.h"

//#define allconnected

double defaultclock_dt = 1*1e-3;	//ms

//double taum_S = 10 * 1e-3; //ms
double Ee = 0 * 1e-3; //mV
double taue = 2 * 1e-3; //ms
double Fon = 50; //Hz
double Foff = 3; //Hz

double s = 100 * 1e-10;

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

double AFBn = 0.1771;
double tau_FBn = 0.0327 * 1e3 * 1e-3;	//ms
double AFBp = 0.1548;
double tau_FBp = 0.2302 * 1e3 * 1e-3;	//ms
double AFFp = 0.0618;
double tau_FFp = 0.0666 * 1e3 * 1e-3;	//ms

double etaU = 0.15;
double etaA = 0.15;

//# Adex Parameters
double C = 281*1e-12;	//pF
double gL = 30*1e-9; //nS
double taum = 281*1e-12 / 30*1e-9;	// C/gL
double EL = -70.6*1e-3;	//mV
double DeltaT = 2*1e-3;	//mV
double vti = -50.4*1e-3;	//mV
double vtrest = -45*1e-3;	//mV
double VTmax = 18*1e-3;	//mV
double tauvt = 50*1e-3;	//ms

double tauw = 144*1e-3;	//ms
double c = 4*1e-9;	//ns
double b = 0.0805*1e-9;	//nA
double Vr = -70.6*1e-3;	//mV

int allconnected = 1;

typedef unsigned long long timestamp;
// Time stamp function for measuring DFE runtime
static timestamp getTimestamp(){

	struct timeval now;
	gettimeofday(&now, NULL);
	return now.tv_usec + (timestamp)now.tv_sec * 1000000;

}
// Kaham Sum is used for more accuracy on the addition of double variables errors
double kahanSum(double input, double tosum, double times)
{
	double c=0.0, sum=input,y,t;
	int count;
	for(count=0; count<times; count++)
	{
		y=tosum-c;
		t=sum+y;
		c=(t-sum)-y;
		sum=t;
	}
	return(sum);
}
void print_neurons(float* x, int M){
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
		printf("%.8e\n",x[i]);
	}
	printf("ID\n");
	for(int i = 5; i < M*6; i+=6){
		printf("%.8e\n",x[i]);
	}
}
void print_synapses(float* syn, uint64_t N_S, uint64_t N_S_pad, uint64_t N_Group_S, uint64_t N_Group_S_pad, uint64_t N_Group_T, uint64_t N_Group_T_pad){
	printf("conn\n");
	for(uint64_t i = 0; i < N_Group_S*12; i+=12){
		for(uint64_t j = 0; j < N_Group_T; j++){
			printf("%d, ", (int)syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
		}
		printf("\n");
	}
	for(uint64_t i = N_Group_S_pad*12; i < (N_Group_S_pad+N_S)*12; i+=12){
		for(uint64_t j = 0; j < N_Group_T; j++){
			printf("%f, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
		}
		printf("\n");
	}
	printf("\nw\n");
	for(uint64_t i = 1; i < N_Group_S*12; i+=12){
		for(uint64_t j = 0; j < N_Group_T; j++){
			printf("%.8f, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
		}
		printf("\n");
	}
	for(uint64_t i = N_Group_S_pad*12+1; i < (N_Group_S_pad+N_S)*12; i+=12){
		for(uint64_t j = 0; j < N_Group_T; j++){
			printf("%.8f, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
		}
		printf("\n");
	}
	printf("\nFFp\n");
	for(uint64_t i = 2; i < N_Group_S*12; i+=12){
		for(uint64_t j = 0; j < N_Group_T; j++){
			printf("%.8f, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
		}
		printf("\n");
	}
	for(uint64_t i = N_Group_S_pad*12+2; i < (N_Group_S_pad+N_S)*12; i+=12){
		for(uint64_t j = 0; j < N_Group_T; j++){
			printf("%.8f, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
		}
		printf("\n");
	}
	printf("\nFBp\n");
	for(uint64_t i = 3; i < N_Group_S*12; i+=12){
		for(uint64_t j = 0; j < N_Group_T; j++){
			printf("%.8f, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
		}
		printf("\n");
	}
	for(uint64_t i = N_Group_S_pad*12+3; i < (N_Group_S_pad+N_S)*12; i+=12){
		for(uint64_t j = 0; j < N_Group_T; j++){
			printf("%.8f, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
		}
		printf("\n");
	}
	printf("\nFBn\n");
	for(uint64_t i = 4; i < N_Group_S*12; i+=12){
		for(uint64_t j = 0; j < N_Group_T; j++){
			printf("%.8f, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
		}
		printf("\n");
	}
	for(uint64_t i = N_Group_S_pad*12+4; i < (N_Group_S_pad+N_S)*12; i+=12){
		for(uint64_t j = 0; j < N_Group_T; j++){
			printf("%.8f, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
		}
		printf("\n");
	}
	printf("\nR\n");
	for(uint64_t i = 5; i < N_Group_S*12; i+=12){
		for(uint64_t j = 0; j < N_Group_T; j++){
			printf("%.8g, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
		}
		printf("\n");
	}
	for(uint64_t i = N_Group_S_pad*12+5; i < (N_Group_S_pad+N_S)*12; i+=12){
		for(uint64_t j = 0; j < N_Group_T; j++){
			printf("%.8g, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
		}
		printf("\n");
	}
	printf("\nu\n");
	for(uint64_t i = 6; i < N_Group_S*12; i+=12){
		for(uint64_t j = 0; j < N_Group_T; j++){
			printf("%.8f, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
		}
		printf("\n");
	}
	for(uint64_t i = N_Group_S_pad*12+6; i < (N_Group_S_pad+N_S)*12; i+=12){
		for(uint64_t j = 0; j < N_Group_T; j++){
			printf("%.8f, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
		}
		printf("\n");
	}
	printf("\nU\n");
	for(uint64_t i = 7; i < N_Group_S*12; i+=12){
		for(uint64_t j = 0; j < N_Group_T; j++){
			printf("%.8g, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
		}
		printf("\n");
	}
	for(uint64_t i = N_Group_S_pad*12+7; i < (N_Group_S_pad+N_S)*12; i+=12){
		for(uint64_t j = 0; j < N_Group_T; j++){
			printf("%.8g, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
		}
		printf("\n");
	}
	printf("\nA\n");
	for(uint64_t i = 8; i < N_Group_S*12; i+=12){
		for(uint64_t j = 0; j < N_Group_T; j++){
			printf("%.8f, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
		}
		printf("\n");
	}
	for(uint64_t i = N_Group_S_pad*12+8; i < (N_Group_S_pad+N_S)*12; i+=12){
		for(uint64_t j = 0; j < N_Group_T; j++){
			printf("%.8f, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
		}
		printf("\n");
	}
	printf("\nlastupdate\n");
	for(uint64_t i = 9; i < N_Group_S*12; i+=12){
		for(uint64_t j = 0; j < N_Group_T; j++){
			printf("%f, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
		}
		printf("\n");
	}
	for(uint64_t i = N_Group_S_pad*12+9; i < (N_Group_S_pad+N_S)*12; i+=12){
		for(uint64_t j = 0; j < N_Group_T; j++){
			printf("%f, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
		}
		printf("\n");
	}
	printf("\ntarget_I\n");
	for(uint64_t i = 10; i < N_Group_S*12; i+=12){
		for(uint64_t j = 0; j < N_Group_T; j++){
			printf("%.8f, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
		}
		printf("\n");
	}
	for(uint64_t i = N_Group_S_pad*12+10; i < (N_Group_S_pad+N_S)*12; i+=12){
		for(uint64_t j = 0; j < N_Group_T; j++){
			printf("%.8f, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
		}
		printf("\n");
	}
	printf("\n");
	printf("\nID\n");
	for(uint64_t i = 11; i < N_Group_S*12; i+=12){
		for(uint64_t j = 0; j < N_Group_T; j++){
			printf("%.8f, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
		}
		printf("\n");
	}
	for(uint64_t i = N_Group_S_pad*12+11; i < (N_Group_S_pad+N_S)*12; i+=12){
		for(uint64_t j = 0; j < N_Group_T; j++){
			printf("%.8f, ", syn[j*12*(N_Group_S_pad+N_S_pad)+i]);
		}
		printf("\n");
	}
	printf("\n");
	fflush(stdout);
}
typedef struct {
	int conn;	/**< Variable that expresses if a given synapse between two neurons exists. */
	double w;	/**< Weight of a synapse. (Was present in BRIAN but never used for the ADEX with STDP simulation.) */
	double FFp;	/**< FFp variable of a synapse. (x+) */
	double FBp;	/**< FBp variable of a synapse. (y+) */
	double FBn;	/**< FBn variable of a synapse. (y-) */
	double R;	/**< R value of a synapse. (r)	*/
	double u;	/**< u value of a synapse. (p)	*/
	double U;	/**< U value of a synapse. (P)	*/
	double A;	/**< A value of a synapse. (q)	*/
	double lastupdate;	/**< Last time a synapse was updated	*/
	double target_I;	/**< The I value for the postsynaptic neuron. */
} Synapse;
typedef struct {
    double vt;	/**< Voltage threshold. */
    double vm;	/**< Membrane potential.	*/
    double I;	/**< Neuron input current.	*/
    double x;	/**< Adaption variable (w).	*/
    int Spike;	/**< Variable to show if the neuron has spiked in a specific moment.	*/
} Neuron;
void resetNeuron(Neuron* neuron) {
    neuron->vm = Vr;
    neuron->x += b;
    neuron->vt = VTmax;
}
void SolveNeurons(Neuron* neurons, int N, int *SpikeArray){
    for(int i = 0; i < N; i++){
        double _vm, _vt, _x;
        _vm = (gL*(EL-neurons[i].vm)+gL*DeltaT*exp((neurons[i].vm-neurons[i].vt)/DeltaT)+neurons[i].I-neurons[i].x)/C;
        _vt = -(neurons[i].vt-vtrest)/tauvt;
        _x = (c*(neurons[i].vm-EL)-neurons[i].x)/tauw;
        neurons[i].vm += _vm * defaultclock_dt;
        neurons[i].vt += _vt * defaultclock_dt;
        neurons[i].x += _x * defaultclock_dt;
        if(neurons[i].vm > neurons[i].vt){
            resetNeuron(&neurons[i]);
            SpikeArray[i] = 1;
            neurons[i].Spike  = 1;
        }
        else{
        	SpikeArray[i] = 0;
        	neurons[i].Spike  = 0;
        }
    }
}
void UpdateSynapses_pre(Synapse** Synapses, Neuron* neurons, int N_S, int N_Group_S, int N_Group_T, int* SpikeArray, double t){
	if(N_S == 0) s = (double)500000;
	for (int i = 0; i < N_S + N_Group_S; i++){
		if (SpikeArray[i+N_Group_T-N_Group_S] > 0){
			for (int j = 0; j < N_Group_T; j++){
				if (Synapses[i][j].conn){
					Synapses[i][j].FFp = Synapses[i][j].FFp * exp(-(-Synapses[i][j].lastupdate + t)/tau_FFp);
					Synapses[i][j].FBn = Synapses[i][j].FBn * exp(-(-Synapses[i][j].lastupdate + t)/tau_FBn);
					Synapses[i][j].u = Synapses[i][j].U + (-Synapses[i][j].U + Synapses[i][j].u) * exp(-(-Synapses[i][j].lastupdate + t)/tau_u);
					Synapses[i][j].FBp = Synapses[i][j].FBp * exp(-(-Synapses[i][j].lastupdate + t)/tau_FBp);
					Synapses[i][j].R = (Synapses[i][j].R - 1) * exp(-(-Synapses[i][j].lastupdate + t)/tau_r) + 1;
					Synapses[i][j].target_I = s * Synapses[i][j].A * Synapses[i][j].R * Synapses[i][j].u;
					Synapses[i][j].U = Synapses[i][j].U + etaU * (-AFBn * Synapses[i][j].FBn * Synapses[i][j].FBp + AFBp * Synapses[i][j].FBp * Synapses[i][j].FFp);
					if (Synapses[i][j].U < Umin) Synapses[i][j].U = Umin;
					else if (Synapses[i][j].U > Umax) Synapses[i][j].U = Umax;
					Synapses[i][j].w = Synapses[i][j].U * Synapses[i][j].A;
					Synapses[i][j].FFp += 1;
					Synapses[i][j].R -= Synapses[i][j].R * Synapses[i][j].u;
					Synapses[i][j].u += Synapses[i][j].U * (1 - Synapses[i][j].u);
					Synapses[i][j].lastupdate = t;
				}
			}
		}
	}
	for (int j = 0; j < N_Group_T; j++){
		for (int i = 0; i < N_S + N_Group_S; i++){
			if (Synapses[i][j].conn && SpikeArray[i+N_Group_T-N_Group_S]){
				neurons[j].I = Synapses[i][j].target_I;
			}
		}
	}
}

void UpdateSynapses_post(Synapse** Synapses, int N_S, int N_Group_S, int N_Group_T, int* SpikeArray, double t){
	for (int i = 0; i < N_Group_T; i++){
	    if (SpikeArray[i] > 0){
	        for (int j = 0; j < N_S + N_Group_S; j++){
	            if (Synapses[j][i].conn){
	                Synapses[j][i].FFp = Synapses[j][i].FFp * exp(-(-Synapses[j][i].lastupdate + t)/tau_FFp);
	                Synapses[j][i].FBn = Synapses[j][i].FBn * exp(-(-Synapses[j][i].lastupdate + t)/tau_FBn);
	                Synapses[j][i].u = Synapses[j][i].U + (-Synapses[j][i].U + Synapses[j][i].u) * exp(-(-Synapses[j][i].lastupdate + t)/tau_u);
	                Synapses[j][i].FBp = Synapses[j][i].FBp * exp(-(-Synapses[j][i].lastupdate + t)/tau_FBp);
	                Synapses[j][i].R = (Synapses[j][i].R - 1) * exp(-(-Synapses[j][i].lastupdate + t)/tau_r) + 1;
	            }
	        }
	    }
	}
	double mean = 0;
	int num = 0;
	for (int l=0; l<N_Group_T; l++){
		for (int k=0; k<N_S+N_Group_S; k++){
			if (Synapses[k][l].conn && SpikeArray[l]){
				mean += AFFp * Synapses[k][l].FFp * Synapses[k][l].FBn;
				num++;
			}
		}
	}
	mean = mean / (double)num;
	for (int i = 0; i < N_Group_T; i++){
	    if (SpikeArray[i] > 0){
	        for (int j = 0; j < N_S + N_Group_S; j++){
	            if (Synapses[j][i].conn){
	            	Synapses[j][i].A = Synapses[j][i].A + etaA * (AFFp * Synapses[j][i].FFp * Synapses[j][i].FBn);
	                Synapses[j][i].A = Synapses[j][i].A - etaA * 0.5 * mean;
	                if (Synapses[j][i].A < Amin) Synapses[j][i].A = Amin;
	                else if (Synapses[j][i].A > Amax) Synapses[j][i].A = Amax;
	            }
	        }
	    }
	}
	for (int i = 0; i < N_Group_T; i++){
	    if (SpikeArray[i] > 0){
	        for (int j = 0; j < N_S + N_Group_S; j++){
	            if (Synapses[j][i].conn){
	            	Synapses[j][i].w = Synapses[j][i].U * Synapses[j][i].A;
	                Synapses[j][i].FBp += 1;
	                Synapses[j][i].FBn += 1;
	                Synapses[j][i].lastupdate = t;
	            }
	        }
	    }
	}


}
void print_neurons_CPU(Neuron* neurons, int N){
	printf("vt\n");
	for(int i = 0; i < N; i++){
		printf("%.8e\n",neurons[i].vt);
	}
	printf("\nvm\n");
	for(int i = 0; i < N; i++){
		printf("%.8e\n",neurons[i].vm);
	}
	printf("\nI\n");
	for(int i = 0; i < N; i++){
		printf("%.8e\n",neurons[i].I);
	}
	printf("\nx\n");
	for(int i = 0; i < N; i++){
		printf("%.8e\n",neurons[i].x);
	}
	printf("\nSpike\n");
	for(int i = 0; i < N; i++){
		printf("%d\n",neurons[i].Spike);
	}
	printf("\n");
	fflush(stdout);
}
void print_synapses_CPU(Synapse** syn, int N_S, int N_T){
	printf("conn\n");
	fflush(stdout);
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++){
			printf("%d ", syn[i][j].conn);
		}
		printf("\n");
	}
	printf("\nw\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++){
			printf("%.8e, ", syn[i][j].w);
		}
		printf("\n");
	}
	printf("\nFFp\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++)
			printf("%.8e, ", syn[i][j].FFp);
		printf("\n");
	}
	printf("\nFBp\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++)
			printf("%.8e, ", syn[i][j].FBp);
		printf("\n");
	}
	printf("\nFBn\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++)
			printf("%.8e, ", syn[i][j].FBn);
		printf("\n");
	}
	printf("\nR\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++)
			printf("%.8e, ", syn[i][j].R);
		printf("\n");
	}
	printf("\nu\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++){
			printf("%.8e, ", syn[i][j].u);
		}
		printf("\n");
	}
	printf("\nU\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++){
			printf("%.8e, ", syn[i][j].U);
		}
		printf("\n");
	}
	printf("\nA\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++){
			printf("%.8e, ", syn[i][j].A);
		}
		printf("\n");
	}
	printf("\nlastupdate\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++)
			printf("%lf, ", syn[i][j].lastupdate);
		printf("\n");
	}
	printf("\ntarget_I\n");
	for(int i =0; i < N_S; i++){
		for(int j = 0; j < N_T; j++){
			printf("%.8e, ", syn[i][j].target_I);
		}
		printf("\n");
	}
	printf("\n");
}
void CPUCode(Neuron *neurons, Synapse **syn, int *SpikeArray, int N_S, int N_Group_S, int N_Group_T, int steps, int spike_interval){
	// Initialization of AdEx Neurons
	for(int i = 0; i<N_Group_T; i++){
		neurons[i].vt = vtrest;
		neurons[i].vm = EL;
		if(N_S == 0) neurons[i].vm = vtrest + 0.005;
		neurons[i].I = 0;
		neurons[i].x = 0;
		neurons[i].Spike = 0;
	}
	// Initialization of Synapses
	if(!allconnected){
		FILE *c = fopen("connections.txt", "r");
		int con = 0;
		for(int i = 0; i < N_Group_S+N_S; i++)
			for(int j = 0; j < N_Group_T; j++){
				fscanf (c, "%d", &con);
				printf("%d ", con);
				syn[j][i].conn = con;//(i+j)%2;//1;
			}
		fclose(c);
	}
	else{
		for(int i = 0; i < N_Group_S+N_S; i++)
				for(int j = 0; j < N_Group_T; j++)
					syn[i][j].conn = 1;//(i+j)%2;//1;
	}
	// Initialization of Synapses for Neurons
	double init_const = 0;
	for(int i = 0; i < N_Group_S; i++){
		for(int j = 0; j < N_Group_T; j++){
			if (syn[i][j].conn) {
			syn[i][j].w = 0;
			syn[i][j].FFp = 0;
				syn[i][j].FBp = init_const;//0;
				syn[i][j].FBn = init_const;//0;
				syn[i][j].R = init_const;//1;
			syn[i][j].u = 0;
			if(N_S == 0) syn[i][j].u = 1;
				syn[i][j].U = init_const;//exp(-(((pow((init_const+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(Umax-Umin)+Umin;//init_const;//exp(-(((pow((init_const+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(Umax-Umin)+Umin;	// takes time
				syn[i][j].A = init_const;//exp(-(((pow((init_const+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(Amax-Amin)+Amin;//init_const;//exp(-(((pow((init_const+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(Amax-Amin)+Amin;	// takes time
				syn[i][j].lastupdate = 0;
			syn[i][j].target_I = 0;
			init_const = (double)((int)(init_const+1)%16777216);
			}
		}
	}
	// Initialization of Synapses for external input (bottom rows)
	for(int i = N_Group_S; i < N_Group_S+N_S; i++){
		for(int j = 0; j < N_Group_T; j++){
			if (syn[i][j].conn) {
			syn[i][j].w = 0;
			syn[i][j].FFp = 0;
				syn[i][j].FBp = init_const;//0;
				syn[i][j].FBn = init_const;//0;
				syn[i][j].R = init_const;//1;
			syn[i][j].u = 0;
				syn[i][j].U = init_const;//exp(-(((pow(((init_const)+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(Umax-Umin)+Umin;	// takes time
				syn[i][j].A = init_const;//exp(-(((pow(((init_const)+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(Amax-Amin)+Amin;	// takes time
			syn[i][j].lastupdate = 0;
			syn[i][j].target_I = 0;
			init_const = (double)((int)(init_const+1)%16777216);
			}
		}
	 }
	for(int t = 0; t < steps; t++){

		SolveNeurons(neurons, N_Group_T, SpikeArray);

		for(int i = N_Group_T; i < N_Group_T+N_S; i++) SpikeArray[i]= (int)(t%spike_interval==0);	// Input Spikes Initialization

		UpdateSynapses_pre(syn, neurons, N_S, N_Group_S, N_Group_T, SpikeArray, (double)(t*defaultclock_dt));

		UpdateSynapses_post(syn, N_S, N_Group_S, N_Group_T, SpikeArray, (double)(t*defaultclock_dt));

	}
}
int main(int argc, char * argv[])
{
	timestamp t0,t1,t2,t3;
	max_file_t *maxfile = Simulation_init();
	int BurstLengthBytes = max_get_burst_size(maxfile, NULL);
	int Unroll_Factor = 2;

	int nruns = 1;

	double realtime = 0;
	double stime = 0.050; //second
	double stime2 = 50; //second

	double resolution_export = 10 * 1e-3; //every x ms

	// Reading parameters of the simulation from the terminal
	int steps = atoi(argv[1]);						// Number of Simulation Steps
	uint64_t N_S = (uint64_t)atoi(argv[2]);			// Number of Input Neurons
	uint64_t N_Group_S = (uint64_t)atoi(argv[3]);	// Number of AdEx Neurons as Source
	uint64_t N_Group_T = (uint64_t)atoi(argv[4]);	// Number of AdEx Neurons as Target
	int SpikeStep = atoi(argv[5]);					// Input Neurons Spike Interval
	uint32_t WriteEn = (uint32_t)atoi(argv[6]);		// Write Enable to copy AdEx Neuron values from DFE to CPU at every timestep
	uint32_t WriteStep = (uint32_t)1;
	allconnected = atoi(argv[7]);					// 0: Connections read from "connections.txt". 1: All connected
	int RunCPU = atoi(argv[8]);						// Run CPU Code with doubles to check acceleration and errors
	int MxM;
	if (N_S == 0) MxM = 1;
	else MxM = 0;
	if (MxM) s = 500000; // For testing

	// Creating padding for the Neurons, which need to be BurstSize alligned
	uint64_t N_S_pad = N_S, N_Group_S_pad = N_Group_S, N_Group_T_pad = N_Group_T;
	if(N_S%BurstLengthBytes != 0) N_S_pad = (N_S/BurstLengthBytes)*BurstLengthBytes+BurstLengthBytes;
	if(N_Group_S%BurstLengthBytes != 0) N_Group_S_pad = (N_Group_S/BurstLengthBytes)*BurstLengthBytes+BurstLengthBytes;
	if(N_Group_T%BurstLengthBytes != 0) N_Group_T_pad = (N_Group_T/BurstLengthBytes)*BurstLengthBytes+BurstLengthBytes;

	printf("N_S_pad = %d, N_Group_S_pad = %d, N_Group_T_pad = %d\n", (int)N_S_pad, (int)N_Group_S_pad, (int)N_Group_T_pad);
	fflush(stdout);
	double input1_pos = 25;
	double input2_pos = 75;
	double rad = 5;

	//Define input 1
	double *F_input1;
	F_input1 = (double*)malloc(sizeof(double)*N_S);
	for(int i = 0; i < N_S; i++){
		F_input1[i] = Foff;
		F_input1[i] = exp(-(pow((i+1)-input1_pos,2)/(2.0*pow(rad,2))))*(Fon-Foff)+Foff; //Define gaussian input
	}

	//Define input 2
	double *F_input2;
	F_input2 = (double*)malloc(sizeof(double)*N_S);
	for(int i = 0; i < N_S; i++){
		F_input2[i] = Foff;
		F_input2[i] = exp(-(pow((i+1)-input2_pos,2)/(2.0*pow(rad,2))))*(Fon-Foff)+Foff; //Define gaussian input
	}
	printf("F_input done.\n");
	fflush(stdout);
	// Creating random Input Arrays
	/*Poisson *input;
	input = (Poisson*)malloc(sizeof(Poisson)*N_S);
	for(int i = 0; i<N_S; i++){				// Initialization of Poisson Neurons
		input[i].GaussArray = F_input1[i];
		input[i].Spike = 0;
	}*/

	// ******** ADEX DE Parameters Initialization *********** //
	int adex_param_size = 16;	// word alligned size of adex parameters array
	int adex_param_array_size = adex_param_size*sizeof(float);	// adex params array bytes
	float *adex_params = malloc(adex_param_array_size);
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
	int stdp_param_array_size = stdp_param_size*sizeof(float);
	float *stdp_params = malloc(stdp_param_array_size);
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
	uint64_t size_d = N_Group_T_pad * 6;
	uint64_t sizeBytes_d_var = size_d * (uint64_t)sizeof(float);
	float *x = malloc(sizeBytes_d_var);
	for(uint64_t i = 0; i<N_Group_T_pad*6; i+=6) {
		x[i] = vtrest;
		if(MxM)	x[i+1] = vtrest + 0.005;
		else x[i+1] = EL;
		x[i+2] = 0;
		x[i+3] = 0;
		x[i+4] = 0;	// Spike
		x[i+5] = i;
	}
	printf("Neuron Initialization done.\n");
	fflush(stdout);
	// ******** Synapses Initialization ******** //
	uint64_t syn_size = (N_S_pad+N_Group_S_pad) * N_Group_T_pad * 12;
	uint64_t syn_sizeBytes = syn_size * (uint64_t)sizeof(float);
	float* syn = malloc(syn_sizeBytes);
	if(!allconnected){
		FILE *in = fopen("connections.txt", "r");
		int con = 0;
		for(uint64_t j = 0; j < N_Group_T_pad; j++)
			for(uint64_t i = 0; i < N_Group_S_pad*12; i+=12)
				if((i < N_Group_S*12) & (j < N_Group_T)){
					fscanf (in, "%d", &con);
					syn[i+j*(N_Group_S_pad+N_S_pad)*12] = con;
				}
				else syn[i+j*(N_Group_S_pad+N_S_pad)*12] = 0;
		for(uint64_t j = 0; j < N_Group_T_pad; j++)
			for(uint64_t i = N_Group_S_pad*12; i < (N_Group_S_pad+N_S_pad)*12; i+=12)
				if((i-N_Group_S_pad*12 < N_S*12) & (j < N_Group_T)){
					fscanf (in, "%d", &con);
					syn[i+j*(N_Group_S_pad+N_S_pad)*12] = con;
				}
				else syn[i+j*(N_Group_S_pad+N_S_pad)*12] = 0;
		fclose(in);
	}
	else{
		for(uint64_t j = 0; j < N_Group_T_pad; j++)
			for(uint64_t i = 0; i < N_Group_S_pad*12; i+=12)
				if((i < N_Group_S*12) & (j < N_Group_T)) syn[i+j*(N_Group_S_pad+N_S_pad)*12] = 1;//(i/12+j)%2;//1;
				else syn[i+j*(N_Group_S_pad+N_S_pad)*12] = 0;
		for(uint64_t j = 0; j < N_Group_T_pad; j++)
			for(uint64_t i = N_Group_S_pad*12; i < (N_Group_S_pad+N_S_pad)*12; i+=12)
				if((i-N_Group_S_pad*12 < N_S*12) & (j < N_Group_T)) syn[i+j*(N_Group_S_pad+N_S_pad)*12] = 1;//(i/12+j)%2;//1;
				else syn[i+j*(N_Group_S_pad+N_S_pad)*12] = 0;
	}
	printf("Synapses Connectivity Initialization done.\n");
	fflush(stdout);
	float init_const = (float)0;
	for(uint64_t i = 0; i < N_Group_S_pad*12; i+=12){
		for(uint64_t j = 0; j < N_Group_T_pad; j++){
			if (syn[i+j*(N_Group_S_pad+N_S_pad)*12]) {
				syn[i+j*(N_Group_S_pad+N_S_pad)*12+3] = init_const;//0;
				syn[i+j*(N_Group_S_pad+N_S_pad)*12+4] = init_const;//0;
				syn[i+j*(N_Group_S_pad+N_S_pad)*12+5] = init_const;//1;
				if (MxM) syn[i+j*(N_Group_S_pad+N_S_pad)*12+6] = 1;	// for testing
				syn[i+j*(N_Group_S_pad+N_S_pad)*12+7] = init_const;//exp(-(((pow((init_const+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(stdp_params[8]-stdp_params[9])+stdp_params[9];	// takes time
				syn[i+j*(N_Group_S_pad+N_S_pad)*12+8] = init_const;//exp(-(((pow((init_const+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(stdp_params[5]-stdp_params[6])+stdp_params[6];	// takes time
				syn[i+j*(N_Group_S_pad+N_S_pad)*12+11] = init_const;
				init_const = (float)((int)(init_const+1)%16777216); //max float number 16777216
			}
		}
	}
	// Initialization of Synapses for external input (bottom rows)
	for(uint64_t i = N_Group_S_pad*12; i < (N_Group_S_pad+N_S_pad)*12; i+=12){
		for(uint64_t j = 0; j < N_Group_T_pad; j++){
			if (syn[i+j*(N_Group_S_pad+N_S_pad)*12]) {
				syn[i+j*(N_Group_S_pad+N_S_pad)*12+3] = init_const;//0;//init_const;//0;
				syn[i+j*(N_Group_S_pad+N_S_pad)*12+4] = init_const;//0;
				syn[i+j*(N_Group_S_pad+N_S_pad)*12+5] = init_const;//1;
				syn[i+j*(N_Group_S_pad+N_S_pad)*12+7] = init_const;//exp(-(((pow(((init_const)+1)-input1_pos,2)))/(2.0*pow(rad+0,2))))*(stdp_params[8]-stdp_params[9])+stdp_params[9];	// takes time
				syn[i+j*(N_Group_S_pad+N_S_pad)*12+8] = init_const;//exp(-(((pow(((init_const)+1)-input1_pos,2)))/(2.0*pow(rad+3,2))))*(stdp_params[5]-stdp_params[6])+stdp_params[6];	// takes time
				syn[i+j*(N_Group_S_pad+N_S_pad)*12+11] = init_const;
				init_const = (float)((int)(init_const+1)%16777216); //max float number 16777216
			}
		}
	}
	printf("Synapses Initialization done.\n");
	fflush(stdout);
	uint32_t timesteps = steps;
	printf("timesteps=%d\n",timesteps);
	fflush(stdout);

	// ********** Input Neurons Spikes Initialization ********** //
	uint64_t SpikeArray_size = N_S_pad/Unroll_Factor * (uint64_t)timesteps;
	uint64_t SpikeArray_bytes = SpikeArray_size * (uint64_t)sizeof(uint8_t);
	uint8_t* SpikeArray = malloc(SpikeArray_bytes);
	uint8_t k;
	for(uint64_t i = 0; i < timesteps; i++){
		for(uint64_t j = 0; j < N_S_pad/Unroll_Factor; j++){
			k = (uint8_t)((i%SpikeStep)==0);
			SpikeArray[i*(N_S_pad/Unroll_Factor)+j] = k | k << 4;
		}
	}
	printf("Input Spikes Initialization done.\n");
	fflush(stdout);
	srand(time(NULL));

	float* NeuronCPU = malloc(4*N_Group_T_pad*(timesteps/WriteStep)*sizeof(float));

	double init_time = 0, load_time = 0, unload_time = 0, memory_reads = 0, memory_writes = 0, DFE_time = 0;

	printf("Writing to LMem.\n");
	fflush(stdout);
	t0 = getTimestamp();
	Simulation_writeLMem(size_d, 0, x);				// Writing AdEx Neurons values to LMem
	printf("Neurons Write Done.\n");
	fflush(stdout);
	Simulation_writeLMem(syn_size, size_d, syn);	// Writing Synapses values to LMem
	printf("Synapses Write Done.\n");
	fflush(stdout);
	t1 = getTimestamp();
	memory_writes = (t1-t0)/1000000.0;
	printf("Running on DFE.\n");
	fflush(stdout);

	// Dividing DFE Simulation into window calculations of 1 hour
	uint64_t max_ticks = 249205632032 + 1; // 53minutes in 4992x4992, MxM, 500 timesteps experiment
	uint64_t ticks = 32 + (2*N_Group_T_pad + 16*(N_S_pad+N_Group_S_pad)*N_Group_T_pad + 4*(N_S_pad+N_Group_S_pad)*N_Group_T_pad)*timesteps;
	int rec = (int)(ticks/max_ticks);
	int step = 0;
	if(rec!=0){
		if((timesteps/rec)%2==0){
			step = timesteps/(rec+1);
		}
		else{
			step = timesteps/rec - 1;
			rec = timesteps/step;
		}
	}

	// Calling Kernels
	long BurstLengthInBytes = BurstLengthBytes;
	t0 = getTimestamp();
	uint64_t spike_start_addr = size_d + syn_size;
	for(int i = 0; i < rec; i++){
		if(!MxM){
			t2 = getTimestamp();
			Simulation_writeLMemInt(N_S_pad/Unroll_Factor * (uint64_t)step, spike_start_addr, SpikeArray);	// Writing Input Neurons Spikes to LMem
			t3 = getTimestamp();
			memory_writes += (t3-t2)/1000000.0;
			SpikeArray += N_S_pad/Unroll_Factor * (uint64_t)step;
		}
		Simulation(N_Group_S_pad, N_Group_T_pad, N_S_pad, adex_param_size, stdp_param_size, (float)(i*step*defaultclock_dt), step, WriteEn, WriteStep, BurstLengthInBytes, adex_params, stdp_params, NeuronCPU);
		NeuronCPU += 4 * N_Group_T_pad * ((i+1)*step/WriteStep-i*step/WriteStep);
	}
	if(rec*step<timesteps){
		int old_rec_step = rec*step;
		step = timesteps-rec*step;
		if(!MxM){
			t2 = getTimestamp();
			Simulation_writeLMemInt(N_S_pad/Unroll_Factor * (uint64_t)step, spike_start_addr, SpikeArray);	// Writing Input Neurons Spikes to LMem
			t3 = getTimestamp();
			memory_writes += (t3-t2)/1000000.0;
			SpikeArray += N_S_pad/Unroll_Factor * (uint64_t)step;
		}
		Simulation(N_Group_S_pad, N_Group_T_pad, N_S_pad, adex_param_size, stdp_param_size, (float)(old_rec_step*defaultclock_dt), step, WriteEn, WriteStep, BurstLengthInBytes, adex_params, stdp_params, NeuronCPU);
	}

	t1 = getTimestamp();
	DFE_time = (t1-t0)/1000000.0;
	printf("DFE Runtime = %9.7f seconds\n", DFE_time);
	fflush(stdout);

	t0 = getTimestamp();
	Simulation_readLMem(size_d, 0, x);				// Reading AdEx Neurons Values
	Simulation_readLMem(syn_size, size_d, syn);		// Reading Synapses Values
	t1 = getTimestamp();
	memory_reads = (t1-t0)/1000000.0;

	printf("memory_reads = %lf, memory_writes = %lf, DFE_time = %lf\n", memory_reads, memory_writes, DFE_time);
	fflush(stdout);

	// Printing Final Values of AdEx Neurons and Synapses
	print_neurons(x, N_Group_T);
	print_synapses(syn,N_S, N_S_pad, N_Group_S, N_Group_S_pad, N_Group_T, N_Group_T_pad);
	//fflush(stdout);

	// Printing to file AdEx Neurons Variables at every timestep
	if(WriteEn){
		FILE *f = fopen("NeuronOut.csv", "w");
		fprintf(f, "Time,ID,vm,vt,x,Spike\n");
		for(int j = 0; j < timesteps; j += WriteStep){
			for(int i = 0; i < N_Group_T*4; i+=4){
				fprintf(f, "%.8e,%d,%.8e,%.8e,%.8e,%.8e\n", j*defaultclock_dt,i/4,NeuronCPU[i],NeuronCPU[i+1],NeuronCPU[i+2],NeuronCPU[i+3]);
			}
		}
		fclose(f);
	}

	// ***************** CPU CODE *************** //
	if(RunCPU){

		// Memory allocation for AdEx Neurons, Synapses and Input Neurons
		Neuron *neurons;
		neurons = (Neuron*)calloc(((int)N_Group_T),sizeof(Neuron));
		printf("Neurons malloc done\n");
		fflush(stdout);
		int *SpikeArray_CPU;
		SpikeArray_CPU = (int*)calloc((int)(N_Group_T+N_S),sizeof(int));
		printf("Spike Array malloc done\n");
		fflush(stdout);
		Synapse *syn_CPU[N_S+N_Group_S];
		for(uint64_t i = 0; i < N_S+N_Group_S; i++) syn_CPU[i] = (Synapse*)calloc((N_Group_T),sizeof(Synapse));
		printf("Synapses malloc done\n");
		fflush(stdout);

		// Calling function to run the CPU Simulation
		t0 = getTimestamp();
		CPUCode(neurons, syn_CPU, SpikeArray_CPU, (int)N_S, (int)N_Group_S, (int)N_Group_T, (int)timesteps, (int)SpikeStep);
		t1 = getTimestamp();

		printf("CPU Time: %lf\n", (t1-t0)/1000000.0);
		printf("Acceleration: %lf\n", (t1-t0)/(1000000.0*DFE_time));
		fflush(stdout);

		// Printing final values of CPU Simulation
		//print_neurons_CPU(neurons, N_Group_T);
		//print_synapses_CPU(syn_CPU,N_S+N_Group_S,N_Group_T);

		// Calculating Differences of DFE and CPU Simulations
		double sum_CPU = 0.0, sum_DFE = 0.0, error = 0.0;
		double average_error = 0.0, relative_average_error = 0.0;
		for(uint64_t i = 0; i < N_Group_T; i++){
			kahanSum(sum_CPU,neurons[i].vt,1);
			kahanSum(sum_DFE,(double)x[i*6],1);
		}
		average_error = kahanSum(sum_CPU,-sum_DFE,1)/N_Group_T;
		relative_average_error = average_error/(sum_CPU/N_Group_T);
		printf("Neurons vt: Average error= %.8e Relative Average Error= %.8e\n", average_error, relative_average_error);
		sum_CPU = 0.0;
		sum_DFE = 0.0;
		error = 0.0;
		for(uint64_t i = 0; i < N_Group_T; i++){
			sum_CPU = kahanSum(sum_CPU,neurons[i].vm,1);
			sum_DFE = kahanSum(sum_DFE,(double)x[i*6+1],1);
		}
		average_error = kahanSum(sum_CPU,-sum_DFE,1)/N_Group_T;
		relative_average_error = average_error/(sum_CPU/N_Group_T);
		printf("Neurons vm: Average error= %.8e Relative Average Error= %.8e\n", average_error, relative_average_error);
		sum_CPU = 0.0;
		sum_DFE = 0.0;
		error = 0.0;
		for(uint64_t i = 0; i < N_Group_T; i++){
			sum_CPU = kahanSum(sum_CPU,neurons[i].I,1);
			sum_DFE = kahanSum(sum_DFE,(double)x[i*6+2],1);
		}
		average_error = kahanSum(sum_CPU,-sum_DFE,1)/N_Group_T;
		relative_average_error = average_error/(sum_CPU/N_Group_T);
		printf("Neurons I: Average error= %.8e Relative Average Error= %.8e\n", average_error, relative_average_error);
		sum_CPU = 0.0;
		sum_DFE = 0.0;
		error = 0.0;
		for(uint64_t i = 0; i < N_Group_T; i++){
			sum_CPU = kahanSum(sum_CPU,neurons[i].x,1);
			sum_DFE = kahanSum(sum_DFE,(double)x[i*6+3],1);
		}
		average_error = kahanSum(sum_CPU,-sum_DFE,1)/N_Group_T;
		relative_average_error = average_error/(sum_CPU/N_Group_T);
		printf("Neurons x: Average error= %.8e Relative Average Error= %.8e\n", average_error, relative_average_error);

		sum_CPU = 0.0;
		sum_DFE = 0.0;
		error = 0.0;
		for(uint64_t i = 0; i < N_Group_T; i++){
			sum_CPU = kahanSum(sum_CPU,neurons[i].Spike,1);
			sum_DFE = kahanSum(sum_DFE,(double)x[i*6+4],1);
		}
		average_error = kahanSum(sum_CPU,-sum_DFE,1)/N_Group_T;
		relative_average_error = average_error/(sum_CPU/N_Group_T);
		printf("Neurons Spikes: Average error= %.8e Relative Average Error= %.8e\n", average_error, relative_average_error);

		sum_CPU = 0.0;
		sum_DFE = 0.0;
		error = 0.0;
		for(uint64_t i = 0; i < N_Group_T; i++){
			for(uint64_t j = 0; j < N_Group_S+N_S; j++){
				sum_CPU = kahanSum(sum_CPU,syn_CPU[j][i].w,1);
			}
		}
		for(uint64_t i = 1; i < N_Group_S*12; i+=12){
			for(uint64_t j = 0; j < N_Group_T; j++){
				sum_DFE = kahanSum(sum_DFE,(double)syn[j*12*(N_Group_S_pad+N_S_pad)+i],1);
			}
		}
		for(uint64_t i = N_Group_S_pad*12+1; i < (N_Group_S_pad+N_S)*12; i+=12){
			for(uint64_t j = 0; j < N_Group_T; j++){
				sum_DFE = kahanSum(sum_DFE,(double)syn[j*12*(N_Group_S_pad+N_S_pad)+i],1);
			}
		}
		average_error = kahanSum(sum_CPU,-sum_DFE,1)/((N_S+N_Group_S)*N_Group_T);
		relative_average_error = average_error/(sum_CPU/((N_S+N_Group_S)*N_Group_T));
		printf("Synapses w: Average error= %.8e Relative Average Error= %.8e\n", average_error, relative_average_error);

		sum_CPU = 0.0;
		sum_DFE = 0.0;
		error = 0.0;
		for(uint64_t i = 0; i < N_Group_T; i++){
			for(uint64_t j = 0; j < N_Group_S+N_S; j++){
				sum_CPU = kahanSum(sum_CPU,syn_CPU[j][i].FFp,1);
			}
		}
		for(uint64_t i = 2; i < N_Group_S*12; i+=12){
			for(uint64_t j = 0; j < N_Group_T; j++){
				sum_DFE = kahanSum(sum_DFE,(double)syn[j*12*(N_Group_S_pad+N_S_pad)+i],1);
			}
		}
		for(uint64_t i = N_Group_S_pad*12+2; i < (N_Group_S_pad+N_S)*12; i+=12){
			for(uint64_t j = 0; j < N_Group_T; j++){
				sum_DFE = kahanSum(sum_DFE,(double)syn[j*12*(N_Group_S_pad+N_S_pad)+i],1);
			}
		}
		average_error = kahanSum(sum_CPU,-sum_DFE,1)/((N_S+N_Group_S)*N_Group_T);
		relative_average_error = average_error/(sum_CPU/((N_S+N_Group_S)*N_Group_T));
		printf("Synapses FFp: Average error= %.8e Relative Average Error= %.8e\n", average_error, relative_average_error);

		sum_CPU = 0.0;
		sum_DFE = 0.0;
		error = 0.0;
		for(uint64_t i = 0; i < N_Group_T; i++){
			for(uint64_t j = 0; j < N_Group_S+N_S; j++){
				sum_CPU = kahanSum(sum_CPU,syn_CPU[j][i].FBp,1);
			}
		}
		for(uint64_t i = 3; i < N_Group_S*12; i+=12){
			for(uint64_t j = 0; j < N_Group_T; j++){
				sum_DFE = kahanSum(sum_DFE,(double)syn[j*12*(N_Group_S_pad+N_S_pad)+i],1);
			}
		}
		for(uint64_t i = N_Group_S_pad*12+3; i < (N_Group_S_pad+N_S)*12; i+=12){
			for(uint64_t j = 0; j < N_Group_T; j++){
				sum_DFE = kahanSum(sum_DFE,(double)syn[j*12*(N_Group_S_pad+N_S_pad)+i],1);
			}
		}
		average_error = kahanSum(sum_CPU,-sum_DFE,1)/((N_S+N_Group_S)*N_Group_T);
		relative_average_error = average_error/(sum_CPU/((N_S+N_Group_S)*N_Group_T));
		printf("Synapses FBp: Average error= %.8e Relative Average Error= %.8e\n", average_error, relative_average_error);

		sum_CPU = 0.0;
		sum_DFE = 0.0;
		error = 0.0;
		for(uint64_t i = 0; i < N_Group_T; i++){
			for(uint64_t j = 0; j < N_Group_S+N_S; j++){
				sum_CPU = kahanSum(sum_CPU,syn_CPU[j][i].FBn,1);
			}
		}
		for(uint64_t i = 4; i < N_Group_S*12; i+=12){
			for(uint64_t j = 0; j < N_Group_T; j++){
				sum_DFE = kahanSum(sum_DFE,(double)syn[j*12*(N_Group_S_pad+N_S_pad)+i],1);
			}
		}
		for(uint64_t i = N_Group_S_pad*12+4; i < (N_Group_S_pad+N_S)*12; i+=12){
			for(uint64_t j = 0; j < N_Group_T; j++){
				sum_DFE = kahanSum(sum_DFE,(double)syn[j*12*(N_Group_S_pad+N_S_pad)+i],1);
			}
		}
		average_error = kahanSum(sum_CPU,-sum_DFE,1)/((N_S+N_Group_S)*N_Group_T);
		relative_average_error = average_error/(sum_CPU/((N_S+N_Group_S)*N_Group_T));
		printf("Synapses FBn: Average error= %.8e Relative Average Error= %.8e\n", average_error, relative_average_error);

		sum_CPU = 0.0;
		sum_DFE = 0.0;
		error = 0.0;
		for(uint64_t i = 0; i < N_Group_T; i++){
			for(uint64_t j = 0; j < N_Group_S+N_S; j++){
				sum_CPU = kahanSum(sum_CPU,syn_CPU[j][i].R,1);
			}
		}
		for(uint64_t i = 5; i < N_Group_S*12; i+=12){
			for(uint64_t j = 0; j < N_Group_T; j++){
				sum_DFE = kahanSum(sum_DFE,(double)syn[j*12*(N_Group_S_pad+N_S_pad)+i],1);
			}
		}
		for(uint64_t i = N_Group_S_pad*12+5; i < (N_Group_S_pad+N_S)*12; i+=12){
			for(uint64_t j = 0; j < N_Group_T; j++){
				sum_DFE = kahanSum(sum_DFE,(double)syn[j*12*(N_Group_S_pad+N_S_pad)+i],1);
			}
		}
		printf("sum_CPU = %.8e, sum_DFE = %.8e\n",sum_CPU, sum_DFE);
		average_error = kahanSum(sum_CPU,-sum_DFE,1)/((N_S+N_Group_S)*N_Group_T);
		relative_average_error = average_error/(sum_CPU/((N_S+N_Group_S)*N_Group_T));
		printf("Synapses R: Average error= %.8e Relative Average Error= %.8e\n", average_error, relative_average_error);

		sum_CPU = 0.0;
		sum_DFE = 0.0;
		error = 0.0;
		for(uint64_t i = 0; i < N_Group_T; i++){
			for(uint64_t j = 0; j < N_Group_S+N_S; j++){
				sum_CPU = kahanSum(sum_CPU,syn_CPU[j][i].u,1);
			}
		}
		for(uint64_t i = 6; i < N_Group_S*12; i+=12){
			for(uint64_t j = 0; j < N_Group_T; j++){
				sum_DFE = kahanSum(sum_DFE,(double)syn[j*12*(N_Group_S_pad+N_S_pad)+i],1);
			}
		}
		for(uint64_t i = N_Group_S_pad*12+6; i < (N_Group_S_pad+N_S)*12; i+=12){
			for(uint64_t j = 0; j < N_Group_T; j++){
				sum_DFE = kahanSum(sum_DFE,(double)syn[j*12*(N_Group_S_pad+N_S_pad)+i],1);
			}
		}
		average_error = kahanSum(sum_CPU,-sum_DFE,1)/((N_S+N_Group_S)*N_Group_T);
		relative_average_error = average_error/(sum_CPU/((N_S+N_Group_S)*N_Group_T));
		printf("Synapses u: Average error= %.8e Relative Average Error= %.8e\n", average_error, relative_average_error);

		sum_CPU = 0.0;
		sum_DFE = 0.0;
		error = 0.0;
		for(uint64_t i = 0; i < N_Group_T; i++){
			for(uint64_t j = 0; j < N_Group_S+N_S; j++){
				sum_CPU = kahanSum(sum_CPU,syn_CPU[j][i].U,1);
			}
		}
		for(uint64_t i = 7; i < N_Group_S*12; i+=12){
			for(uint64_t j = 0; j < N_Group_T; j++){
				sum_DFE = kahanSum(sum_DFE,(double)syn[j*12*(N_Group_S_pad+N_S_pad)+i],1);
			}
		}
		for(uint64_t i = N_Group_S_pad*12+7; i < (N_Group_S_pad+N_S)*12; i+=12){
			for(uint64_t j = 0; j < N_Group_T; j++){
				sum_DFE = kahanSum(sum_DFE,(double)syn[j*12*(N_Group_S_pad+N_S_pad)+i],1);
			}
		}
		printf("sum_CPU = %.8e, sum_DFE = %.8e\n",sum_CPU, sum_DFE);
		average_error = kahanSum(sum_CPU,-sum_DFE,1)/((N_S+N_Group_S)*N_Group_T);
		relative_average_error = average_error/(sum_CPU/((N_S+N_Group_S)*N_Group_T));
		printf("Synapses U: Average error= %.8e Relative Average Error= %.8e\n", average_error, relative_average_error);

		sum_CPU = 0.0;
		sum_DFE = 0.0;
		error = 0.0;
		for(uint64_t i = 0; i < N_Group_T; i++){
			for(uint64_t j = 0; j < N_Group_S+N_S; j++){
				sum_CPU = kahanSum(sum_CPU,syn_CPU[j][i].A,1);
			}
		}
		for(uint64_t i = 8; i < N_Group_S*12; i+=12){
			for(uint64_t j = 0; j < N_Group_T; j++){
				sum_DFE = kahanSum(sum_DFE,(double)syn[j*12*(N_Group_S_pad+N_S_pad)+i],1);
			}
		}
		for(uint64_t i = N_Group_S_pad*12+8; i < (N_Group_S_pad+N_S)*12; i+=12){
			for(uint64_t j = 0; j < N_Group_T; j++){
				sum_DFE = kahanSum(sum_DFE,(double)syn[j*12*(N_Group_S_pad+N_S_pad)+i],1);
			}
		}
		average_error = kahanSum(sum_CPU,-sum_DFE,1)/((N_S+N_Group_S)*N_Group_T);
		relative_average_error = average_error/(sum_CPU/((N_S+N_Group_S)*N_Group_T));
		printf("Synapses A: Average error= %.8e Relative Average Error= %.8e\n", average_error, relative_average_error);
		fflush(stdout);
	}
	return 0;
}
