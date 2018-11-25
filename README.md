# AdexSimMaxeler
Simulation of the Adaptive Exponent Integrate-and-Fire Neuron model with Spike-timing Dependent Plasticity accelerated on Maxeler Dataflow Engines

## Getting Started

The following instructions help build and execute the simulation on a Maxeler DFE.

### Prerequisites

The program has been compiled with MaxCompiler 2017.2.1 and targeted on a MAIA DFE Model.

## Compilation

The following instructions are needed to compile the program.

### For Simulation

```
export PATH=$HOME/bin:/opt/altera/13.1/quartus/bin/:$PATH
source "/opt/maxcompiler2017.2.1/settings.sh"
make RUNRULES="Simulation"
```

### For DFE

```
cd ~/CPUCode
export PATH=$HOME/bin:/opt/altera/13.1/quartus/bin/:$PATH
source "/opt/maxcompiler2017.2.1/settings.sh"
make RUNRULES="DFE"
```

## Running

There are two options to run the program. The first one is running on the machine's CPU simulating the DFE operation and the second one is running on DFEs.

### For Simulation

```
source "/opt/maxcompiler2017.2.1/settings.sh"
make RUNRULES="Simulation" startsim // And then follow the instructions of the output of this instruction
./../RunRules/Simulation/Simulation Number_of_Timesteps N_S N_Group_S N_Group_T SpikeInt WriteEn AllConnected RunCPU
```
### For DFE

```
./../RunRules/Simulation/Simulation Number_of_Timesteps N_S N_Group_S N_Group_T SpikeInt WriteEn RunCPU AllConnected
```
Number of Timesteps: Number of Simulation Steps\
N_S: Number of Input Neurons\
N_Group_S: Number of AdEx Neurons as Source\
N_Group_T: Number of AdEx Neurons as Target\
SpikeInt: Input Neurons Spike Interval\
WriteEn: 1: Neuron Values output for every timestep, 0: Output only of the final step\
AllConnected: 1: All connected, 0: Connections read from "connections.txt"\
RunCPU: 1: Run CPU Code with doubles to check acceleration and errors, 0: Only DFE run\


