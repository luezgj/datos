// ************************************************************
// Header file for Izhikevich Neural network
// Madhavun Candadai Vasu
// May 12, 2016 - created
//
// March 11, 2020 - modified by Luciano Gimenez
//
// Other headers from Randall Beer
// ************************************************************

# pragma once

#include "VectorMatrix.h"
#include "random.h"
#include "params.h"
#include "IzhiNNArch.h"
#include <iostream>
#include <math.h>

class IzhiNN {
private:
	int t;						// current time step of the network
	int T;						//current global time step
	int size;                   // size
	//Borré a b c d (tengo a y d en la arquitectura)
	TVector<double> v, u;      // states and temp buffers

	IzhiNNArch* arch;
	Delays* del;

	TVector<double> inputs;             // external input
	TVector<int> outputs, prevOutputs;   // outputs are 0 or 1

	float input_bias;
	
	float	S[N];
	float	I[N];

	float	s[N][M], sd[N][M];		// matrix of synaptic weights and their derivatives
	float* s_pre[N][3 * M], * sd_pre[N][3 * M];		// presynaptic weights

	int		N_firings;				// the number of fired neurons 
	int		firings[N_firings_max][2]; // indeces and timings of spikes

	//int		N_firings_global;				// the number of fired neurons 
//	int		firings_global[N_firings_max_global][2]; // indeces and timings of spikes

	double	LTP[N][STDP_stepSize + 1 + D], LTD[N];		  //

	bool STDPEnabled=true;
public:
	//constructor
	IzhiNN();
	//destructor
	~IzhiNN();

	//getters and setters - size
	void setNetworkVectors(); // defined in cpp file
	void setPreWeigths(); 
	void setArch(IzhiNNArch * architecture){arch=architecture;};
	void setDelays(Delays * delays){del=delays;};

	void setInputBias(double inputBias) { input_bias = inputBias; };

	//getters and setters - states
	void setStateV(int neuronIndex, double value) { v[neuronIndex] = value; };
	double getStateV(int neuronIndex) { return v[neuronIndex]; };
	void setStateU(int neuronIndex, double value) { u[neuronIndex] = value; };
	double getStateU(int neuronIndex) { return u[neuronIndex]; };

	double getNeuronalInput(int n) { return I[n]; };
	double getSensoryInput(int n) { return S[n]; };

	//getters and setters - inputs/outputs
	void setInput(int neuronIndex, double inputVal) { inputs[neuronIndex] = inputVal; };
	double getInput(int neuronIndex) { return inputs[neuronIndex]; };
	// cannot set output directly; only via v
	//void setOutput(int neuronIndex, double outputVal){outputs[neuronIndex] = outputVal;};
	double getOutput(int neuronIndex) { return outputs[neuronIndex]; };

	//getters and setters - weights
	void setWeights(int from, int nConection, double val) { s[from][nConection] = val; };
	double getWeight(int from, int nConection) { return s[from][nConection]; };

	// network dynamics
	void eulerStep(double stepSize, double (*distanceInputs)[numDistSensors]);
	void randomInitNeuronStates(RandomState& rs);
	void refractoryInitNeuronStates();
	void randomInitNetwork();

	//reseters
	void resetInputs();
	void resetDerivative();
	void resetTime(){t=0; T=0;};
	void resetPlasticity();

	void saveSpikes(const char * filename);


	float getPreWeight(int neuron, int connection){return *s_pre[neuron][connection];};

	// set STDP enabled (true) or disabled(false)
	void changeSTDPStatus(bool status) { STDPEnabled=status; };
};