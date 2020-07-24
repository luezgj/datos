// ************************************************************
// Header file for Agent's nervous system
// Madhavun Candadai Vasu
// Jun 07, 2016 - created
//
// March 11, 2020 - modified by Luciano Gimenez
//
// TODO: Make atrributes private
// ************************************************************

#pragma once

#include "IzhiNN.h"
#include "IzhiNNArch.h"
#include <queue>

class NervousSystem {
private:
	IzhiNN cns;                                     //central nervous sytem
	IzhiNNArch arch;
	int cnsSize, windowSize;
	queue<double> LMWindow, RMWindow;       // window of activity for each motor
	double LMWinsum, RMWinsum;                      // sum of activity in queue

	double outputBias;

public:
	NervousSystem() {};
	~NervousSystem() {};

	void setWindowSize(int winSize) { windowSize = winSize; };

	void setInputBias(double inputBias) { cns.setInputBias(inputBias); };

	// init cns
	void initCNS(genotype gen,Delays * delays);
	double step(double stepSize, TVector<double> distanceInputs);
	void updateRate();
	void reset();
	double getStateV(int n) { return cns.getStateV(n); };
	double getInput(int n) { return cns.getNeuronalInput(n); };
	double getSensoryInput(int n) { return cns.getSensoryInput(n); };
	void saveSpikes(const char* filename) { cns.saveSpikes(filename); };

	void setOutputBias(double bias) { outputBias=bias; };
}; 