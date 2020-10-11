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
	double step(double stepSize, double (*distanceInputs)[numDistSensors]);
	void updateRate();
	void reset();
	double getStateV(int n) { return cns.getStateV(n); };
	double getInput(int n) { return cns.getNeuronalInput(n); };
	double getSensoryInput(int n) { return cns.getSensoryInput(n); };
	void saveSpikes(const char* filename) { cns.saveSpikes(filename); };

	void setOutputBias(double bias) { outputBias=bias; };

	int getPost(int from, int conNumber) { return arch.getPost(from, conNumber); };
	int getPre(int from, int conNumber){return arch.getPre(from,conNumber);};
	int getNPre(int from){return arch.getNPre(from);};
	int getDPre(int from, int conNumber){return arch.getDPre(from,conNumber);};

	double getA(int nNumber) { return arch.getA(nNumber); };
	double getD(int nNumber) { return arch.getD(nNumber); };

	double getWeight(int from, int nConection) { return cns.getWeight(from,nConection); };

	float getPreWeight(int neuron, int connection){return cns.getPreWeight(neuron,connection);};

	void changeSTDPStatus(bool status) { cns.changeSTDPStatus(status); };
}; 