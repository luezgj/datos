// ************************************************************
// Definition file for Agent's nervous system
// Madhavun Candadai Vasu
// Jun 07, 2016 - created
//
// March 11, 2020 - modified by Luciano Gimenez
//
// ************************************************************

#include "NervousSystem.h"
#define SYM

//TODO: setters for peripheral nervous system
// init cns
void NervousSystem::initCNS(genotype gen,Delays * delays) {
	int from;
	arch.initArch(gen.ints,delays);

	cns.setNetworkVectors();
	cns.setArch(&arch);
	cns.setDelays(delays);

	int genIndex=1;

	for (from = 0; from < (Ne/2); from++) //excitatory connections
	{
		for (int connection = 0; connection < M; connection++)
		{
			cns.setWeights(from, connection, gen.reals[genIndex]);
			int simetrico = IzhiNNArch::sim(from);
			cns.setWeights(simetrico, connection, gen.reals[genIndex]);
			genIndex++;
		}
	}

	for (from = Ne; from < Ne+(Ni/2); from++) //inhibitory connections
	{
		for (int connection = 0; connection < M; connection++)
		{
			cns.setWeights(from, connection, inhWeigth);
			cns.setWeights(IzhiNNArch::sim(from), connection, inhWeigth);
		}
	}

	for (from = Ne+Ni; from < (Ne+Ni+(numDistSensors/2)); from++) //sensory to hidden connections
	{
		for (int connection = 0; connection < inConnections; connection++)
		{
			cns.setWeights(from, connection, gen.reals[genIndex]);
			cns.setWeights(IzhiNNArch::sim(from), connection, gen.reals[genIndex]);
			genIndex++;
		}
	}
	//neurona sensora del medio
	if (numDistSensors % 2 == 1) for (int connection = 0; connection < inConnections; connection++) {
		from = Ne + Ni + (numDistSensors / 2);
		cns.setWeights(from, connection, gen.reals[genIndex]);
		cns.setWeights(from, ++connection, gen.reals[genIndex]);
		genIndex++;
	}
	//las motoras son siempre dos, entonces no hace falta manejar distinto la del medio (que no existe)
	for (from = Ne+Ni+numDistSensors; from < Ne+Ni+numDistSensors+(numMotorNeurons/2) ; from++) //hidden to motors connections
	{
		for (int connection = 0; connection < outConnections; connection++)
		{
			cns.setWeights(from, connection, gen.reals[genIndex]);
			cns.setWeights(IzhiNNArch::sim(from), connection, gen.reals[genIndex]);
			genIndex++;
		}
	}

	// init states
	cns.refractoryInitNeuronStates();

	// init windows
	for (int i = 1; i <= windowSize; i++) LMWindow.push(0.);
	LMWinsum = 0;
	for (int i = 1; i <= windowSize; i++) RMWindow.push(0.);
	RMWinsum = 0;
}


double NervousSystem::step(double stepSize, TVector<double> distanceInputs) {
	//cout << "\t\t in NervousSystem step" << endl;

	// step izhiNN
	cns.eulerStep(stepSize, distanceInputs);

	// update windows and return difference in rate
	updateRate();
	//cout << "\t\tRMWINSUM " << RMWinsum << " LMWINSUM " << LMWinsum;
	return ((RMWinsum - LMWinsum) / windowSize)+outputBias;
}

void NervousSystem::updateRate() {
	double outSumL, outSumR;

	// compute weighted sum from inter
	outSumL = 0.; outSumR = 0.;
	int leftIndex = Ne + Ni + numDistSensors;   //index of 
	int rigthIndex = leftIndex + 1;				//the neurons
	for (int mm = 0; mm < M; mm++) { //all conections to motors
		int LInterneuron = arch.getPre(leftIndex, mm);
		int RInterneuron = arch.getPre(rigthIndex, mm);
		double LWeight = cns.getWeights(leftIndex, mm);
		double RWeight = cns.getWeights(rigthIndex, mm);
		if (LInterneuron>=Ne) LWeight=(-LWeight);
		if (RInterneuron>=Ne) RWeight=(-RWeight);
		outSumL += cns.getOutput(LInterneuron) * LWeight;
		outSumR += cns.getOutput(RInterneuron) * RWeight;
	}

	//update LM
	LMWinsum -= LMWindow.front();
	LMWindow.pop();
	LMWinsum += outSumL;
	LMWindow.push(outSumL);

	// update RM
	RMWinsum -= RMWindow.front();
	RMWindow.pop();
	RMWinsum += outSumR;
	RMWindow.push(outSumR);
}

void NervousSystem::reset() {
	// reset neuron states
	cns.refractoryInitNeuronStates();
	cns.resetInputs();
	cns.resetDerivative();
	cns.resetTime();

	// reset RateCode windows
	int s = LMWindow.size();
	for (int i = 1; i <= s; i++) {
		LMWindow.pop();
		LMWindow.push(0.);
		RMWindow.pop();
		RMWindow.push(0.);
	}
	LMWinsum = 0.;
	RMWinsum = 0.;
}