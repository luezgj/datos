// ************************************************************
// Header file that contains all params and consts
// Madhavun Candadai Vasu
// Jun 30, 2016 - created
//
// ************************************************************
# pragma once
# define SYM
#include <cmath>

// Simulation Parameters
const double integStepSize = 0.1;           // integration step size
const int TIME_CONSTANT = 1;

const int numTrials = 12;

struct categEvalParams {
	const int numTrials = 12;                     // number of trials per category
	const double maxDistance = 45;               // maxDistance between agent and object when oY = 0
};
const int presentationsNum=4; //number of presenting set of trials (including the set of evaluation)

//Izhikevich model:
const 	float	AExc= 0.02;   	// a parameter- izh neuron model
const 	float 	AInh= 0.1;		// a parameter- izh neuron model

const 	float	DExc= 8.0;		// d parameter- izh neuron model
const 	float 	DInh= 2.0;		// d parameter- izh neuron model

const	int		M =100/*5*/;		// the number of synapses per neuron   (10% de N)
const	int		D =20/*20*/;		// maximal axonal conduction delay    

const 	int 	inConnections= 100;   //number of connections from each sensor
const 	int 	outConnections= 100;   //number of connections to each motor

const double inhWeigth = -5.0;    //fixed inhibitory wheigth

const	float	sm = 10.0;		// maximal synaptic strength	

// Agent Params

const	int		Ne =800;		// excitatory neurons			
const	int		Ni =200;		// inhibitory neurons

const int numDistSensors = 7;

const int cnsSize = Ne + Ni + numDistSensors;                      // size of central nervous system - number of interneuronsss

const int numInterNeurons = cnsSize;

const int numMotorNeurons = 2;		//motor neurons

const int evolvableWeightNeurons = ceil(float(Ne) / 2) + ceil(float(numDistSensors) / 2) + ceil(float(numMotorNeurons) / 2);

const	int	N = Ne + Ni + numDistSensors + numMotorNeurons;		// total number of neurons	


//real genotype      weigths and bias/rates
const int genotypeSize = evolvableWeightNeurons * M + 1 + 1 +1 +1 /*ONLY FOR ACKLEY +3*/;
// exitatoryConnectionWeigths+distSensorToInter+interToLM+interToRM+rateCodeGain+windowSize+inputBias+outputBias

//integer genotype      architecture
const int archGenotypeSize=ceil(float(cnsSize) / 2) * M + ceil(float(numDistSensors) / 2)* inConnections + ceil(float(numMotorNeurons) / 2) * outConnections;

//const int N_firings_max_global = 100 * N;	// upper limit on the number of fired neurons per trial

const int N_firings_max = 100 * N;	// upper limit on the number of fired neurons per trial

const int STDP_stepSize = 100; // steps to update weigths by STDP

// Environment parms.
const double DIAMETER = 30;
const double ENVWIDTH = 400;
const double OBJVEL = 3.;
const double ENVHEIGHT = 265;
const double SENSORDISTANCE = 250;

// Evolution Params
const int popSize = 200;//200;//200;//100;
const int maxGens = 850;//1500;// 2000;//1200;//2500;//3000;
const double mutationVariance = 5 /* Lo aumenté xq tengo muchos más parametros, y si lo reparto me queda casi nada en cada uno0.5*/;
const double crossoverProb = 0.65;

const int PARENTS_TOURNAMENT_SIZE = 12;
const int TOURNAMENT_SIZE = 12;

// Architecture evolution
const double randomMutProb=0.05;

// IzhiNN
const int SPIKEV = 30; 