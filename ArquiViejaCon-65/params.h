// ************************************************************
// Header file that contains all params and consts
// Madhavun Candadai Vasu
// Jun 30, 2016 - created
//
// ************************************************************
# pragma once
# define SYM

// Simulation Parameters
const double integStepSize = 0.1;           // integration step size
const int TIME_CONSTANT = 1;

struct categEvalParams {
	const int numTrials = 24;                     // number of trials per category
	const double maxDistance = 45;               // maxDistance between agent and object when oY = 0
};

const	int		M = 2/*5*/;		// the number of synapses per neuron   (10% de N)
const	int		D = 2/*20*/;		// maximal axonal conduction delay    CAPAZ CONVENGA BAJARLO UN POCO (si no simulo en cada segundo los ms)

const double inhWeigth = -5.0;    //fixed inhibitory wheigth

// Agent Params

const	int		Ne = 8;		// excitatory neurons			
const	int		Ni = 2;		// inhibitory neurons

const int cnsSize = Ne + Ni;                      // size of central nervous system - number of interneuronsss

const int numDistSensors = 7;
const int numInterNeurons = cnsSize;

const int numMotorNeurons = 2;		//motor neurons

const int evolvableNeurons = ceil(float(Ne) / 2) + ceil(float(numDistSensors) / 2) + ceil(float(numMotorNeurons) / 2);

const	int	N = Ne + Ni + numDistSensors + numMotorNeurons;		// total number of neurons	

const int genotypeSize = evolvableNeurons * M + 1 + 1 +1 /*ONLY FOR ACKLEY +3*/;
// exitatoryConnectionWeigths+distSensorToInter+interToLM+interToRM+rateCodeGain+windowSize+inputBias
//const int windowSize = 15/integStepSize;

const int N_firings_max = 100 * N;	// upper limit on the number of fired neurons per trial

// Environment parms.
const double DIAMETER = 30;
const double ENVWIDTH = 400;
const double OBJVEL = 3.;
const double ENVHEIGHT = 265;
const double SENSORDISTANCE = 250;

// Evolution Params
const int popSize = 100;//100;
const int maxGens = 5;//3000;
const double mutationVariance = 0.5;
const double crossoverProb = 0.5;

// IzhiNN
const int SPIKEV = 30; 