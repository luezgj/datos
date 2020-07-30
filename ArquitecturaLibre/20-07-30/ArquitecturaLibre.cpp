// ************************************************************
// Main for evolving categorizing Agent
// Madhavun Candadai Vasu
// Jun 20, 2016 - created
//
// ************************************************************
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <math.h>
#include "Agent.h"
#include "TSearch.h"
#include "Environment.h"
#include "IzhiNNArch.h"
#include "params.h"

categEvalParams cEvalParams;

//Flags
int writeFlag = 0;

Delays* delays; //or load from a file

//************
// Fitness
//************
void convertGenotypeToPhenotype(genotype g, Agent& a) {
	//cout << "Converting genotype to phenotype" << endl;

	double input_bias = MapSearchParameter(g.reals[genotypeSize - 2], -4, -2);
	a.setInputBias(input_bias);

	double output_bias = MapSearchParameter(g.reals[genotypeSize - 3], -100, 100);
	//a.setOutputBias(output_bias);

	// penultimate gene is rateCode
	double rateCodeGain = MapSearchParameter(g.reals[genotypeSize - 1], 1,500/*200, 1200*/);
	a.setRateCodeGain(rateCodeGain);

	// last gene is windowSize
	double windowSize = round(MapSearchParameter(g.reals[genotypeSize],/* 1, 8*/5, 25)) / integStepSize;
	a.setWindowSize(windowSize);

	TVector<double> nsParams;
	nsParams.SetBounds(1, genotypeSize);

	int index = 1;
	for (int evNeuron = 0; evNeuron < evolvableWeightNeurons; evNeuron++) {
		for (int connection = 0; connection < M; connection++) {
			nsParams[index] = MapSearchParameter(g.reals[index], 0.0, 20.0);   //Creo que izquierdo busca entre 0 y 2
			index++;
		}
	}

	genotype genConverted;
	genConverted.reals = nsParams;
	genConverted.ints = g.ints;

	//cout << "\tFinished scaling weights" << endl;
	a.initNS(genConverted, delays);

	if (writeFlag == 1) {
		for (int i = 1; i <= g.ints.Size(); i++)
			cout << g.ints[i] << endl;
		for (int i = 1; i <= nsParams.Size(); i++)
			cout << nsParams[i] << endl;
		cout << output_bias << endl;
		cout << input_bias << endl;
		cout << rateCodeGain << endl;
		cout << windowSize << endl;
	}
	//cout << "Finished converting genotype to phenotype.." << endl;
}

double evaluateCateg(Agent& a) {
	//cout << "\tIn evaluateCateg" << endl;
	Environment env(DIAMETER, ENVWIDTH, OBJVEL, ENVHEIGHT, SENSORDISTANCE);
	double fitness = 0., fit1 = 0., fit2 = 0.;

	// vector of starting positions for object x-coordinate
	//para agente asimetrico
	/*TVector<double> startingPositions;
	startingPositions.SetBounds(1, cEvalParams.numTrials * 2);
	for (int i = 1; i <= cEvalParams.numTrials; i++) {
		double posOffset = (i - 1) * (100. / (cEvalParams.numTrials - 1)) - 50.;
		startingPositions[i] = posOffset;                           // for isCircle trials
		startingPositions[i + cEvalParams.numTrials] = posOffset;   // for !isCircle trials
	}*/
	//cEvalParams.maxDistance = std::abs(startingPositions[1]);

	// vector of starting positions for object x-coordinate
	//para agente simetrico
	//para agentes simétricos
	TVector<double> startingPositions;
	startingPositions.SetBounds(1, cEvalParams.numTrials * 2);
	for (int i = 1; i <= cEvalParams.numTrials; i++) {
		double posOffset = (i - 1) * (100. / (cEvalParams.numTrials * 2 - 1)) - 50.;
		startingPositions[i] = posOffset;                           // for isCircle trials
		startingPositions[i + cEvalParams.numTrials] = posOffset;   // for !isCircle trials
	}
	//cEvalParams.maxDistance = std::abs(startingPositions[1]);

	// externalInputs
	TVector<double> distanceInputs;
	distanceInputs.SetBounds(1, numDistSensors);
	distanceInputs.FillContents(0.);

	double agentPos = 0.;
	// object params
	double oY, oX;
	bool isCircle = 1;

	for (int trial = 1; trial <= cEvalParams.numTrials * 2; trial++) {
		//cout << "\tTrial-" << trial << endl;
		if (trial == cEvalParams.numTrials + 1) isCircle = !isCircle;

		// reset inputs
		distanceInputs.FillContents(0.);
		// reset agent
		a.reset();

		double t = integStepSize;
		oY = ENVHEIGHT;
		oX = startingPositions[trial];
		agentPos = 0.;
		// present object and step agent
		while (oY > 0.) {
			env.getDistanceInputs(distanceInputs, isCircle, agentPos, oX, oY);

			double offset = 0.;

			offset = a.step(integStepSize, distanceInputs);
			agentPos += offset;
			//cout << "\t\toY=" << oY << " oX=" << oX << " agentPos=" << agentPos << endl;
			if (writeFlag == 3) {
				cout << t << " " << isCircle << " " << oX << " " << oY << " " << agentPos << " " << offset / integStepSize << " ";
				for (int d = 1; d <= numDistSensors; d++) {
					cout << distanceInputs[d] << " ";
				}
				/*
				for (int n = 0; n < cnsSize; n++) {
					cout << a.getNeuronStateV(n) << " ";
				}*/
				cout << endl;
			}
			// Update Object's position
			oY = oY - OBJVEL * integStepSize;
			t += integStepSize;
		}

		// Finished presenting object
		double distance = std::abs(agentPos - oX);

		if (distance > cEvalParams.maxDistance) { distance = cEvalParams.maxDistance; }
		distance = distance / cEvalParams.maxDistance;
		//cout << "\tFi1 " << fit1 << "fit2 " << fit2 << endl;
			// Compute Fitness
		if (isCircle) { fit1 += (1 - distance); }
		else { fit2 += distance; }

		if (writeFlag == 3) {
			string spikesName;
			if (isCircle)
				spikesName = "spikes_c" + to_string(trial) + ".dat";
			else
				spikesName = "spikes_l" + to_string(trial- cEvalParams.numTrials) + ".dat";

			a.saveSpikes(spikesName.c_str());
		}
	}

	//fitness = (fit1/cEvalParams.numTrials)*(fit2/cEvalParams.numTrials);
	fitness = (fit1 + fit2) / (cEvalParams.numTrials * 2);
	// Average over trials
	return fitness;///(cEvalParams.numTrials*2);
}

double evaluateFitness(genotype& g, RandomState& rs) {
	//cout << "evaluateFitness" << endl;
	// init
	double fitCateg = 0.;

	Agent a;
	convertGenotypeToPhenotype(g, a);
	//cout << "converted Genotype To Phenotype" << endl;
	a.reset();
	fitCateg = evaluateCateg(a);

	return fitCateg;
}

//*************
// Display
//*************

void EvolutionaryRunDisplay(int generation, double bestPerf, double avgPerf, double perfVar, double executionTime) {
	cout << generation << " " << bestPerf << " " << avgPerf << " " << perfVar << " ";
	cout << executionTime << "s";
	cout << endl;
}

void ResultsDisplay(TSearch& s) {

	genotype bestGen;
	ofstream BestIndividualFile;

	bestGen = s.BestIndividual();
	BestIndividualFile.open("./best.gen.dat");
	BestIndividualFile << setprecision(32);
	BestIndividualFile << bestGen.reals << endl;
	BestIndividualFile << bestGen.ints << endl;
	BestIndividualFile.close();

	// write out best phenotype
	ofstream bestPhenotypeFile;
	bestPhenotypeFile.open("./best.phen.dat");
	cout.rdbuf(bestPhenotypeFile.rdbuf());
	writeFlag = 1;
	Agent a;
	convertGenotypeToPhenotype(bestGen, a);
	bestPhenotypeFile.close();

	/*ofstream bestBehaviorFile;
	bestBehaviorFile.open("./bestBehavior.dat");
	cout.rdbuf(bestBehaviorFile.rdbuf());
	writeFlag = 2;
	RandomState rs;
	evaluateFitness(bestVector, rs);
	bestBehaviorFile.close();*/

	ofstream bestCategFile;
	bestCategFile.open("./bestCateg.dat");
	cout.rdbuf(bestCategFile.rdbuf());
	//cEvalParams.numTrials = 24;
	writeFlag = 3;
	RandomState rs;
	evaluateFitness(bestGen, rs);
	bestCategFile.close();
}

int main() {

#ifdef SYM
	cout << "Symmetric network!" << endl;
#endif
	TSearch s(genotypeSize, archGenotypeSize);

	delays = new Delays();

	cout << "Genotype Size - " << genotypeSize << endl;
	cout << "Arch Genotype Size - " << archGenotypeSize << endl;
	// config search
	long randomseed = static_cast<long>(time(NULL));
	//long randomseed = 233456;
	s.SetRandomSeed(randomseed);
	s.SetEvaluationFunction(evaluateFitness);
	////cout << "RandomSeed - " << randomseed << endl;

	SetRandomSeed(randomseed);

	//write params to file
	ofstream paramsFile;
	paramsFile.open("./paramsFile.dat");
	paramsFile << randomseed << endl;
	paramsFile << popSize << endl;
	paramsFile << maxGens << endl;
	paramsFile << mutationVariance << endl;
	paramsFile << crossoverProb << endl;
	paramsFile << cnsSize << endl;
	//paramsFile << windowSize << endl;
	paramsFile.close();

	// display functions
	s.SetPopulationStatisticsDisplayFunction(EvolutionaryRunDisplay);
	s.SetSearchResultsDisplayFunction(ResultsDisplay);

	s.SetSelectionMode(RANK_BASED);
	s.SetReproductionMode(GENETIC_ALGORITHM);
	s.SetParentSelectionMode(TOURNAMENT);
	s.SetPopulationSize(popSize);
	s.SetMaxGenerations(maxGens);
	s.SetMutationVariance(mutationVariance);
	s.SetCrossoverProbability(crossoverProb);
	s.SetCrossoverMode(UNIFORM);
	s.SetMaxExpectedOffspring(1.9);
	s.SetElitistFraction(0.1);
	s.SetSearchConstraint(1);
	s.SetCheckpointInterval(0);
	s.SetReEvaluationFlag(0);

	// redirect standard output to a file
	ofstream evolfile;
	evolfile.open("./fitness.dat");
	cout.rdbuf(evolfile.rdbuf());

	// start evolution
	s.ExecuteSearch();

	evolfile.close();

	// back to old buf
	cout.rdbuf();
	//cout << "Finished Execution" << endl;





	/*
	arch = new IzhiNNArch();
	arch->loadFile("./arch.dat");
	writeFlag = 3;
	Agent ag;

	TVector<double> bestVector;
	ifstream BestIndividualFile;
	bestVector.SetSize(genotypeSize);


	BestIndividualFile.open("./best.gen.dat");
	BestIndividualFile >> setprecision(32);
	for (int i = 1; i <= genotypeSize; i++)
		BestIndividualFile >> bestVector[i];
	BestIndividualFile.close();

	convertGenotypeToPhenotype(bestVector, ag);

	// redirect standard output to a file
	ofstream evolfile;
	evolfile.open("./categorizationsym.dat");
	cout.rdbuf(evolfile.rdbuf());

	trySymmetry(ag);


	evolfile.close();

	// back to old buf
	cout.rdbuf();
	*/

}