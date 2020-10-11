// ************************************************************
// Main for evolving categorizing Agent
// Madhavun Candadai Vasu
// Jun 20, 2016 - created
//
// ************************************************************
#include <iostream>
#include <fstream>
#include <stdlib.h>
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

static string ExperimentSufix="";

Delays* delays; 


const	int	W=2;	// initial width of polychronous groups
int		min_group_path = 5;		// minimal length of a group
int		min_group_time = 5;	// minimal duration of a group (ms)


//************
// Fitness
//************
void convertGenotypeToPhenotype(genotype g, Agent& a) {
	//cout << "Converting genotype to phenotype" << endl;

	double input_bias = MapSearchParameter(g.reals[genotypeSize - 2], -4, -2);
	a.setInputBias(input_bias);

	double output_bias = MapSearchParameter(g.reals[genotypeSize - 3], -10, 10);
	a.setOutputBias(output_bias);

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
			nsParams[index] = MapSearchParameter(g.reals[index], 0.0, 10.0);   //Creo que izquierdo busca entre 0 y 2
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

	for(int presentation= 1; presentation<= presentationsNum; presentation++){
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
			while (oY > DIAMETER / 2) {
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

			if (presentation == presentationsNum) {
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
						spikesName = ExperimentSufix +"_spikes_c" + to_string(trial) + ".dat";
					else
						spikesName = ExperimentSufix +"_spikes_l" + to_string(trial- cEvalParams.numTrials) + ".dat";

					a.saveSpikes(spikesName.c_str());
				}
			}
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
	BestIndividualFile.open("./"+ ExperimentSufix +"_best.gen.dat");
	BestIndividualFile << setprecision(32);
	BestIndividualFile << bestGen.reals << endl;
	BestIndividualFile << bestGen.ints << endl;
	BestIndividualFile.close();

	// write out best phenotype
	ofstream bestPhenotypeFile;
	bestPhenotypeFile.open("./"+ ExperimentSufix +"_best.phen.dat");
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
	bestCategFile.open("./"+ ExperimentSufix +"_bestCateg.dat");
	cout.rdbuf(bestCategFile.rdbuf());
	//cEvalParams.numTrials = 24;
	writeFlag = 3;
	RandomState rs;
	evaluateFitness(bestGen, rs);
	bestCategFile.close();
}






int			N_polychronous;


double		C_rel = 0.95*C_max;
const int	polylenmax = N;
int			N_postspikes[polylenmax], I_postspikes[polylenmax][N], J_postspikes[polylenmax][N], D_postspikes[polylenmax][N], L_postspikes[polylenmax][N];
double		C_postspikes[polylenmax][N];
int			N_links, links[2*W*polylenmax][4];
int			group[polylenmax], t_fired[polylenmax], layer[polylenmax];
int			gr3[W], tf3[W];
int			I_my_pre[3*M], D_my_pre[3*M], N_my_pre;
int			N_fired;


FILE		*fpoly;

const	int	latency = D; // maximum latency 


void	polychronous(int nnum, Agent & a)
{
	int	i,j, t, p, k;
	int npre[W];
	int dd;
	int	t_last, timing;
	int	Dmax, L_max; 
	int	used[W], discard;

	double v[N],u[N],I[N];
	

	N_my_pre = 0;
	for (i=0;i<a.getNPre(nnum);i++)
	if (a.getPreWeight(nnum,i) > C_rel) 
	{
		I_my_pre[N_my_pre]=a.getPre(nnum,i);
		D_my_pre[N_my_pre]=a.getDPre(nnum,i);
		N_my_pre++;
	}
	if (N_my_pre<W) return;

	for (i=0;i<W;i++)	npre[i]=i;

	while (0==0) 
	{
		Dmax=0;
		for (i=0;i<W;i++) if (Dmax < D_my_pre[npre[i]]) Dmax=D_my_pre[npre[i]];
		
		for (i=0;i<W;i++)
		{
			group[i]=I_my_pre[npre[i]];
			t_fired[i]= Dmax-D_my_pre[npre[i]];
			layer[i]=1;
			
			for (dd=0; dd<D; dd++)		 
			for (j=0; j<delays->getDelayLength(group[i],dd); j++)         
			{
				p = a.getPost(group[i],delays->getConnNumberByDelay(group[i],dd,j));  
				if ((a.getWeigth(group[i],delays->getConnNumberByDelay(group[i],dd,j)) > C_rel) & (dd>=D_my_pre[npre[i]]))
				{
					timing = t_fired[i]+dd+1;
					J_postspikes[timing][N_postspikes[timing]]=group[i];				// presynaptic
					D_postspikes[timing][N_postspikes[timing]]=dd;						// delay
					C_postspikes[timing][N_postspikes[timing]]=a.getWeigth(group[i],delays->getConnNumberByDelay(group[i],dd,j));	// syn weight
					I_postspikes[timing][N_postspikes[timing]++]=p;						// index of post target	
				}
			}
		}

		for (i=0;i<N;i++) {v[i]=-70; u[i]=0.2*v[i]; I[i]=0;};

		N_links = 0;
		N_fired=W;
		t_last = D+D+latency+1;
		t=-1;
		while ((++t<t_last) & (N_fired < polylenmax))
		{    
			for (p=0;p<N_postspikes[t];p++) 
			  I[I_postspikes[t][p]]+=C_postspikes[t][p]; 
 
		  	for (i=0;i<N;i++)
			{
				v[i]+=0.5*((0.04*v[i]+5)*v[i]+140-u[i]+I[i]);
				v[i]+=0.5*((0.04*v[i]+5)*v[i]+140-u[i]+I[i]);
				u[i]+=a[i]*(0.2*v[i]-u[i]);
				I[i]=0;
			}

			for (i=0;i<N;i++) 
			if (v[i]>=30)
			{
				v[i] = -65;
				u[i]+=d[i];

				if (N_fired < polylenmax)
				{
					t_fired[N_fired]= t;
					group[N_fired++]=i;
					for (dd=0; dd<D; dd++)
					for (j=0; j<delays->getDelayLength(i,dd); j++)
					if (( a.getWeigth(i,delays->getConnNumberByDelay(i,dd,j)) > C_rel) | (i>=Ne)) 
					{
						timing = t+dd+1;
						J_postspikes[timing][N_postspikes[timing]]=i;				// presynaptic
						D_postspikes[timing][N_postspikes[timing]]=dd;				// delay
//						L_postspikes[timing][N_postspikes[timing]]=NL+1;			// layer
						C_postspikes[timing][N_postspikes[timing]]=a.getWeigth(i,delays->getConnNumberByDelay(i,dd,j));	   // syn weight
						I_postspikes[timing][N_postspikes[timing]++]=a.getPost(i,delays->getConnNumberByDelay(i,dd,j));// index of post target	
					}
					if (t_last < timing+1) 
					{
						t_last = timing+1;
						if (t_last > polylenmax-D-1) t_last = polylenmax-D-1;
					}
				}
			}
		}
		
		if (N_fired>2*W)
		{
			N_links=0;
			L_max=0;
			for (i=W;i<N_fired;i++)
			{
				layer[i]=0;
				for (p=t_fired[i]; (p>t_fired[i]-latency) & (p>=0); p--)
				for (j=0;j<N_postspikes[p];j++)
				if ((I_postspikes[p][j]==group[i]) & (J_postspikes[p][j]<Ne)) 
				{
				   for (k=0;k<i;k++)
				   if ((group[k]==J_postspikes[p][j]) & (layer[k]+1>layer[i])) layer[i]=layer[k]+1;
				   {
					   links[N_links][0]=J_postspikes[p][j];
					   links[N_links][1]=I_postspikes[p][j];
					   links[N_links][2]=D_postspikes[p][j];
					   links[N_links++][3]=layer[i];
					   if (L_max < layer[i]) L_max = layer[i]; 
				   }
				}
			}
										 
			discard = 0;
			for (i=0;i<W;i++)
			{
				used[i]=0;
				for (j=0;j<N_links;j++) if ((links[j][0] == group[i]) & (links[j][1] < Ne)) used[i]++;
				if (used[i] == 1) discard = 1;
			}

//			if ((discard == 0) & (t_fired[N_fired-1] > min_group_time) )  // (L_max >= min_group_path))
			if ((discard == 0) & (L_max >= min_group_path))
			{

				for (i=0;i<W;i++) {gr3[i]=group[i]; tf3[i]=t_fired[i];};

				N_polychronous++;
				cout << "\ni= " << nnum << ", N_polychronous= " << N_polychronous << ", N_fired = " << N_fired << ", L_max = " << L_max << ", T=" << t_fired[N_fired-1];
//				fprintf(fpoly, " %d  %d,       ", N_fired, L_max);
//				for (i=0; i<N_fired; i++)
//					fprintf(fpoly, " %d %d, ", group[i], t_fired[i]);
//				fprintf(fpoly, "        ");
//				for (j=0;j<N_links;j++)
//				   fprintf(fpoly, " %d %d %d %d,  ", links[j][0], links[j][1], links[j][2], links[j][3]);
//				fprintf(fpoly, "\n");
			}
		}

  		for (dd=Dmax;dd<t_last;dd++) N_postspikes[dd]=0;
		if (t_last == polylenmax-D) for (dd=t_last;dd<polylenmax;dd++) N_postspikes[dd]=0;

		i=1;
		while (++npre[W-i] > N_my_pre-i) if (++i > W) return; 
		while (i>1) {npre[W-i+1]=npre[W-i]+1; i--;}
	}	
}


void	all_polychronous()
{
	int	i;
	N_polychronous=0;
	fpoly = fopen("..//polyall.dat","w");
   	for (i=0;i<polylenmax;i++) N_postspikes[i]=0;

	for (i=0;i<Ne;i++) polychronous(i);

	cout << "\nN_polychronous=" << N_polychronous << "\n";
	fclose(fpoly);
}







int main(int argc, char* argv[]) {

#ifdef SYM
	cout << "Symmetric network!" << endl;
#endif

	delays = new Delays();

	float crossoverProbs[2] = {0.60, 0.70};
	float mutationVariances[2] = {5 , 10};
	int  TourSizes[2] = {10, 15};
	TCrossoverArchWeight crossArcWeights[2] = {TOGETHER, APART};
	TCrossoverMode crossModes[2] = {UNIFORM ,TWO_POINT};
	
	int creepSteps[2] = {5, 10};

	if (argc != 7) {
		cout << "error. Tiene que haber 6 argumentos";
		return -1;
	}

	int crossoverProbsI= stoi(argv[1]);
	int mutationVariancesI = stoi(argv[2]);
	int TourSizesI = stoi(argv[3]);
	int crossArcWeightsI = stoi(argv[4]);
	int crossModesI = stoi(argv[5]);
	int experiment = stoi(argv[6]);
	
	//int creepStepI = stoi(argv[7]);
	
	ExperimentSufix = "CP" + to_string(crossoverProbsI) +
		"_MV" + to_string(mutationVariancesI) +
		"_TO" + to_string(TourSizesI) +
		"_AW" + to_string(crossArcWeightsI) +
		"_CM" + to_string(crossModesI) +
		"_EX" + to_string(experiment);// +
		//"_CP" + to_string(creepStepI);

	cout << ExperimentSufix << endl;

	
	TSearch s(genotypeSize, archGenotypeSize);

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
	paramsFile.open("./" + ExperimentSufix + "_paramsFile.dat");
	paramsFile << randomseed << endl;
	paramsFile << popSize << endl;
	paramsFile << maxGens << endl;
	paramsFile << mutationVariances[mutationVariancesI] << endl;
	paramsFile << crossoverProbs[crossoverProbsI] << endl;
	paramsFile << cnsSize << endl;
	//paramsFile << windowSize << endl;
	paramsFile.close();

	// display functions
	s.SetPopulationStatisticsDisplayFunction(EvolutionaryRunDisplay);
	s.SetSearchResultsDisplayFunction(ResultsDisplay);

	s.SetSelectionMode(RANK_BASED);
	s.SetReproductionMode(GENETIC_ALGORITHM);
	s.SetParentSelectionMode(TOURNAMENT);
	s.SetTournamentSize(TourSizes[TourSizesI]);
	s.SetParentsTournamentSize(TourSizes[TourSizesI]);

	s.SetArchMutationMode(RANDOM_RESET);
	//s.setCreepMutationStep(crossModes[creepStepI]);

	s.SetCrossoverArchWeight(crossArcWeights[crossArcWeightsI]);
	s.SetPopulationSize(popSize);
	s.SetMaxGenerations(maxGens);
	s.SetMutationVariance(mutationVariances[mutationVariancesI]);
	s.SetArchRandMutationProb(mutationVariances[mutationVariancesI] * 0.01);
	s.SetCrossoverProbability(crossoverProbs[crossoverProbsI]);
	s.SetCrossoverMode(crossModes[crossModesI]);
	s.SetMaxExpectedOffspring(1.9);
	s.SetElitistFraction(0.1);
	s.SetSearchConstraint(1);
	s.SetCheckpointInterval(0);
	s.SetReEvaluationFlag(0);

	// redirect standard output to a file
	ofstream evolfile;
	evolfile.open("./" + ExperimentSufix + "_fitness.dat");
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