// ************************************************************
// Source file for Izhikevich Neural network
// Madhavun Candadai Vasu
// May 12, 2016 - created
//
// March 11, 2020 - modified by Luciano Gimenez
//
// Other headers from Randall Beer
// ************************************************************


#include "IzhiNN.h"
#include "random.h"

#define rand01 (0.9999999*double(rand())/RAND_MAX) 
#define getrandom(max1) ((rand()%(int)((max1)))) // random integer between 0 and max-1

//Constructor 
IzhiNN::IzhiNN() {
	setNetworkVectors();
	N_firings=0;
	t=0;
	T=0;
	resetInputs();
	resetDerivative();
}
//destructor
IzhiNN::~IzhiNN() {
}

//*********
//UTILS
//*********
void IzhiNN::setNetworkVectors() {
	//cout << "\t Setting izhiNN size " << ne << " " << ni << endl;
	size = cnsSize;
	//set and init other vars that depend on size
	v.SetBounds(0, size-1);
	v.FillContents(0.0);
	u.SetBounds(0, size-1);
	u.FillContents(0.0);

	inputs.SetBounds(0, size-1);
	inputs.FillContents(0.0);
	outputs.SetBounds(0, size-1);
	outputs.FillContents(0.0);
	prevOutputs.SetBounds(0, size-1);
	prevOutputs.FillContents(0.);
}

void IzhiNN::resetInputs(){
	for (int neuron = 0; neuron < N; neuron++){
		S[neuron]=0.0f;
		I[neuron]=0.0f;
	}
}

void IzhiNN::resetDerivative(){
	for (int neuron = 0; neuron < N; neuron++){
		for (int connection = 0; connection < M; connection++){
			sd[neuron][connection]=0.0f;
		}
	}	
}

void IzhiNN::setPreWeigths(){
	//Busca las conexiones presinapticas.  ignora las de las motoras xq ya son presinapticas
	for (int i = 0; i < N - numMotorNeurons; i++) {
		int conexion=0;
		for (int j = 0; j < Ne; j++)
			for (int k = 0; k < M; k++)
				if (arch->getPost(j,k) == i) {		// find all presynaptic neurons 
					s_pre[i][conexion] = &s[j][k];	// pointer to the synaptic weight	
					sd_pre[i][conexion] = &sd[j][k];// pointer to the derivative
					conexion++;
				}
	}
}

//*********
//INIT
//*********
void IzhiNN::refractoryInitNeuronStates() {
	for (int i = 0; i < cnsSize; i++) { 
		// refractory init states of spiking neurons
		v[i] = -65;
		u[i] = 0.2/*Izhikevich b[i] fijo, Candadi usaba b[i] */ * v[i];
	}
	N_firings=0;
}

void IzhiNN::randomInitNeuronStates(RandomState& rs) {
	for (int i = 0; i < size; i++) {
		// random init states
		v[i] = rs.UniformRandom(-65, 20); // refractory state to almost spiking
		u[i] = /*b[i]*/0.2 * v[i];
	}
}

void IzhiNN::resetPlasticity() {
	//	LTP = zeros(N,1001+D);
	for (int i=0;i<N;i++)
		for (int j=0;j<1+D;j++)
			LTP[i][j]=0;

	//	LTD = zeros(N,1);
	for (int i=0;i<N;i++)	LTD[i]=0;
}

/*void IzhiNN::randomInitNetwork() {
	//Do nothing
}*/


//*********
//CONTROL
//*********
void IzhiNN::eulerStep(double stepSize, double (*distanceInputs)[numDistSensors]) {
	//cout << "\t\t\tIn euler step for IzhiNN " << size << endl;
	// check for neurons that may spike
	for (int i = 0; i < size; i++) {
		////cout << "\t\t\t\t| i = " << i << "|v[i] = " << v[i];
		if (v[i] > SPIKEV) {
			// spike
			//cout << "spike de " << i;
			outputs[i] = 1;
			v[i] = -65.0/*c[i]*/;
			u[i] += arch->getD(i);
			
			if (STDPEnabled)
			{
				LTP[i][t+D]= 0.1;		
				LTD[i]=0.12;
				if (i<Ne+Ni)
				for (int j=0;j<arch->getNPre(i);j++) *sd_pre[i][j]+=LTP[arch->getPre(i,j)][t+D-arch->getDPre(i,j)-1];// this spike was after pre-synaptic spikes
			}
				
			firings[N_firings][0] = T;
			firings[N_firings++][1] = i;

			if (N_firings == N_firings_max)
			{
				cout << "*** Two many spikes, t=" << T << "*** (ignoring)";
				N_firings = 0;
			}

			
		}
		else {
			// no spike
			outputs[i] = 0;
		}
		////cout << "| o[i] = " << outputs[i] << "||";
	}
	////cout << endl;

	for (int i = 0; i < N; i++) {
		I[i] = 0.0;			// reset the input
		S[i] = 0.0;			// reset the sensory input
	}

	int talamicInNeuron = int(floor((Ne + Ni) * rand01));
	I[talamicInNeuron] = 20;

	//cout << "entrada talámica a " << talamicInNeuron << endl;

	int k = N_firings-1;
	while ((k>0) && (T - firings[k][0] < D)) {   //por todos los spikes que pasaron hace tan poco que importan
		for (int j = 0; j < del->getDelayLength(firings[k][1],T - firings[k][0]); j++)   //Por todas las conexiones que tienen un delay igual a (t actual- t spike)
		{
			int i = arch->getPost(firings[k][1],del->getConnNumberByDelay(firings[k][1],T - firings[k][0],j));
			//cout << "Hubo spike de " << firings[k][1] << " en el tiempo " << firings[k][0] << endl;
			//cout << "En el tiempo  " << t << "activa a " << i << endl;
			
			//llegó el spike a la neurona postsinaptica
			double adition = s[firings[k][1]][del->getConnNumberByDelay(firings[k][1], T - firings[k][0], j)];

			I[i] += adition;
			//cout << "Le suma " << adition<< endl;

			//STDP
			if (STDPEnabled && ( (firings[k][1] < Ne) || ( (firings[k][1] >= Ne+Ni) && (firings[k][1] < Ne+Ni+numDistSensors )) ) ) // this spike is before postsynaptic spikes
				sd[firings[k][1]][ del->getConnNumberByDelay(firings[k][1], T - firings[k][0] ,j) ] -= LTD[i];
		}
		k--;
	}	


	//ACÁ HAY QUE ACTUALIZAR LAS ENTRADAS EXTERNAS
	//for (int j = 0; j < numDistSensors; j++) {
		//Chequear que esto esté bien
		//for (int mm = 0; mm < inConnections; mm++) {			//Actualización de las neuronas conectadas al sensor
			//cout << "quiero acceder a S sub " << arch->getPost(Ne + Ni + j, mm) - 1 << endl;
			//cout << "que vale" << S[arch->getPost(Ne + Ni + j, mm) - 1] << endl;
			//cout << "quiero acceder a s sub " << Ne + Ni + j - 1 <<","<< mm - 1 << endl;
			//cout << "que vale" << s[Ne + Ni + j - 1][mm - 1] << endl;

			//S[arch->getPost(Ne + Ni + j, mm)] += s[Ne + Ni + j][mm] *  (1.0 / (1.0 + exp(-( ((*distanceInputs)[j] + input_bias)))));

			//no uso la función de activación en la entrada
			//S[arch->getPost(Ne + Ni + j, mm)] += s[Ne + Ni + j][mm] * ((*distanceInputs)[j] + input_bias);
		//}
	//}

	for (int j = 0; j < numDistSensors; j++) {
		S[Ne + Ni + j] += ( ( (*distanceInputs)[j] *2 ) + input_bias);
	}

	for (int i = 0; i < cnsSize; i++) {
		/*
		if (I[i] != 0) {
			cout << "I[ " << i << " ] = " << I[i] << endl;
		}*/

		v[i] +=/* stepSize **/ 0.5 * ((0.04 * v[i] + 5) * v[i] + 140 - u[i] + I[i] + S[i]); // for numerical stability
		v[i] += /*stepSize **/ 0.5 * ((0.04 * v[i] + 5) * v[i] + 140 - u[i] + I[i] + S[i]); // time step is 0.5 ms
		u[i] += /*stepSize **/ arch->getA(i) * (0.2 * v[i] - u[i]);

		//STDP
		if (STDPEnabled){
			LTP[i][t + D + 1] = 0.95 * LTP[i][t + D];
			LTD[i] *= 0.95;
		}
	}

	if ((T != 0) && ((T + 1) % STDP_stepSize == 0))
	{
		if (STDPEnabled) {
			for (int i = 0; i < N; i++)		// prepare for the next sec
				for (int j = 0; j < D + 1; j++)
					LTP[i][j] = LTP[i][STDP_stepSize + j];

			for (int i = 0; i < Ne; i++)	// modify only exc connections 
				for (int j = 0; j < M; j++)
				{
					sd[i][j] *= 0.9;
					s[i][j] += 0.0005 + sd[i][j];
					if (s[i][j] > sm) s[i][j] = sm;
					if (s[i][j] < 0) s[i][j] = 0.0;
				}

			for (int i = Ne + Ni; i < Ne + Ni + numDistSensors; i++) // also modify sensor connections 
				for (int j = 0; j < M; j++) {
					sd[i][j] *= 0.9;
					s[i][j] += 0.001 + sd[i][j];
					if (s[i][j] > sm) s[i][j] = sm;
					if (s[i][j] < 0) s[i][j] = 0.0;
				}
		}
		t = -1; //start new STDP step
	}

	t++;
	T++;
}


void IzhiNN::saveSpikes(const char* filename) {
	FILE* fs;
	errno_t errorCode = fopen_s(&fs, filename, "w");
	if (errorCode == 0) {
		for (int i = 0; i < N_firings; i++)
			if (firings[i][0] >= 0)
				fprintf(fs, "%d  %d\n", firings[i][0], firings[i][1]);
		fclose(fs);
	}
}

//*********
//INPUT/OUTPUT
//*********
/*
#include <iomanip>
ostream& operator<<(ostream& os, IzhiNN& inn){
	// Set the precision
	os << setprecision(32);
	int size = inn.getNumInhibitory()+inn.getNumExcitatory();
	// Write size
	os << inn.getNumExcitatory() << endl << endl;
	os << inn.getNumInhibitory() << endl << endl;
	// Write params
	for (int i = 1; i <= size; i++)
		os << inn.getParamA(i) << " ";
	os << endl << endl;
	for (int i = 1; i <= size; i++)
		os << inn.getParamB(i) << " ";
	os << endl << endl;
	for (int i = 1; i <= size; i++)
		os << inn.getParamC(i) << " ";
	os << endl << endl;
	for (int i = 1; i <= size; i++)
		os << inn.getParamD(i) << " ";
	os << endl << endl;
	// Write states
	for (int i = 1; i <= size; i++)
		os << inn.getStateV(i) << " ";
	os << endl << endl;
	for (int i = 1; i <= size; i++)
		os << inn.getStateU(i) << " ";
	os << endl << endl;
	//Write weights
	for (int i = 1; i <= size; i++) {
		for (int j = 1; j <= size; j++)
			os << inn.getWeights(i,j) << " ";
		os << endl << endl;
	}
	//Write outputs
	for(int i=1; i<=size; i++){
		os << inn.getOutput(i) << " ";
	os << endl;
	}
	// Return the ostream
	return os;
}
istream& operator>>(istream& is, IzhiNN& inn){
	//Read size
	int ne,ni,size;
	is >> ne;
	is >> ni;
	inn.setNetworkSize(ne,ni);
	size = ne+ni;
	// Read params
	for (int i = 1; i <= size; i++){
		double paramAi;
		is >> paramAi;
		inn.setParamA(i,paramAi);
	}
	for (int i = 1; i <= size; i++){
		double paramBi;
		is >> paramBi;
		inn.setParamB(i,paramBi);
	}
	for (int i = 1; i <= size; i++){
		double paramCi;
		is >> paramCi;
		inn.setParamC(i,paramCi);
	}
	for (int i = 1; i <= size; i++){
		double paramDi;
		is >> paramDi;
		inn.setParamD(i,paramDi);
	}
	// Read states
	for (int i = 1; i <= size; i++){
		double stateVi;
		is >> stateVi;
		inn.setStateV(i,stateVi);
	}
	for (int i = 1; i <= size; i++){
		double stateUi;
		is >> stateUi;
		inn.setStateU(i,stateUi);
	}
	// Read weights
	for (int i=1; i<=size; i++){
		for(int j=1; j<=size; j++){
			double weightij;
			is >> weightij;
			inn.setWeights(i,j,weightij);
		}
	}
}*/