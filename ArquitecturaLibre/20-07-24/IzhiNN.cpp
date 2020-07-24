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

//Constructor 
IzhiNN::IzhiNN() {
	setNetworkVectors();
	N_firings=0;
	t=0;
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
		S[N]=0.0f;
		I[N]=0.0f;
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
	//Busca las conexiones presinapticas
	for (int i = 0; i < N; i++) {
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

/*void IzhiNN::randomInitNetwork() {
	//Do nothing
}*/


//*********
//CONTROL
//*********
void IzhiNN::eulerStep(double stepSize, TVector<double> distanceInputs) {
	//cout << "\t\t\tIn euler step for IzhiNN " << size << endl;
	// check for neurons that may spike
	for (int i = 0; i < size; i++) {
		////cout << "\t\t\t\t| i = " << i << "|v[i] = " << v[i];
		if (v[i] > SPIKEV) {
			// spike
			outputs[i] = 1;
			v[i] = -65.0/*c[i]*/;
			u[i] += arch->getD(i);
			firings[N_firings][0] = t;
			firings[N_firings++][1] = i;
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


	int k = N_firings;
	while ((k>0) && (t - firings[--k][0] < D)) {   //por todos los spikes que pasaron hace tan poco que importan
		for (int j = 0; j < del->getDelayLength(firings[k][1],t - firings[k][0]); j++)   //Por todas las conexiones que tienen un delay igual a (t actual- t spike)
		{
			int i = arch->getPost(firings[k][1],del->getConnNumberByDelay(firings[k][1],t - firings[k][0],j));
			//cout << "Hubo spike de " << firings[k][1] << " en el tiempo " << firings[k][0] << endl;
			//cout << "En el tiempo  " << t << "activa a " << i << endl;
			
			//llegó el spike a la neurona postsinaptica
			int adition = s[firings[k][1]][del->getConnNumberByDelay(firings[k][1], t - firings[k][0], j)];

			I[i] += adition;
			//cout << "Le suma " << adition<< endl;

			/*STDP
			if (firings[k][1] < Ne) // this spike is before postsynaptic spikes
				sd[firings[k][1]][delays[firings[k][1]][t - firings[k][0]][j]] -= LTD[i];*/
		}
	}	

	//ACÁ HAY QUE ACTUALIZAR LAS ENTRADAS EXTERNAS
	for (int j = 0; j < numDistSensors; j++) {
		//Chequear que esto esté bien
		for (int mm = 0; mm < inConnections; mm++) {			//Actualización de las neuronas conectadas al sensor
			//cout << "quiero acceder a S sub " << arch->getPost(Ne + Ni + j, mm) - 1 << endl;
			//cout << "que vale" << S[arch->getPost(Ne + Ni + j, mm) - 1] << endl;
			//cout << "quiero acceder a s sub " << Ne + Ni + j - 1 <<","<< mm - 1 << endl;
			//cout << "que vale" << s[Ne + Ni + j - 1][mm - 1] << endl;
			S[arch->getPost(Ne + Ni + j, mm)] += s[Ne + Ni + j][mm] * (1.0 / (1.0 + exp(-(distanceInputs[j+1] + input_bias))));
		}
	}

	for (int i = 0; i < cnsSize; i++) {

		v[i] += stepSize * 0.5 * ((0.04 * v[i] + 5) * v[i] + 140 - u[i] + I[i] + S[i]); // for numerical stability
		v[i] += stepSize * 0.5 * ((0.04 * v[i] + 5) * v[i] + 140 - u[i] + I[i] + S[i]); // time step is 0.5 ms
		u[i] += stepSize * arch->getA(i) * (0.2 * v[i] - u[i]);

		//STDP
		//LTP[i][t + D + 1] = 0.95 * LTP[i][t + D];
		//LTD[i] *= 0.95;
	}

	t++;
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