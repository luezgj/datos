// ************************************************************
// Source file for Izhikevich Neural network Architecture
// Luciano Giménez
// Nov 13, 2019 - created
//
// ************************************************************
# pragma once
#include "IzhiNNArch.h"

//return the number of the symetric neuron
int IzhiNNArch::sim(int neuronNumber) {
	int number;
	int translate;
	int tam;
	if (neuronNumber < 0) return -1;
	
	if (neuronNumber < Ne) {
		translate = 0;
		tam = Ne;
	}
	else if (neuronNumber < Ne + Ni) {
		translate = Ne;
		tam = Ni;
	}
	else if (neuronNumber < Ne + Ni + numDistSensors) {
		translate = Ne + Ni;
		tam = numDistSensors;
	}
	else {
		translate = Ne + Ni + numDistSensors;
		tam = numMotorNeurons;
	}
	neuronNumber = neuronNumber - translate;
	if (neuronNumber < floor(float(tam) / 2)) number = neuronNumber + ((floor(float(tam) / 2) - neuronNumber) * 2);
	else number = neuronNumber - ((neuronNumber - (floor(double(tam) / 2))) * 2);
	return ((tam % 2 == 0) ? number - 1 : number) + translate;
	

	/*
	if (neuronNumber < Ne/2) {
		translate = 0;
		tam = Ne/2;
	}
	else if (neuronNumber < Ne) {
		translate = Ne / 2;
		tam = Ne/2;
	}
	else if (neuronNumber < Ne + (Ni/2)) {
		translate = Ne;
		tam = Ni/2;
	}
	else if (neuronNumber < Ne + Ni) {
		translate = Ne + (Ni / 2);
		tam = Ni/2;
	}
	else if (neuronNumber < Ne + Ni + numDistSensors) {
		translate = Ne + Ni;
		tam = numDistSensors;
	}
	else {
		translate = Ne + Ni + numDistSensors;
		tam = numMotorNeurons;
	}
	neuronNumber = neuronNumber - translate;
	if (neuronNumber < floor(float(tam) / 2)) number = neuronNumber + ((floor(float(tam) / 2) - neuronNumber) * 2);
	else number = neuronNumber - ((neuronNumber - (floor(double(tam) / 2))) * 2);
	return ((tam % 2 == 0) ? number - 1 : number) + translate;
	*/
}

//Inicia una red con la arquitectura de gen
void IzhiNNArch::initArch(TVector<int> gen, Delays * d) {
	int i, j, k, jj, dd, postNeuron;
	
	delays=d;
	/* no modifico ni a ni d, por lo que los voy a manejar como constantes
	for (i = 0; i < Ne; i++) a[i] = 0.02;// RS type
	for (i = Ne; i < N; i++) a[i] = 0.1; // FS type

	for (i = 0; i < Ne; i++) d[i] = 8.0; // RS type
	for (i = Ne; i < N; i++) d[i] = 2.0; // FS type
	*/

	int index=1;
	//Conecta las neuronas excitatorias
	for (i = 0; i < (Ne/2); i++) for (j = 0; j < M; j++) {
		postNeuron= gen[index];
		post[i][j] = postNeuron;
		post[sim(i)][j] = sim(postNeuron);
		index++;
	}


	//Conecta las neuronas inhibitorias
	for (i = Ne; i < Ne+(Ni/2); i++) for (j = 0; j < M; j++) {
		postNeuron= gen[index];
		post[i][j] = postNeuron;
		post[sim(i)][j] = sim(postNeuron);
		index++;
	}

	//Conecta las neuronas sensoras
	for (i = Ne + Ni ; i < Ne+ Ni + (numDistSensors/2); i++) for (j = 0; j < inConnections; j++) {
		postNeuron= gen[index];
		post[i][j] = postNeuron;
		post[sim(i)][j] = sim(postNeuron);
		index++;
	}
	//sensora del medio
	if (numDistSensors%2==1) for (j = 0; j < inConnections; j++) {
		i = Ne + Ni + (numDistSensors / 2);
		postNeuron= gen[index];
		post[i][j] = postNeuron;
		post[i][++j] = sim(postNeuron);
		index++;
	}

	//Conectar interneuronas a las salidas,  solo conecta exhitatorias->salida
	//se guardan presinapsis en la matriz de post
	//solo hay dos motoras, no hace falta tratar la del medio
	for (i = Ne + Ni + numDistSensors; i < Ne + Ni + numDistSensors + ceil(float(numMotorNeurons)/2); i++) for (j = 0; j < outConnections; j++) {
		postNeuron= gen[index];
		post[i][j] = postNeuron;
		post[sim(i)][j] = sim(postNeuron);
		index++;
	}
	

	//Busca las conexiones presinapticas,   ignora las de las motoras xq ya son presinapticas
	for (i = 0; i < N-numMotorNeurons; i++) {
		N_pre[i] = 0;
		for (j = 0; j < Ne; j++)
			for (k = 0; k < M; k++){
				//cout << "La conexion [" << j <<","<<k <<"] que es "<<  post[j][k] << "es igual a la neurona"<< i <<"?" ;
				if (post[j][k] == i) {		// find all presynaptic neurons 
					I_pre[i][N_pre[i]] = j;	// add this neuron to the list
					for (dd = 0; dd < D; dd++)	// find the delay
						for (jj = 0; jj < delays->getDelayLength(j,dd); jj++)
							if (post[j][delays->getConnNumberByDelay(j,dd,jj)] == i) D_pre[i][N_pre[i]++] = dd;
				}
			}
	}

	//Copia las conexiones presinapticas de las neu. motoras a las otras matrices
	for (i = N - numMotorNeurons; i < N; i++) {
		N_pre[i] = outConnections;
		for (int m = 0; m < outConnections;m++) {
			I_pre[i][m] = post[i][m];	// add this neuron to the list
			D_pre[i][m] = 0;			//D_pre no tiene sentido porque no tiene delays
		}
	}
}
//destructor
IzhiNNArch::~IzhiNNArch() {
	
}