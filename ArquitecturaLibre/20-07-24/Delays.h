// ************************************************************
// Header file for Izhikevich Neural network Delays
// Luciano Gimenez
// Jun 23, 2020 - created
//
// Other headers from Randall Beer
// ************************************************************
#pragma once
#include "params.h"

class Delays {
private:
	//por ahora los delays son iguales para todas las redes
	short	delays_length[N][D];	// distribution of delays
	short	delays[N][D][M];		// arrangement of delays	

public:
	//constructor
	Delays();
	//destructor
	~Delays();

	int getDelayLength(int neuron, int delay){return delays_length[neuron][delay];};
	//int getDelay(int from, int conNumber){return I_pre[from-1][conNumber-1];}; //creo que está mal
	int getConnNumberByDelay(int from, int delay, int order/*numero de conexion con ese delay*/){return (delays[from][delay][order]);};
};