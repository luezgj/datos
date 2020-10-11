// ************************************************************
// Header file for Izhikevich Neural network Architecture
// Luciano Gimenez
// Nov 13, 2019 - created
//
// Other headers from Randall Beer
// ************************************************************

# pragma once

#include "VectorMatrix.h"
#include "random.h"
#include "params.h"
#include <iostream>
#include <math.h>
#include "Delays.h"
#include "TSearch.h"

class IzhiNNArch {
private:
	int		post[N][M];				// indeces of postsynaptic neurons
								// En las conexiones postsinapticas de las salidas (que no tienen xq no se conecta nada para adelante), guardo las conexiones presinapticas
	
	int		N_pre[N], I_pre[N][3 * M]; // presynaptic information
	short 	D_pre[N][3 * M];	

	Delays* delays;
	// los borro porque van a ocupar mucho espacio y no los estoy usando: float	a[N], d[N];				// neuronal dynamics parameters
public:
	//constructor
	IzhiNNArch(){};
	//destructor
	~IzhiNNArch();

	void initArch(TVector<int> gen, Delays * d);
	//return the symetric neuron number of neuronNumber
	static int sim(int neuronNumber);

	float getA(int neuron){return ((neuron<Ne)?AExc:AInh); /*a[neuron-1];*/};
	float getD(int neuron){return ((neuron<Ne)?DExc:DInh);/*d[neuron-1];*/};

	int getPost(int from, int conNumber){return post[from][conNumber];};
	int getPre(int from, int conNumber){return I_pre[from][conNumber];};

	int getNPre(int from){return N_pre[from];};
	int getDPre(int from, int conNumber){return D_pre[from][conNumber];};
};