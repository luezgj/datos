// ************************************************************
// Source file for Izhikevich Neural network Delays
// Luciano Giménez
// Jun 23, 2020 - created
//
// ************************************************************
# pragma once
#include "Delays.h"

//Constructor de una red aleatoria
Delays::Delays() {
	//Setea distribución de delays y setea la cantidad en delays[][][]    Solo de la capa oculta, los demás no tienen delays
	for (int i = 0; i < cnsSize; i++) {
		short ind = 0;
		if (i < Ne) {
			bool impar= (Ne%2 ==1) && (i == Ne/2);
			int delaysImpar = (M / D) / 2;
			for (int j = 0; j < D; j++) {
				if (impar) {
					if (j == 0 && delaysImpar == 0) delays_length[i][j] = 1;
					else delays_length[i][j] = delaysImpar;	// uniform distribution of exc. synaptic delays
					for (int k = 0; k < delays_length[i][j]; k++)
						delays[i][j][k] = ind++;
				}
				else{
					delays_length[i][j] = M / D;	// uniform distribution of exc. synaptic delays
					for (int k = 0; k < delays_length[i][j]; k++)
						delays[i][j][k] = ind++;
				}
			}
		}
		else {
			for (int j = 0; j < D; j++) delays_length[i][j] = 0;
			delays_length[i][0] = M;			// all inhibitory delays are 1 ms
			for (int k = 0; k < delays_length[i][0]; k++)
				delays[i][0][k] = ind++;
		}
	}
}
//destructor
Delays::~Delays() {
	
}