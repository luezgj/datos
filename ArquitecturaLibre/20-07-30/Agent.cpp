// ************************************************************
// Agent class definition
// Madhavun Candadai Vasu
// Jun 06, 2016 - created
//
// March, 09, 2020- modified by Luciano Gimenez
//
// TODO: in initNS use setters to set private attributes
// ************************************************************

//************
// Constructor
//************
#include "Agent.h"
#define SYM

//************
// Initialize
//************
void Agent::initNS(genotype genotype,Delays * delays){
    //cout << "Initializing NS from Agent class" << endl;
    
    //cout << "\tcalling initCNS" << endl;
    ns.initCNS(genotype,delays);

    //cout << "\tDone initing Agent" << endl;
}

//************
// Control
//************
double Agent::step(double stepSize, TVector<double> distanceInputs){
    //cout << "\t\tIn Agent step" << endl;
    double motorNeuronsDiff = ns.step(stepSize, distanceInputs);
    double velocity = (motorNeuronsDiff*rateCodeGain)/       /*ver que onda esto (antes estaba bien pero ahora ya no se que onda, lo dejo con Ne+Ni)-->  lo cambié por outConnections para probar*/ outConnections;


    double offset = stepSize*velocity;
    //cout << "\t\tVelocity = " << velocity << " offset = " << offset << endl;
    return offset;
}

void Agent::reset(){
    ns.reset();
}