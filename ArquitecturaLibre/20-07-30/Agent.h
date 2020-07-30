// ************************************************************
// Header file for Agent
// Madhavun Candadai Vasu
// Jun 06, 2016 - created
//
// ************************************************************

#include "NervousSystem.h"
#include "TSearch.h"

class Agent{
    private:
        double rateCodeGain;
        NervousSystem ns; // nervous system
    public:
        Agent(){rateCodeGain = 750;};
        ~Agent(){};
        
        void initNS(genotype gen,Delays * delays);
        double step(double stepSize, TVector<double> distanceInputs); // one step of the agent's NS and body (same time scale?) and return offset in position
        void reset();
        double getNeuronStateV(int n){return ns.getStateV(n);};
		double getInput(int n) { return ns.getInput(n); };
		double getSensoryInput(int n) { return ns.getSensoryInput(n); };
        void setRateCodeGain(double rg){rateCodeGain = rg;};
        void setWindowSize(double ws){ns.setWindowSize(ws);};
		void setInputBias(double inputBias) { ns.setInputBias(inputBias); };
		void setOutputBias(double bias) { ns.setOutputBias(bias); };

        void saveSpikes(const char* filename) { ns.saveSpikes(filename); };

};