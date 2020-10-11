// ************************************************************
// Header file for Agent
// Madhavun Candadai Vasu
// Jun 06, 2016 - created
//
// ************************************************************

#include "NervousSystem.h"
#include "TSearch.h"
#include "params.h"

class Agent{
    private:
        double rateCodeGain;
        NervousSystem ns; // nervous system
    public:
        Agent(){rateCodeGain = 750;};
        ~Agent(){};
        
        void initNS(genotype gen,Delays * delays);
        double step(double stepSize, double (*distanceInputs)[numDistSensors]); // one step of the agent's NS and body (same time scale?) and return offset in position
        void reset();
        double getNeuronStateV(int n){return ns.getStateV(n);};
		double getInput(int n) { return ns.getInput(n); };
		double getSensoryInput(int n) { return ns.getSensoryInput(n); };
        void setRateCodeGain(double rg){rateCodeGain = rg;};
        void setWindowSize(double ws){ns.setWindowSize(ws);};
		void setInputBias(double inputBias) { ns.setInputBias(inputBias); };
		void setOutputBias(double bias) { ns.setOutputBias(bias); };

        void saveSpikes(const char* filename) { ns.saveSpikes(filename); };

		int getPost(int from, int conNumber) { return ns.getPost(from,conNumber); };
        int getPre(int from, int conNumber){return ns.getPre(from,conNumber);};
        int getNPre(int from){return ns.getNPre(from);};
        int getDPre(int from, int conNumber){return ns.getDPre(from,conNumber);};


        double getWeight(int from, int nConection) { return ns.getWeight(from,nConection); };

		double getA(int nNumber) { return ns.getA(nNumber); };
		double getD(int nNumber) { return ns.getD(nNumber); };

        float getPreWeight(int neuron, int connection){return ns.getPreWeight(neuron,connection);};

		void changeSTDPStatus(bool status) { ns.changeSTDPStatus(status); };
};