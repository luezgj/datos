// ************************************************************
// Header file for Environment functions
// Madhavun Candadai Vasu
// Jun 22, 2016 - created
//
// ************************************************************
# pragma once

#include <iostream>
#include <cmath>
#include "VectorMatrix.h"
#include "random.h"

// Constants for distance sensors
const double COTAN0 = 3.73205080756888; 	// cotan(pi/12)
const double COTAN1 = 5.67128181961771;	    // cotan(pi/18)
const double COTAN2 = 11.43005230276134;    // cotan(pi/36)

class Environment{
    private:
        double DIAMETER;
        double ENVWIDTH;
        double OBJVEL;
        double ENVHEIGHT;
        double SENSORDISTANCE;
    public:
        Environment(){};
        Environment(double dia, double envWidth, double objVel, double envHt, double sensDist);
        ~Environment(){};

        void getDistanceInputs(TVector<double>& distanceInputs, bool isCircle, double aX, double oX, double oY);
        double getConcentration(double gradientPeak, int agentPos);
};