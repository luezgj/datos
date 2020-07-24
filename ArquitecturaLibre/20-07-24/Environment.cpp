// ************************************************************
// Environment class definition
// Madhavun Candadai Vasu
// Jun 22, 2016 - created
//
//
// ************************************************************
#include "Environment.h"

Environment::Environment(double dia, double envWidth, double objVel, double envHt, double sensDist){
	DIAMETER = dia;               // diameter of cirle and agent, length of horizontal line
    	ENVWIDTH = envWidth;
	OBJVEL = objVel;                  // rate at which object falls
 	ENVHEIGHT = envHt;
 	SENSORDISTANCE = sensDist;
//	cout << "Environment params - " << DIAMETER <<  " " << ENVWIDTH << " " <<  OBJVEL << " " << ENVHEIGHT << " " << SENSORDISTANCE << endl;
 
}

void Environment::getDistanceInputs(TVector<double>& distanceInputs, bool isCircle, double aX, double oX, double oY){
    double aY = 0;
    double numDistSensors = distanceInputs.Size();
//	cout << "From func - " << DIAMETER <<  " " << ENVWIDTH << " " <<  OBJVEL << " " << ENVHEIGHT << " " << SENSORDISTANCE << endl;
//	cout << "From func - " << oX << " " << oY << " " << aX << endl;
 
    //TVector<double> distanceInputs;
    //distanceInputs.SetBounds(1,numDistSensors);
    //distanceInputs.FillContents(0.);
    double m, b;
    for (int i=0;i<numDistSensors;i++){
        double input_value;
        // Find m, b for y = mx + b equation describing ray
        switch(i) {
            case 0:
                m = -COTAN0;
                break;
            case 1:
                m  = -COTAN1;
                break;
            case 2:
                m  = -COTAN2;
                break;
            case 3:
                m  = 0;
                break;
            case 4:
                m  = COTAN2;
                break;
            case 5:
                m  = COTAN1;
                break;
            case 6:
                m  = COTAN0;
                break;
        }

        b = aY - m*aX;

        // Line Segment
        if (!isCircle) {
            double x_intersect;
            if (m == 0) {
                x_intersect = aX;
            } else {
                x_intersect = (oY - b)/m;
            }
            // No intersection
            if ((fabs(oX-x_intersect) > DIAMETER/2.0) || (oY < aY)) {
                input_value = 0;
            } else {
                double distance = sqrt(pow(x_intersect-aX, 2.0) + pow(oY-aY, 2.0)) - DIAMETER/2.0;
                input_value = 10*(SENSORDISTANCE-distance)/SENSORDISTANCE;
                if (input_value < 0) {
                    input_value = 0;
                } else if (input_value > 10) {
                    input_value = 10;
                }
            }
        // Circle
        }
        else {
            double A, B, C, discriminant;
            if (m == 0) {
                A = 1;
                B = -2*oY;
                C = pow(oY,2) + pow((aX-oX),2) - pow(DIAMETER/2,2);
                discriminant = pow(B, 2.0) - 4*A*C;
                if (discriminant < 0) {
                    input_value = 0;
                } else {
                    double y_intersect = (-B-sqrt(pow(B,2)-4*A*C))/(2*A);
                    double distance = y_intersect - (DIAMETER/2.0);
                    input_value = 10*(SENSORDISTANCE-distance)/SENSORDISTANCE;
                    if (input_value < 0) {
                        input_value = 0;
                    } else if (input_value > 10) {
                        input_value = 10;
                    }
                }
            }
            else {
                // Values for quadratic equation, to solve for x_intersect
                A = pow(m, 2) + 1;
                B = 2*(m*b - oX - oY*m);
                C = pow(oX, 2.0) + pow(oY, 2.0) + pow(b, 2.0) - 2*oY*b - pow(DIAMETER/2.0, 2.0);
                discriminant = pow(B, 2.0) - 4*A*C;
                // No intersection
                if (discriminant < 0) {
                    input_value = 0;
                } else {
                    double x_intersect1 = (-B+sqrt(discriminant))/(2*A);
                    double y_intersect1 = m*x_intersect1 + b;
                    double x_intersect2 = (-B-sqrt(discriminant))/(2*A);
                    double y_intersect2 = m*x_intersect2 + b;
                    if ((y_intersect1 < y_intersect2) && (y_intersect1 > aY)) {
                        double distance = sqrt(pow(x_intersect1-aX, 2.0) + pow(y_intersect1-aY, 2.0)) - DIAMETER/2.0;
                        input_value = 10*(SENSORDISTANCE-distance)/SENSORDISTANCE;
                        if (input_value < 0) {
                            input_value = 0;
                        } else if (input_value > 10) {
                            input_value = 10;
                        }
                    } else if (y_intersect2 > aY) {
                        double distance = sqrt(pow(x_intersect2-aX, 2.0) + pow(y_intersect2-aY, 2.0)) - DIAMETER/2.0;
                        input_value = 10*(SENSORDISTANCE-distance)/SENSORDISTANCE;
                        if (input_value < 0) {
                            input_value = 0;
                        } else if (input_value > 10) {
                            input_value = 10;
                        }
                    } else {
                        input_value = 0; // in case point of intersection is below sensor
                    }
                }
            }
        }
//	cout << "input value " << input_value << "||";
        distanceInputs[i+1] = input_value;
    }
//	cout << endl;
    //return distanceInputs;
}

double Environment::getConcentration(double gradientPeak, int agentPos){
    return 0.;
    // equation for concetration distribution
    double stdDev = 5.;
    static const float inv_sqrt_2pi = 0.3989422804014327;
    float a = (agentPos - gradientPeak) / stdDev;

    return 100*(inv_sqrt_2pi / stdDev * std::exp(-0.5f * a * a));
}