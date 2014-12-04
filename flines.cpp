// This file is part of MEP.
//
//    Copyright 2014 Arto Sandroos.
//
//    MEP is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    MEP is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with MEP.  If not, see <http://www.gnu.org/licenses/>.

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>

#include "constants.h"
#include "mathopers.h"
#include "paramreader.h"
#include "emfield.h"

using namespace std;

struct Fline {
	vector<double> x;
	vector<double> y;
	vector<double> z;
};

/** A very barbaric field line calculator. First it determines the field direction at the given 
 * start point, and then starts to calculate the magnetic field line using an Eulerian iterator.
 * There is a magic number SAVEINTERVAL here, that denotes the interval the field line coordinates 
 * are written out.
 * Modification by RV on 25.2.09: we first take a half-step to the mid point of the step, then 
 * evaluate the field and then take the full step.
 */
void calcline(EMfield& field,double x0,double y0,double z0,vector<double>& x,vector<double>& y,vector<double>& z) {
	double t = 0.0;
	double pos[3];
	double E[3];
	double B[3];
	double Bunit[3];
	pos[0] = x0;
	pos[1] = y0;
	pos[2] = z0;
	field.values(t,pos,E,B);
	
	x.clear();
	y.clear();
	z.clear();
	// Determine whether we should propagate to the parallel or antiparallel field direction,
	// so that we do not immediately encounter the surface of the planet. 
	double dir = 1.0;
	if (dotProduct(B,pos) < 0.0) dir = -1.0;
	const double STEP = 1000.0; // The step size, i.e., 1 km.
	
	unsigned int counter = 0;
	unsigned int SAVEINTERVAL = 50;
	do {
		if (counter % SAVEINTERVAL == 0) {
			x.push_back(pos[0] / MERCURY_EQUAT_RADIUS);
			y.push_back(pos[1] / MERCURY_EQUAT_RADIUS);
			z.push_back(pos[2] / MERCURY_EQUAT_RADIUS);
			counter = 0;
		}
		++counter;
			
		unitVector(B,Bunit);
		pos[0] += dir*STEP/2*Bunit[0];
		pos[1] += dir*STEP/2*Bunit[1];
		pos[2] += dir*STEP/2*Bunit[2];		
		if (field.values(t,pos,E,B)) {
		  pos[0] += -dir*STEP/2*Bunit[0];
		  pos[1] += -dir*STEP/2*Bunit[1];
		  pos[2] += -dir*STEP/2*Bunit[2];		
		  unitVector(B,Bunit);
		  pos[0] += dir*STEP*Bunit[0];
		  pos[1] += dir*STEP*Bunit[1];
		  pos[2] += dir*STEP*Bunit[2];		
		}
	} while (field.values(t,pos,E,B) == true);
}

int main(int argn,char* args[]) {
	if (argn != 5) {
		cerr << endl;
		cerr << "USAGE: ./flines <tilt angle(deg)> <N_theta> <N_phi> <parameter file>" << endl;
		cerr << endl;
		return 1;
	}
	
	double tilt = atof(args[1])*M_PI/180.0;
	int N_theta = atoi(args[2]);
	int N_phi = atoi(args[3]);
	ParamReader paramReader(args[4]);
	EMfield emfield;
	
	if (N_theta < 1) {
		cerr << "ERROR: N_theta must be greater than 0." << endl;
		return 1;
	}
	if (N_phi < 1) {
		cerr << "ERROR: N_phi must be greater than 0." << endl;
		return 1;
	}
	
	double theta,phi;
	double theta_min,theta_max,d_theta;
	double phi_min,phi_max,d_phi;
	// Determine appropriate min/max value for theta and phi angles. 
	theta_min = M_PI/2.0/(N_theta+1);
	theta_max = M_PI/2.0 - theta_min;
	if (N_theta == 1) d_theta = 0.0;
	else d_theta = (theta_max-theta_min) / (N_theta - 1);
	
	phi_min = 0.0;
	phi_max = 2.0*M_PI - 2.0*M_PI/N_phi;
	if (N_phi == 1) d_phi = 0.0;
	else d_phi = (phi_max-phi_min) / (N_phi -1);

	// Calculate starting points for field line calculation, using a tilted 
	// spherical coordinate system. 
	vector<double> xcrds;
	vector<double> ycrds;
	vector<double> zcrds;
	vector<Fline> fieldLines;
	const double R = 1.01*MERCURY_EQUAT_RADIUS;
	for (int j=0; j<N_phi; ++j) {
		phi = phi_min + j*d_phi;
		for (int i=0; i<N_theta; ++i) {
			theta = theta_min + i*d_theta;
			double x = R*sin(theta)*cos(phi);
			double y = R*sin(theta)*sin(phi);
			double z = R*cos(theta);
			
			xcrds.push_back( x );
			ycrds.push_back( y*cos(tilt) - z*sin(tilt) );
			zcrds.push_back( y*sin(tilt) + z*cos(tilt) );
			unsigned int index = xcrds.size()-1;
			Fline f;
			calcline(emfield,xcrds[index],ycrds[index],zcrds[index],f.x,f.y,f.z);
			fieldLines.push_back(f);
		}
	}
	// Write out the calculated field lines.
	cout.precision(3);
	cout << fixed;
	for (unsigned int i=0; i<fieldLines.size(); ++i) {
		for (unsigned int j=0; j<fieldLines[i].x.size(); ++j) 
			cout << fieldLines[i].x[j] << '\t' << fieldLines[i].y[j] << '\t' << fieldLines[i].z[j] << endl;
		cout << endl;
	}
	
	return 0;
}














