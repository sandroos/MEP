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

#include "constants.h"

using namespace std;

int main(int argn,char* args[]) {
	
	if (argn != 2) {
		cerr << endl;
		cerr << "USAGE: ./writesphere <radius>" << endl;
		cerr << "\t Where radius is in hermean radii, e.g., 1.0" << endl;
		cerr << endl;
		return 1;
	}
	double rad = atof(args[1]);
	
	unsigned int N_theta = 7;
	unsigned int N_phi = 10;
	double theta_min = 5.0;
	double theta_max = 175.0;
	
	theta_min = theta_min*M_PI/180.0;
	theta_max = theta_max*M_PI/180.0;
	
	double d_theta = (theta_max-theta_min)/(N_theta-1);
	double d_phi = 2.0*M_PI/(N_phi-1);

	double theta,phi;
	cout.precision(2);
	cout << fixed;
	
	for (unsigned int i=0; i<N_theta; ++i) {
		theta = theta_min + i*d_theta;
		for (unsigned int j=0; j<N_phi; ++j) {
			phi = j*d_phi;
			cout << rad*sin(theta)*cos(phi) << '\t' << rad*sin(theta)*sin(phi) << '\t' << rad*cos(theta) << endl;
		}
		cout << endl;
	}
	/*
	for (unsigned int j=0; j<N_phi; ++j) {
		phi = j*d_phi;
		for (unsigned int i=0; i<N_theta; ++i) {
			theta = theta_min + i*d_theta;
			cout << rad*sin(theta)*cos(phi) << '\t' << rad*sin(theta)*sin(phi) << '\t' << rad*cos(theta) << endl;
		}
		cout << endl;
	}
	*/
	return 0;
}

