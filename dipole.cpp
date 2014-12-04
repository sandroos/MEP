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
#include <map>

#include "constants.h"
#include "emfield.h"
#include "paramreader.h"
#include "mathopers.h"

using namespace std;

static double pos[3];         // Position of the dipole.
static double ori[3];         // Orientation of the dipole
static double B_pole;         // Magnitude of B at the poles
static double R_pole;         // Distance from the origin to the poles.
static double m;              // Magnetic moment of the dipole.
static double r_min;
static double r_max;

EMfield::EMfield() {
	map<string,string> params = ParamReader::getParameters("Dipole");
	map<string,string> needed;
	needed["x_ori"];
	needed["y_ori"];
	needed["z_ori"];
	needed["B_pole_reference(T)"];
	needed["R_pole(km)"];
	needed["offset_angle(deg)"];
	needed["offset(km)"];
	needed["minradius(hermean)"];
	needed["maxradius(hermean)"];
	
	// Read parameters
	bool success = true;
	for (map<string,string>::iterator i = needed.begin(); i != needed.end(); ++i) {
		if (params.find(i->first) == params.end()) {
			cerr << "ERROR(Dipole): Parameter '" << i->first << "' was not found!" << endl;
			success = false;
			continue;
		}
		i->second = params[i->first];
	}
	if (success == false) {exit(1);}
	ori[0] = atof(needed["x_ori"].c_str());
	ori[1] = atof(needed["y_ori"].c_str());
	ori[2] = atof(needed["z_ori"].c_str());
	B_pole = atof(needed["B_pole_reference(T)"].c_str());
	if (needed["R_pole(km)"] == "MERCURY") R_pole = MERCURY_EQUAT_RADIUS/1000.0;
	else R_pole = atof(needed["R_pole(km)"].c_str());
	double offset = atof(needed["offset(km)"].c_str());
	double offangle = atof(needed["offset_angle(deg)"].c_str());
	r_min = atof(needed["minradius(hermean)"].c_str())*MERCURY_EQUAT_RADIUS;
	r_max = atof(needed["maxradius(hermean)"].c_str())*MERCURY_EQUAT_RADIUS;
	
	// Scale input parameters.
	unitVector(ori);
	R_pole *= 1000.0;
	m = B_pole*R_pole*R_pole*R_pole/2.0;
	offset *= 1000.0;
	offangle = offangle*M_PI/180.0;
	
	pos[0] = offset*cos(offangle);
	pos[1] = offset*sin(offangle);
	pos[2] = 0.0;
	
	#ifdef DEBUG
	  cout << "DIPOLE:" << endl;
	  cout << "\t B_pole_reference = " << B_pole << endl;
	  cout << "\t R_pole           = " << R_pole << endl;
 	  cout << "\t m                = " << m << endl;
	  cout << "\t orientation      = " << ori[0] << '\t' << ori[1] << '\t' << ori[2] << endl;
	  cout << "\t offset angle     = " << offangle*180.0/M_PI << endl;
	  cout << "\t offset distance  = " << offset << endl;
	#endif
}

EMfield::~EMfield() { }

bool EMfield::values(const double& t,const double r[3],double E[3],double B[3]) {
	static double rdot[3];    // Vector from dipole to particle position.
	static double r_unit[3];  // Unit vector from dipole to particle position.
	static double R;          // Magnitude of rdot.
	static double R3;         // R to the power of 3.
	static double m_dot_r;    // Dot product between dipole orientation and r_unit. 
	
	rdot[0] = r[0] - pos[0];
	rdot[1] = r[1] - pos[1];
	rdot[2] = r[2] - pos[2];
	unitVector(rdot,r_unit);
	m_dot_r = dotProduct(ori,r_unit);
	R = vectorMagnitude(rdot);
	R3 = R*R*R;

	E[0] = 0.0;
	E[1] = 0.0;
	E[2] = 0.0;
	B[0] = m*(3*m_dot_r*r_unit[0] - ori[0])/R3;
	B[1] = m*(3*m_dot_r*r_unit[1] - ori[1])/R3;
	B[2] = m*(3*m_dot_r*r_unit[2] - ori[2])/R3;
	
	R = vectorMagnitude(r);
	if (R <= r_min || R >= r_max) return false;
	return true;
}

double EMfield::maxB() {return fabs(B_pole);}














