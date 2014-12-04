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

// A mirrored dipole model for Mercury. The planetary dipole is, in this case, 
// assumed to lie in the center of the planet, i.e., the planetary dipole offset = 0.
// The planetary dipole can be tilted in the yz-plane (x points towards the Sun). 
// The mirror dipole is somewhere in the x-axis, and it has the same tilt as the planetary 
// dipole.

// The following parameters relate to the planetary dipole.
//static double pos0[3];         /**< Position of the dipole. */
static double ori[3];         /**< Orientation of the dipole. */
static double B0_pole;         /**< Magnitude of B at the poles. */
static double R0_pole;         /**< Distance from the origin to the poles. */
static double m0;              /**< Magnetic moment of the dipole. */

// The following parameters relate to the mirror dipole.
static double pos[3];
static double B1_pole;
static double R1_pole;
static double m1;

static double maxpos[3];
static double r_min;
static double r_max;

EMfield::EMfield() {
	map<string,string> params = ParamReader::getParameters("MirrorDipole");
	map<string,string> needed;
	needed["mirrordipole_distance(hermean)"];
	needed["tiltangle(deg)"];
	needed["B0_pole_reference(T)"];
	needed["R0_pole(hermean)"];
	needed["B1_pole_reference(T)"];
	needed["R1_pole(hermean)"];
	needed["minradius(hermean)"];
	needed["maxradius(hermean)"];
	needed["x_maxradius(hermean)"];
	needed["y_maxradius(hermean)"];
	needed["z_maxradius(hermean)"];
	
	// Read parameters, and check that all the needed parameters were given.
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
	
	// Convert the parameter values (strings) into floating points and scale values.
	B0_pole = atof(needed["B0_pole_reference(T)"].c_str());
	R0_pole = atof(needed["R0_pole(hermean)"].c_str())*MERCURY_EQUAT_RADIUS;
	B1_pole = atof(needed["B1_pole_reference(T)"].c_str());
	R1_pole = atof(needed["R1_pole(hermean)"].c_str())*MERCURY_EQUAT_RADIUS;
	double mirrordist = atof(needed["mirrordipole_distance(hermean)"].c_str())*MERCURY_EQUAT_RADIUS;
	double tiltangle = atof(needed["tiltangle(deg)"].c_str());
	
	r_min = atof(needed["minradius(hermean)"].c_str())*MERCURY_EQUAT_RADIUS;
	r_max = atof(needed["maxradius(hermean)"].c_str())*MERCURY_EQUAT_RADIUS;
	maxpos[0] = atof(needed["x_maxradius(hermean)"].c_str())*MERCURY_EQUAT_RADIUS;
	maxpos[1] = atof(needed["y_maxradius(hermean)"].c_str())*MERCURY_EQUAT_RADIUS;
	maxpos[2] = atof(needed["z_maxradius(hermean)"].c_str())*MERCURY_EQUAT_RADIUS;
	
	// Scale input parameters.
	m0 = B0_pole*R0_pole*R0_pole*R0_pole/2.0;
	m1 = B1_pole*R1_pole*R1_pole*R1_pole/2.0;
	tiltangle = M_PI/180.0*tiltangle;

	ori[0] = 0.0;
	ori[1] = sin(tiltangle);
	ori[2] = -cos(tiltangle);

	pos[0] = mirrordist;
	pos[1] = 0.0;
	pos[2] = 0.0;
	
	#ifdef DEBUG
	  cout << "DIPOLE:" << endl;
	  cout << "\t B0_pole_reference = " << B0_pole << endl;
	  cout << "\t R0_pole           = " << R0_pole << endl;
 	  cout << "\t m0                = " << m0 << endl;
	  cout << "\t tilt angle (rad)  = " << tiltangle << endl;
	  cout << "\t orientation       = " << ori[0] << '\t' << ori[1] << '\t' << ori[2] << endl;
	  cout << "\t B1_pole_reference = " << B1_pole << endl;
	  cout << "\t R1_pole           = " << R1_pole << endl;
	  cout << "\t m1                = " << m1 << endl;
	  cout << "\t Mirror dipole pos = " << pos[0] << '\t' << pos[1] << '\t' << pos[2] << endl;
	#endif
}

EMfield::~EMfield() { }

bool EMfield::values(const double& t,const double r[3],double E[3],double B[3]) {
	static double rdot[3];    // Vector from dipole to particle position.
	static double r_unit[3];  // Unit vector from dipole to particle position.
	static double R;          // Magnitude of rdot.
	static double R3;         // R to the power of 3.
	static double m_dot_r;    // Dot product between dipole orientation and r_unit. 
	bool rvalue = true;
	
	// First calculate the field due to the planetary dipole.
	unitVector(r,r_unit);
	m_dot_r = dotProduct(ori,r_unit);
	R = vectorMagnitude(r);
	R3 = R*R*R;
	B[0] = m0*(3*m_dot_r*r_unit[0])/R3;
	B[1] = m0*(3*m_dot_r*r_unit[1] - ori[1])/R3;
	B[2] = m0*(3*m_dot_r*r_unit[2] - ori[2])/R3;
	if (R <= r_min) rvalue = false;
	
	// Then add the contribution from the mirror dipole.
	rdot[0] = r[0] - pos[0];
	rdot[1] = r[1];
	rdot[2] = r[2];
	unitVector(rdot,r_unit);
	m_dot_r = dotProduct(ori,r_unit);
	R = vectorMagnitude(rdot);
	R3 = R*R*R;
	B[0] += m1*(3*m_dot_r*r_unit[0])/R3;
	B[1] += m1*(3*m_dot_r*r_unit[1] - ori[1])/R3;
	B[2] += m1*(3*m_dot_r*r_unit[2] - ori[2])/R3;	

	// Check that the particle is still inside the "exit sphere". 
	rdot[0] = r[0] - maxpos[0];
	rdot[1] = r[1] - maxpos[1];
	rdot[2] = r[2] - maxpos[2];
	if (vectorMagnitude(rdot) >= r_max) rvalue = false;	
	return rvalue;
}

double EMfield::maxB() {return fabs(B0_pole);}














