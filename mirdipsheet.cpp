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

// A mirrored dipole model for Mercury with a current sheet in the tail. The planetary dipole is, in this case, 
// assumed to lie in the center of the planet, i.e., the planetary dipole offset = 0.
// The planetary dipole can be tilted in the yz-plane (x points towards the Sun). 
// The mirror dipole is somewhere in the x-axis, and it has the same tilt as the planetary 
// dipole.

// The following parameters relate to the planetary dipole.
static double tilt;            /**< Tilt angle. */
static double B0_pole;         /**< Magnitude of B at the poles. */
static double R0_pole;         /**< Distance from the origin to the poles. */
static double m0;              /**< Magnetic moment of the dipole. */

// The following parameters relate to the mirror dipole.
static double mirror_dipole_position;
static double B1_pole;
static double R1_pole;
static double m1;

// The following parameters relate to the current sheet.
static double Bs;              /**< Magnitude of the sheet magnetic field. */
static double L;               /**< Width of the current sheet in x-direction. */
static double D;               /**< Half-thickness of the current sheet in z-direction, i.e., thickness is 2D. */
static double x_s;             /**< Position of the front boundary of the current sheet in x-direction. */
static double D2;              /**< D squared. */
static double xs_minus_L;      /**< x_s - L. */

static double maxpos[3];
static double r_min;
static double r_max;

static double sin_tilt;        /**< sin of tilt angle. */
static double cos_tilt;        /**< cos of tilt angle. */
static double isin_tilt;       /**< sin of -1*tilt angle (bad naming, but anyway). */
static double icos_tilt;       /**< cos of -1*tilt angle. */

// Rotate by an angle tilt in yz-plane. 
inline void rotate(double& y,double& z) {
	static double ytmp;
	static double ztmp;
	ytmp = y*cos_tilt - z*sin_tilt;
	ztmp = y*sin_tilt + z*cos_tilt;
	y = ytmp;
	z = ztmp;
}

// Do an inverse rotation in yz-plane by an angle tilt.
inline void irotate(double& y,double& z) {
	static double ytmp;
	static double ztmp;
	ytmp = y*icos_tilt - z*isin_tilt;
	ztmp = y*isin_tilt + z*icos_tilt;
	y = ytmp;
	z = ztmp;
}

EMfield::EMfield() {
	map<string,string> params = ParamReader::getParameters("MirDipSheet");
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
	needed["B_sheet(T)"];
	needed["L_sheet_x-width(hermean)"];
	needed["D_sheet_z-thickness(hermean)"];
	needed["x_sheet_front_boundary(hermean)"];
	
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
	mirror_dipole_position = atof(needed["mirrordipole_distance(hermean)"].c_str())*MERCURY_EQUAT_RADIUS;
	tilt = -1.0*M_PI/180.0*atof(needed["tiltangle(deg)"].c_str());
	
	r_min = atof(needed["minradius(hermean)"].c_str())*MERCURY_EQUAT_RADIUS;
	r_max = atof(needed["maxradius(hermean)"].c_str())*MERCURY_EQUAT_RADIUS;
	maxpos[0] = atof(needed["x_maxradius(hermean)"].c_str())*MERCURY_EQUAT_RADIUS;
	maxpos[1] = atof(needed["y_maxradius(hermean)"].c_str())*MERCURY_EQUAT_RADIUS;
	maxpos[2] = atof(needed["z_maxradius(hermean)"].c_str())*MERCURY_EQUAT_RADIUS;
	
	Bs = atof(needed["B_sheet(T)"].c_str());
	L = atof(needed["L_sheet_x-width(hermean)"].c_str())*MERCURY_EQUAT_RADIUS;
	D = atof(needed["D_sheet_z-thickness(hermean)"].c_str())*MERCURY_EQUAT_RADIUS;
	x_s = atof(needed["x_sheet_front_boundary(hermean)"].c_str())*MERCURY_EQUAT_RADIUS;
	
	// Scale input parameters.
	m0 = B0_pole*R0_pole*R0_pole*R0_pole/2.0;
	m1 = B1_pole*R1_pole*R1_pole*R1_pole/2.0;
	D2 = D*D;
	xs_minus_L = x_s - L;
	sin_tilt = sin(tilt);
	cos_tilt = cos(tilt);
	isin_tilt = sin(-1.0*tilt);
	icos_tilt = cos(-1.0*tilt);
	
	//ori[0] = 0.0;
	//ori[1] = sin(tiltangle);
	//ori[2] = -cos(tiltangle);

	//pos[0] = mirrordist;
	//pos[1] = 0.0;
	//pos[2] = 0.0;
	
	#ifdef DEBUG
	  cerr << "DIPOLE:" << endl;
	  cerr << "\t B0_pole_reference = " << B0_pole << endl;
	  cerr << "\t R0_pole           = " << R0_pole << endl;
 	  cerr << "\t m0                = " << m0 << endl;
	  cerr << "\t tilt angle (rad)  = " << tilt << endl;
	  cerr << "\t B1_pole_reference = " << B1_pole << endl;
	  cerr << "\t R1_pole           = " << R1_pole << endl;
	  cerr << "\t m1                = " << m1 << endl;
	  cerr << "\t Mirror dipole pos = " << mirror_dipole_position << endl;
	  cerr << "\t Sheet field Bs    = " << Bs << endl;
	  cerr << "\t Sheet front bound = " << x_s << endl;
	  cerr << "\t Sheet x-width     = " << L << endl;
	  cerr << "\t Sheet z-thickness = " << D << endl;
	#endif
}

EMfield::~EMfield() { }

bool EMfield::values(const double& t,const double r[3],double E[3],double B[3]) {
	static double r_rot[3];   // r rotated into the tilted coordinate system.
	static double rdot[3];    // Vector from dipole to particle position.
	//static double r_unit[3];  // Unit vector from dipole to particle position.
	static double R;          // Magnitude of rdot.
	static double R2;         // rdot squared;
	static double R3;         // rdot to the power of 3.
	static double z2D2;       // Square root of z^2 + D^2.
	static double XsLx;       // x_s - L - x.
	static double Xsx;        // x_s - x;
	
	static bool rvalue;
	rvalue = true;

	r_rot[0] = r[0];
	r_rot[1] = r[1];
	r_rot[2] = r[2];
	rotate(r_rot[1],r_rot[2]);

	//cerr << r[0] << ' ' << r[1] << ' ' << r[2] << ' ';
	//cerr << r_rot[0] << ' ' << r_rot[1] << ' ' << r_rot[2] << endl;
	
	// First calculate the field due to the planetary dipole.
	rdot[0] = r_rot[0];
	rdot[1] = r_rot[1];
	rdot[2] = r_rot[2];
	R = vectorMagnitude(rdot);
	R2 = R*R;
	R3 = R*R2;
	B[0] = 3*m0*rdot[0]*rdot[2]/R2/R3;
	B[1] = 3*m0*rdot[1]*rdot[2]/R2/R3;
	B[2] = m0/R3*(3*rdot[2]*rdot[2]/R2 - 1.0);	
	if (R <= r_min) rvalue = false;

	// Then add the contribution from the tail field.
	
	z2D2 = sqrt(rdot[2]*rdot[2] + D2);
	XsLx = xs_minus_L - rdot[0];
	Xsx  = x_s - rdot[0];
	B[0] += rdot[2]*Bs/z2D2*(atan(XsLx/z2D2) - atan(Xsx/z2D2));
	B[2] += 0.5*Bs*log((XsLx*XsLx+z2D2*z2D2)/(Xsx*Xsx+z2D2*z2D2));
	
	// Finally, add the contribution from the mirror dipole.
	rdot[0] = r_rot[0] - mirror_dipole_position;
	rdot[1] = r_rot[1];
	rdot[2] = r_rot[2];
	R = vectorMagnitude(rdot);
	R2 = R*R;
	R3 = R*R2;
	B[0] += 3*m1*rdot[0]*rdot[2]/R2/R3;
	B[1] += 3*m1*rdot[1]*rdot[2]/R2/R3;
	B[2] += m1/R3*(3*rdot[2]*rdot[2]/R2 - 1.0);
	
	// Check that the particle is still inside the "exit sphere". 
	rdot[0] = r_rot[0] - maxpos[0];
	rdot[1] = r_rot[1] - maxpos[1];
	rdot[2] = r_rot[2] - maxpos[2];
	if (vectorMagnitude(rdot) >= r_max) rvalue = false;
	
	// Rotate the magnetic field back into original coordinate system.
	irotate(B[1],B[2]);
	
	return rvalue;
}

double EMfield::maxB() {return fabs(B0_pole);}














