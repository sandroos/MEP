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
#include <fstream>

#include "constants.h"
#include "injector.h"
#include "paramreader.h"
#include "mathopers.h"
#include "randomnumbers.h"

using namespace std;

static const double GAMMALIMIT = 1.0001;

Injector::Injector() {
	N_injected = 0;
	
	// Define the parameters that need to be read.
	map<string,string> params = ParamReader::getParameters("Injector");
	map<string,string> needed;
	needed["q(elementary)"] = "";
	needed["m(protons)"] = "";
	needed["pwrlawindex"] = "";
	needed["energy_min(keV)"] = "";
	needed["energy_max(keV)"] = "";
	needed["particles(int)"] = "";
	needed["x_spacecraft(km)"] = "";
	needed["y_spacecraft(km)"] = "";
	needed["z_spacecraft(km)"] = "";
	needed["x_velocity"] = "";
	needed["y_velocity"] = "";
	needed["z_velocity"] = "";
	needed["conehalfwidth(deg)"] = "";
	needed["instrument"] = "";
	
	// Read the values of the needed parameters.
	bool success = true;
	for (map<string,string>::iterator i = needed.begin(); i != needed.end(); ++i) {
		if (params.find(i->first) == params.end()) {
			cerr << "ERROR(Injector): Parameter '" << i->first << "' was not found!" << endl;
			success = false;
			continue;
		}
		i->second = params[i->first];
	}
	if (success == false) {exit(1);}
	q_inj = atof(needed["q(elementary)"].c_str());
	m_inj = atof(needed["m(protons)"].c_str());
	r_inj[0] = atof(needed["x_spacecraft(km)"].c_str());
	r_inj[1] = atof(needed["y_spacecraft(km)"].c_str());
	r_inj[2] = atof(needed["z_spacecraft(km)"].c_str());
	unitvel[0] = atof(needed["x_velocity"].c_str());
	unitvel[1] = atof(needed["y_velocity"].c_str());
	unitvel[2] = atof(needed["z_velocity"].c_str());
	double cone = atof(needed["conehalfwidth(deg)"].c_str());
	pwrlawindex = atof(needed["pwrlawindex"].c_str());
	W_inj_min = atof(needed["energy_min(keV)"].c_str());
	W_inj_max = atof(needed["energy_max(keV)"].c_str());
	N_particles = atoi(needed["particles(int)"].c_str());
	instr = atoi(needed["instrument"].c_str())-1;
	
	// Scale variables.
	unitVector(unitvel);
	if (fabs(pwrlawindex) < 1.0e-6) pwrlawindex = 0.0;
	W_inj_min = W_inj_min*1000.0*ELEMENTARY_CHARGE;
	W_inj_max = W_inj_max*1000.0*ELEMENTARY_CHARGE;
	q_inj = q_inj * ELEMENTARY_CHARGE;
	m_inj = m_inj * PROTON_MASS;
	for (int i=0; i<3; ++i) r_inj[i] *= 1000.0;
	
	y_min = 1.0;
	y_max = W_inj_max / W_inj_min;
	norm_energy = (pwrlawindex + 1.0) / (pow(y_max,1.0+pwrlawindex) - pow(y_min,1.0+pwrlawindex));
	
	// Create an orthonormal basis where the z-axis points to the direction of the nadir,
	// x-axis to the direction of the spacecraft velocity, and y-axis completes the set. 
	unitnadir[0] = -r_inj[0];
	unitnadir[1] = 0.0;
	unitnadir[2] = -1.0*unitnadir[0]*unitvel[0]/unitvel[2];
	unitVector(unitnadir);
	crossProduct(unitnadir,unitvel,unity);

	instruments.resize(5);
	instruments[0].e[0] = -1.0*unitnadir[0]; // Instrument #1
	instruments[0].e[1] = -1.0*unitnadir[1];
	instruments[0].e[2] = -1.0*unitnadir[2];
	instruments[1].e[0] = cos(135*M_PI/180.0)*unitvel[0] + sin(135*M_PI/180.0)*unity[0]; // Instrument #2
	instruments[1].e[1] = cos(135*M_PI/180.0)*unitvel[1] + sin(135*M_PI/180.0)*unity[1];
	instruments[1].e[2] = cos(135*M_PI/180.0)*unitvel[2] + sin(135*M_PI/180.0)*unity[2];
	instruments[2].e[0] = cos(45*M_PI/180.0)*unitvel[0] + sin(45*M_PI/180.0)*unity[0]; // Instrument #3
	instruments[2].e[1] = cos(45*M_PI/180.0)*unitvel[1] + sin(45*M_PI/180.0)*unity[1];
	instruments[2].e[2] = cos(45*M_PI/180.0)*unitvel[2] + sin(45*M_PI/180.0)*unity[2];
	instruments[3].e[0] = cos(315*M_PI/180.0)*unitvel[0] + sin(315*M_PI/180.0)*unity[0]; // Instrument #4
	instruments[3].e[1] = cos(315*M_PI/180.0)*unitvel[1] + sin(315*M_PI/180.0)*unity[1];
	instruments[3].e[2] = cos(315*M_PI/180.0)*unitvel[2] + sin(315*M_PI/180.0)*unity[2];
	instruments[4].e[0] = cos(225*M_PI/180.0)*unitvel[0] + sin(225*M_PI/180.0)*unity[0]; // Instrument #5
	instruments[4].e[1] = cos(225*M_PI/180.0)*unitvel[1] + sin(225*M_PI/180.0)*unity[1];
	instruments[4].e[2] = cos(225*M_PI/180.0)*unitvel[2] + sin(225*M_PI/180.0)*unity[2];
	for (unsigned int i=0; i<instruments.size(); ++i) {
		instruments[i].coneWidth = M_PI*cone/180.0;
		instruments[i].cos_cone = cos(M_PI*cone/180.0);
		for (int j=0; j<3; ++j) instruments[i].e[j] = -1*instruments[i].e[j];
	}

	if (instr < 0 || instr >= static_cast<int>(instruments.size())) {
		cerr << "ERROR(Injector): Chosen instrument #" << instr << " does not exists!" << endl;
		exit(1);
	}
	writeGnuplot();
	#ifdef DEBUG
	  cout << "INJECTOR:" << endl;
	  cout << "\t N_particles     = " << N_particles << endl;
	  cout << "\t q               = " << q_inj << endl;
	  cout << "\t m               = " << m_inj << endl;
	  cout << "\t W(min/max)      = " << W_inj_min << "\t" << W_inj_max << endl;
	  cout << "\t pwrlawindex     = " << pwrlawindex << endl;
	  cout << "\t Spacecraft pos  = " << r_inj[0] << '\t' << r_inj[1] << '\t' << r_inj[2] << endl;
	  cout << "\t Nadir direction = " << unitnadir[0] << '\t' << unitnadir[1] << '\t' << unitnadir[2] << endl;
	  cout << "\t Velocity direc  = " << unitvel[0] << '\t' << unitvel[1] << '\t' << unitvel[2] << endl;
	  for (unsigned int i=0; i<instruments.size(); ++i) {
		  cout << "\t Instrument #" << i << "  = ";
		  for (int j=0; j<3; ++j) cout << instruments[i].e[j] << '\t';
		  cout << endl;
	  }
	#endif
}

Injector::~Injector() { }

void Injector::writeGnuplot() {
	fstream out("sat.txt", fstream::out);
	out.precision(2);
	out << fixed;
	// Write the nadir vector
	out << r_inj[0]/MERCURY_EQUAT_RADIUS << '\t' << r_inj[1]/MERCURY_EQUAT_RADIUS << '\t';
	out << r_inj[2]/MERCURY_EQUAT_RADIUS << '\t';
	for (int i=0; i<3; ++i) out << 0.3*unitnadir[i] << '\t';
	out << endl;
	// Write the velocity vector
	out << r_inj[0]/MERCURY_EQUAT_RADIUS << '\t' << r_inj[1]/MERCURY_EQUAT_RADIUS << '\t';
	out << r_inj[2]/MERCURY_EQUAT_RADIUS << '\t';
	for (int i=0; i<3; ++i) out << 0.3*unitvel[i] << '\t';
	out << endl;
	// Write instrument orientations
	for (unsigned int j=0; j<instruments.size(); ++j) {
		out << r_inj[0]/MERCURY_EQUAT_RADIUS << '\t' << r_inj[1]/MERCURY_EQUAT_RADIUS << '\t';
		out << r_inj[2]/MERCURY_EQUAT_RADIUS << '\t';
		for (int i=0; i<3; ++i) out << -0.2*instruments[j].e[i] << '\t';
		out << endl;
	}
	out.close();
	
	out.open("gp_sat", fstream::out);
	for (unsigned int i=0; i<instruments.size(); ++i) {
		out << "set label " << i+1 << " \""<<i+1<<"\" at ";
		for (int j=0; j<3; ++j) {
			out << r_inj[j]/MERCURY_EQUAT_RADIUS - 0.25*instruments[i].e[j];
			if (j != 2) out << ',';
		}
		out << endl;
	}
	out.close();
}

bool Injector::getNewParticle(double& t,double& q,double& m,double& g,double r[3],double v[3]) {
	static double energy;
	static double speed;
	if (N_injected == N_particles) return false;
	
	// Calculate the (random) injection energy
	if (pwrlawindex == 0.0) {
		energy = W_inj_min + RandomNumbers::uniform()*(W_inj_max-W_inj_min);
	} else {
		energy = pow(y_min,1.0+pwrlawindex) + (pwrlawindex+1.0)/norm_energy * RandomNumbers::uniform();
		energy = pow(energy,1.0/(1.0+pwrlawindex));
		energy = energy * W_inj_min;
	}
	
	// Calculate the gamma-factor and particle speed.
	g = 1.0 + energy/m_inj/SPEED_OF_LIGHT/SPEED_OF_LIGHT;
	if (g > GAMMALIMIT) {
		speed = SPEED_OF_LIGHT*sqrt(g*g-1.0)/g;
	} else {
		speed = sqrt(2.0*energy/m_inj);
	}

	// Calculate a random velocity vector until we get one that is inside
	// the instruments detection cone.
	bool accept = false;
	do {
		double cos_theta = 1.0 - 2.0*RandomNumbers::uniform();
		double theta = acos(cos_theta);
		double phi = 2.0*M_PI*RandomNumbers::uniform();
		v[0] = speed*sin(theta)*cos(phi);
		v[1] = speed*sin(theta)*sin(phi);
		v[2] = speed*cos(theta);
		double cos_angle = dotProduct(v,instruments[instr].e) / speed;
		if (cos_angle >= instruments[instr].cos_cone) accept = true;
	} while (accept == false);
	t = 0.0;
	q = q_inj;
	m = m_inj;
	r[0] = r_inj[0];
	r[1] = r_inj[1];
	r[2] = r_inj[2];
	++N_injected;
	//cerr << "Injected" << endl;
	/*
	for (int i=0; i<3; ++i) cout << instruments[instr].e[i] << ' ';
	for (int i=0; i<3; ++i) cout << v[i] << ' ';
	cout << acos(dotProduct(v,instruments[instr].e)/speed)*180.0/M_PI << ' ';
	cout << acos(instruments[instr].cos_cone)*180.0/M_PI;
	cout << endl;
	 */
	return true;
}



