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
#include "mover.h"
#include "emfield.h"

using namespace std;

static double ddt;

Mover::Mover() {
	r = NULL;
	v = NULL;
	E = NULL;
	B = NULL;
}

Mover::~Mover() { }

void Mover::initialize(const bool& backw,double& dt,const double& q,const double& m,double* const r,double* const v,double* const E,double* const B) {
	static const double C2 = SPEED_OF_LIGHT*SPEED_OF_LIGHT;
	this->r = r;
	this->v = v;
	this->E = E;
	this->B = B;
	// Calculate a value for the time step.
	gamma = sqrt(1.0 + (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]) / C2);
	dt = EMfield::maxB();
	dt = m/fabs(q)/dt/10.0;
	if (backw == true) dt = -1.0*dt;
	ddt = dt;
	CNST = 0.5*q/m*dt;
	behind();
}

void Mover::finalize(double& gamma) {
	ahead();
	gamma = this->gamma;
}

void Mover::ahead() {
	static double t_mag;
	static double v_dot[3];
	static double v_minus[3];
	static double t[3];
	static double s[3];
	static double v_plus[3];
	static const double C2 = SPEED_OF_LIGHT*SPEED_OF_LIGHT;
	v_minus[0] = v[0] + CNST*E[0];
	v_minus[1] = v[1] + CNST*E[1];
	v_minus[2] = v[2] + CNST*E[1];
	gamma = sqrt(1.0 + (v_minus[0]*v_minus[0]+v_minus[1]*v_minus[1]+v_minus[2]*v_minus[2]) / C2);
	t[0] = CNST*B[0]/gamma;
	t[1] = CNST*B[1]/gamma;
	t[2] = CNST*B[2]/gamma;
	t_mag = t[0]*t[0] + t[1]*t[1] + t[2]*t[2];
	s[0] = t[0] * 2 / (1 + t_mag);
	s[1] = t[1] * 2 / (1 + t_mag);
	s[2] = t[2] * 2 / (1 + t_mag);
	v_dot[0] = v_minus[0] + v_minus[1]*t[2] - v_minus[2]*t[1];
	v_dot[1] = v_minus[1] + v_minus[2]*t[0] - v_minus[0]*t[2];
	v_dot[2] = v_minus[2] + v_minus[0]*t[1] - v_minus[1]*t[0];
	v_plus[0] = v_minus[0] + v_dot[1]*s[2] - v_dot[2]*s[1];
	v_plus[1] = v_minus[1] + v_dot[2]*s[0] - v_dot[0]*s[2];
	v_plus[2] = v_minus[2] + v_dot[0]*s[1] - v_dot[1]*s[0];
	v[0] = v_plus[0];
	v[1] = v_plus[1];
	v[2] = v_plus[2];
	gamma = sqrt(1.0 + (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]) / C2);
}

void Mover::behind() {
	static double t_mag;
	static double v_dot[3];
	static double v_minus[3];
	static double t[3];
	static double s[3];
	static double v_plus[3];
	static const double C2 = SPEED_OF_LIGHT*SPEED_OF_LIGHT;
	v_minus[0] = v[0];
	v_minus[1] = v[1];
	v_minus[2] = v[2];
	gamma = sqrt(1.0 + (v_minus[0]*v_minus[0]+v_minus[1]*v_minus[1]+v_minus[2]*v_minus[2]) / C2);
	t[0] = -CNST*B[0]/gamma;
	t[1] = -CNST*B[1]/gamma;
	t[2] = -CNST*B[2]/gamma;
	t_mag = t[0]*t[0] + t[1]*t[1] + t[2]*t[2];
	s[0] = t[0] * 2 / (1 + t_mag);
	s[1] = t[1] * 2 / (1 + t_mag);
	s[2] = t[2] * 2 / (1 + t_mag);
	v_dot[0] = v_minus[0] + v_minus[1]*t[2] - v_minus[2]*t[1];
	v_dot[1] = v_minus[1] + v_minus[2]*t[0] - v_minus[0]*t[2];
	v_dot[2] = v_minus[2] + v_minus[0]*t[1] - v_minus[1]*t[0];
	v_plus[0] = v_minus[0] + v_dot[1]*s[2] - v_dot[2]*s[1];
	v_plus[1] = v_minus[1] + v_dot[2]*s[0] - v_dot[0]*s[2];
	v_plus[2] = v_minus[2] + v_dot[0]*s[1] - v_dot[1]*s[0];
	v[0] = v_plus[0] - CNST*E[0];
	v[1] = v_plus[1] - CNST*E[1];
	v[2] = v_plus[2] - CNST*E[2];
	gamma = sqrt(1.0 + (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]) / C2);
}

void Mover::push() {
	static double t_mag;
	static double v_dot[3];
	static double v_minus[3];
	static double t[3];
	static double s[3];
	static double v_plus[3];
	static const double C2 = SPEED_OF_LIGHT*SPEED_OF_LIGHT;
	
	v_minus[0] = v[0] + CNST*E[0];
	v_minus[1] = v[1] + CNST*E[1];
	v_minus[2] = v[2] + CNST*E[2];
	gamma = sqrt(1.0 + (v_minus[0]*v_minus[0]+v_minus[1]*v_minus[1]+v_minus[2]*v_minus[2]) / C2);
	t[0] = CNST*B[0]/gamma;
	t[1] = CNST*B[1]/gamma;
	t[2] = CNST*B[2]/gamma;
	t_mag = t[0]*t[0] + t[1]*t[1] + t[2]*t[2];
	s[0] = t[0] * 2 / (1 + t_mag);
	s[1] = t[1] * 2 / (1 + t_mag);
	s[2] = t[2] * 2 / (1 + t_mag);
	v_dot[0] = v_minus[0] + v_minus[1]*t[2] - v_minus[2]*t[1];
	v_dot[1] = v_minus[1] + v_minus[2]*t[0] - v_minus[0]*t[2];
	v_dot[2] = v_minus[2] + v_minus[0]*t[1] - v_minus[1]*t[0];
	v_plus[0] = v_minus[0] + v_dot[1]*s[2] - v_dot[2]*s[1];
	v_plus[1] = v_minus[1] + v_dot[2]*s[0] - v_dot[0]*s[2];
	v_plus[2] = v_minus[2] + v_dot[0]*s[1] - v_dot[1]*s[0];
	v[0] = v_plus[0] + CNST*E[0];
	v[1] = v_plus[1] + CNST*E[1];
	v[2] = v_plus[2] + CNST*E[2];
	gamma = sqrt(1.0 + (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]) / C2);
	r[0] += v[0]*ddt/gamma;
	r[1] += v[1]*ddt/gamma;
	r[2] += v[2]*ddt/gamma;
}












