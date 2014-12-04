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

#ifndef SCALING_H
#define SCALING_H

class Scaling {
public:
	Scaling();
	~Scaling();
	
	static real J_to_dim(const real& J) {return J * DIST*DIST * TIME / CHARGE;}
	static real E_to_dim(const real& E) {return E * TIME*TIME * CHARGE / MASS / DIST;}
	static real B_to_dim(const real& B) {return B * TIME * CHARGE / MASS;}
	static real V_to_dim(const real& V) {return V * TIME / DIST;}
	static real L_to_dim(const real& L) {return L / DIST;}
	static real T_to_dim(const real& T) {return T / TIME;}
	static real Q_to_dim(const real& Q) {return Q / CHARGE;}
	static real M_to_dim(const real& M) {return M / MASS;}
	static real RHOM_to_dim(const real& RHOM) {return RHOM * DIST*DIST*DIST / MASS;}
	static real RHOQ_to_dim(const real& RHOQ) {return RHOQ * DIST*DIST*DIST / CHARGE;}
	static real W_to_dim(const real& W) {return W * TIME*TIME / MASS / DIST/DIST;}
	static real N_to_dim(const real& n) {return n * DIST*DIST*DIST;}
	static real P_to_dim(const real& P) {return P * DIST * TIME*TIME / MASS;}
	static real I_to_dim(const real& I) {return I * TIME / CHARGE;}
	
	static real J_to_phy(const real& J) {return J / DIST/DIST / TIME * CHARGE;}
	static real E_to_phy(const real& E) {return E / TIME/TIME / CHARGE * MASS * DIST;}
	static real B_to_phy(const real& B) {return B / TIME / CHARGE * MASS;}
	static real V_to_phy(const real& V) {return V / TIME * DIST;}
	static real L_to_phy(const real& L) {return L * DIST;}
	static real T_to_phy(const real& T) {return T * TIME;}
	static real Q_to_phy(const real& Q) {return Q * CHARGE;}
	static real M_to_phy(const real& M) {return M * MASS;}
	static real RHOM_to_phy(const real& RHOM) {return RHOM * MASS / DIST/DIST/DIST;}
	static real RHOQ_to_phy(const real& RHOQ) {return RHOQ * CHARGE / DIST/DIST/DIST;}
	static real W_to_phy(const real& W) {return W * MASS * DIST*DIST / TIME/TIME;}
	static real N_to_phy(const real& n) {return n / DIST/DIST/DIST;}
	static real P_to_phy(const real& P) {return P * MASS / DIST / TIME/TIME;}
	static real I_to_phy(const real& I) {return I * CHARGE / TIME;}
	
private:
	double DIST;
	double MASS;
	double TIME;
	double CHARGE;
};

#endif
