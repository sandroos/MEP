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

#ifndef INJECTOR_H
#define INJECTOR_H
#include <vector>

struct Instrument {
	double e[3];             /**< A unit vector showing the instrument orientation. Points 
							  * to the inside of the instrument. */
	double coneWidth;        /**< Width of the observation cone.*/
	double cos_cone;         /**< Cosine(observation cone width). */
};

class Injector {
public:
	Injector();
	~Injector();
	
	bool getNewParticle(double& t,double& q,double& m,double& g,double r[3],double v[3]);
private:
	double q_inj;              /**< Charge of the injected particles. */
	double m_inj;              /**< Mass of the injected particles. */
	double W_inj_min;          /**< Minimum energy of the injected particles. */
	double W_inj_max;          /**< Maximum energy of the injected particles. */
	double pwrlawindex;
	unsigned int N_injected;
	unsigned int N_particles;
	double r_inj[3];           /**< The coordinates of the injection (spacecraft) position. */
	//double ori[3];             /**< A unit vector to the direction of the exit hole of the spacecraft. */

	//double spacecraft[3];      /**< Position of the spacecraft. */
	double unitnadir[3];       /**< A unit vector to the direction of the nadir. */
	double unitvel[3];         /**< A unit vector to the direction of the spacecraft velocity. */
	double unity[3];
	std::vector<Instrument> instruments; 
	int instr;
	
	//double e1[3];
	//double e2[3];
	//double e3[3];
	double y_min;
	double y_max;
	double norm_energy;
	
	void writeGnuplot();
};

#endif
