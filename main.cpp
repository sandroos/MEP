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
#include <fstream>
#include <map>

#include "constants.h"
#include "mathopers.h"
#include "paramreader.h"
#include "randomnumbers.h"
#include "mover.h"
#include "injector.h"
#include "emfield.h"

using namespace std;

#ifdef TRAJ
  const double RM = MERCURY_EQUAT_RADIUS;
#endif

int main(int argn,char* args[]) {
	if (argn != 3) {
		cerr << "USAGE: ./mercury <parameter file> <output file>" << endl;
		cerr << endl;
		return 1;
	}
	ParamReader paramReader(args[1]);
	RandomNumbers randomNumbers;
	EMfield field;
	Mover mover;
	Injector injector;
	fstream out(args[2],fstream::out);
	if (out.good() == false) {
		cerr << "ERROR(main): Could not open '" << args[2] << "' for output!" << endl;
		return(1);
	}
	
	/** The idea here is that we get a particle from an injector (r_ini,v_ini,q,m,gamma). The 
	 * particle trajectory is then calculated forward and backward in time using temporary 
	 * variables (r,v,t). When the forward iteration ends (as indicated by iterate == false), 
	 * the forward exit state is stored in (t_fw,r_fw,v_fw), and the iteration is restarted 
	 * using the same injection state (r_ini,v_ini,...), but this time we're iterating backwards in time. 
	 * The "past" exit state is stored in (t_bw,r_bw,v_bw). 
	 * 
	 * However, some of the injected particles may be trapped in the magnetic field, and we also need 
	 * an exit condition with respect to simulation time -- otherwise the simulation will just 
	 * continue indefinitely. The maximum simulation time is stored in variable maxtime_per_dt, 
	 * which basically is just the maximum number of iterations allowed. If the simulation time exceeds 
	 * the maximum allowed time, the particle will be rejected.
	 */
	
	double t;         /**< Current simulation time for each particle. t=0 is the injection time. */
	double dt;        /**< Timestep used (and calculated by) the particle mover, positive or negative. */
	double q;         /**< Charge (C). */
	double m;         /**< Mass (kg). */
	double gamma;     /**< Gamma-factor. */
	double r[3];      /**< Position of the particle in cartesian coordinates (m), Mercury is at the origin. */
	double v[3];      /**< Velocity of the particle (m/s). */
	double E[3];      /**< Electric field at position r (N/C). */
	double B[3];      /**< Magnetic field at position r (T). */
	bool iterate;     /**< When set to false, the particle iteration should be stopped. */
	
	double r_ini[3]; /**< Injection position of the particle. */
	double v_ini[3]; /**< Injection velocity of the particle. */
	double t_fw=0.0; /**< Exit time when iterating forward in time. */
	double r_fw[3];  /**< Exit position -- "" -- .*/
	double v_fw[3];  /**< Exit velocity -- "" -- .*/
	double gamma_fw; /**< Exit gamma factor      .*/
	double t_bw=0.0; /**< Exit time when iterating backward in time. */
	double r_bw[3];  /**< Exit position -- "" -- .*/
	double v_bw[3];  /**< Exit velocity -- "" -- .*/
	double gamma_bw; /**< Exit gamma factor      .*/
	
	bool acceptParticle;              /**< If true, the particle has exited the simulation box both when calculating 
									   * forward and backward in time, and thus should be written into output file. */

	unsigned int N_discarded = 0;     /**< Count the number of particle that we're not accepted, in practice these 
									   * are trapped particles. */

	double maxtime_per_dt;            /**< Maximum number of iterations allowed, i.e., an exit condition with respect to time. */
	double MAXTIME;                   /**< abs(maximum simulation time), in seconds. */
	
	#ifdef TRAJ
	  unsigned int saveinterval;      /**< How many iterations between consecutive writes to the output file? */
	  unsigned int t_counter;         /**< This is used to calculate the number of iterations taken. */
	  unsigned int fw_saves = 0;      /**< The number of "future" states written into the output file. */
	  unsigned int bw_saves = 0;      /**< The number of "past" states written into the output file. */
	#endif
	{
		/** Read the parameter values from the "Simulation" brackets in the parameter file. */
		map<string,string> params = ParamReader::getParameters("Simulation");
		if (params.find("maxtime_per_dt") == params.end()) {
			cerr << "ERROR(main): Parameter 'maxtime_per_dt' was not found!" << endl;
			return 1;
		}
		maxtime_per_dt = atof(params["maxtime_per_dt"].c_str());
		#ifdef TRAJ
		  if (params.find("saveinterval_dt") == params.end()) {
			  cerr << "ERROR(main): Parameter 'saveinterval_dt' was not found." << endl;
			  return 1;
		  }
		  saveinterval = atoi(params["saveinterval_dt"].c_str());
		  out.precision(3);
		  out << scientific;
		#endif
	}
	
	while (injector.getNewParticle(t,q,m,gamma,r_ini,v_ini) == true) {
		// We have now obtained a new particle from the injector. The trajectory should now 
		// be calculated forward and backward in time.
		acceptParticle = true;
		for (int dir=0; dir<2; ++dir) {
			t = 0.0;
			for (int i=0; i<3; ++i) {
				r[i] = r_ini[i];
				v[i] = gamma*v_ini[i]; // Relativistic mover uses gamma*velocity. 
			}
			nullVector(E);
			nullVector(B);
			field.values(t,r,E,B);
			iterate = true;
			switch (dir) {
			case 0: // Forward in time.
				#ifdef TRAJ
				  t_counter = 0;
				#endif
				mover.initialize(false,dt,q,m,r,v,E,B);
				MAXTIME = fabs(maxtime_per_dt*dt);
				while (iterate == true) {
					#ifdef TRAJ
					  if (t_counter % saveinterval == 0) {
						  out << t << ' ' << r[0]/RM << ' ' << r[1]/RM << ' ' << r[2]/RM << ' ';
						  out << v[0] << ' ' << v[1] << ' ' << v[2] << ' ';
						  out << B[0] << ' ' << B[1] << ' ' << B[2] << endl;
						  ++fw_saves;
					  }
					  ++t_counter;
					#endif
					mover.push();
					t += dt;
					iterate = field.values(t,r,E,B);
					if (t > MAXTIME) {
						iterate = false;
						acceptParticle = false;
					}
				}
				mover.finalize(gamma_fw);
				t_fw = t;
				for (int i=0; i<3; ++i) {
					r_fw[i] = r[i];
					v_fw[i] = v[i];
				}
				#ifdef TRAJ
				  out << endl;
				#endif
				break;
			case 1: // Backward in time.
				#ifdef TRAJ
				  t_counter = 0;
				#endif
				mover.initialize(true,dt,q,m,r,v,E,B);
				MAXTIME = fabs(maxtime_per_dt*dt);
				while (iterate == true) {
					#ifdef TRAJ
					  if (t_counter % saveinterval == 0) {
						  out << t << ' ' << r[0]/RM << ' ' << r[1]/RM << ' ' << r[2]/RM << ' ';
						  out << v[0] << ' ' << v[1] << ' ' << v[2] << ' ';
						  out << B[0] << ' ' << B[1] << ' ' << B[2] << endl;
						  ++bw_saves;
					  }
					  ++t_counter;
					#endif
					mover.push();
					t += dt;
					iterate = field.values(t,r,E,B);
					if (t < -MAXTIME) {
						iterate = false;
						acceptParticle = false;
					}
				}
				mover.finalize(gamma_bw);
				t_bw = t;
				for (int i=0; i<3; ++i) {
					r_bw[i] = r[i];
					v_bw[i] = v[i];
				}
				break;
			}
		}
		if (acceptParticle == false) {
			++N_discarded;
			#ifdef TRAJ
			  break;
			#endif
			continue;
		}
		#ifdef TRAJ
		  // This one is here to prevent gnuplot from connecting all points with lines 
		  // in the case that fw_saves == bw_saves.
		  if (fw_saves % bw_saves == 0) {
			  out << t << ' ' << r[0]/RM << ' ' << r[1]/RM << ' ' << r[2]/RM << ' ';
			  out << v[0] << ' ' << v[1] << ' ' << v[2] << ' ';
			  out << B[0] << ' ' << B[1] << ' ' << B[2] << endl;
		  }
		  break;
		#endif
		// Forward and backward t,r,v have now been calculated. The states are now written into 
		// the output files.
		out << gamma << ' ';
		for (int i=0; i<3; ++i) out << v_ini[i] << ' ';
		out << t_fw << ' ';
		for (int i=0; i<3; ++i) out << r_fw[i]/MERCURY_EQUAT_RADIUS << ' ';
		out << gamma_fw << ' ';
		for (int i=0; i<3; ++i) out << v_fw[i] << ' ';
		out << t_bw << ' ';
		for (int i=0; i<3; ++i) out << r_bw[i]/MERCURY_EQUAT_RADIUS << ' ';
		out << gamma_bw << ' ';
		for (int i=0; i<3; ++i) out << v_bw[i] << ' ';
		out << endl;
	}	
	out.close();
	
	cout << "Simulation complete." << endl;
	cout << "\t Discarded " << N_discarded << " particles." << endl;
	cout << endl;
}





