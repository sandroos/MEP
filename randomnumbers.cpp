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
#include <ctime>

#include "randomnumbers.h"
#include "paramreader.h"

using namespace std;

// Initialize static member variables.
#ifdef DEBUG
  bool RandomNumbers::initialized = false;
#endif
int RandomNumbers::seed = -87623;

double RandomNumbers::ran2(int &idum) {
   return (1.0*rand())/RAND_MAX;
}

RandomNumbers::RandomNumbers() {
	// Get the seed value from the parameter file.
	map<string,string> params = ParamReader::getParameters("RandomNumbers");
	if (params.find("seed") == params.end()) {
		cerr << "WARNING(RandomNumbers): Could not find parameter 'seed' in parameter file!" << endl;
      seed = time(NULL);
		cerr << "\t Using value " << seed << " from time(NULL)" << endl;
	} else {
		seed = atoi(params["seed"].c_str());
      if (seed < 0) seed = -1*seed;
	}
   srand(seed);
   #ifdef DEBUG
		  initialized = true;
	#endif
}

RandomNumbers::RandomNumbers(const int& seed) {
	#ifdef DEBUG
	  if (initialized == true) return;
	  initialized = true;
	#endif
	this->seed = seed;
   if (this->seed < 0) this->seed = -1*this->seed;
   srand(seed);
}

RandomNumbers::~RandomNumbers() { }

double RandomNumbers::uniform() {
	#ifdef DEBUG
	  if (initialized == false) {
		  cerr << "ERROR(RandomNumbers): Class not initialized, uniform() called." << endl;
		  exit(1);
	  }
	#endif
	return ran2(seed);
}
