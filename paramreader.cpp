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

#include <fstream>
#include <vector>
#include "paramreader.h"
#include "asciireaders.h"

using namespace std;

map<string,map<string,string> > ParamReader::params;

ParamReader::ParamReader(const std::string& paramfile) {
	fstream in(paramfile.c_str(),fstream::in);
	if (in.good() == false) {
		cerr << "ERROR(ParamReader): Could not open file '" << paramfile << "' for reading!" << endl;
		exit(1);
	}
	
	// Go through all lines in the parameter file. Insert the option-value pairs for 
	// each identifier into map params. 
	string current_id;
	map<string,string> tmp_map;
	vector<string> line;
	int depth = 0;
	while (readline(in,line) == true) {
      if (line[0] == "}" && depth == 0) {
         cerr << "ERROR: Parameter file contained '}' with no opening braces." << endl;
         exit(1);
      } else if (line[0] == "//") {
         // Skip comment line
         line.clear();
         continue;
      } else if (line.size() == 1 && line[0] != "}") {
         cerr << "ERROR: Each line in parameter file must contain 2 strings, except for the ending brackets." << endl;
         exit(1);
      }
		
		if (line[0] == "}") {
			depth = 0;
			params[current_id] = tmp_map;
			tmp_map.clear();
		} else if (line[1] == "{") {
			current_id = line[0];
			++depth;
		} else {
			tmp_map[line[0]] = line[1];
		}
		line.clear();
	}
	
	in.close();
}

ParamReader::~ParamReader() { }

std::map<std::string,std::string> ParamReader::getParameters(const std::string& id) {
	if (params.find(id) == params.end()) return map<string,string>(); // ID not found, return empty map.
	else return params[id];
}






