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

#ifndef PARAMREADER_H
#define PARAMREADER_H

#include <cstdlib>
#include <iostream>
#include <map>

class ParamReader {
public:
	ParamReader(const std::string& paramfile);
	~ParamReader();
	
	static std::map<std::string,std::string> getParameters(const std::string& id);
	
private:
	ParamReader();
	
	static std::map<std::string,std::map<std::string,std::string> > params;
};

#endif
