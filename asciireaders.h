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

#ifndef ASCIIREADERS_H
#define ASCIIREADERS_H

#include <cstdlib>
#include <fstream>
#include <vector>

//void parseValues(const std::string& line,std::vector<float>& values);
//bool readline(std::fstream& in,std::vector<float>& values);
//bool readline(std::fstream& in,std::vector<double>& values);
//bool readline(std::fstream& in,std::vector<std::string>& values);

/** DECLARATIONS */

template<typename T> void parseValues(const std::string& line,std::vector<T>& values);
template<typename T> bool readline(std::iostream& in,std::vector<T>& values);
bool readline(std::iostream& in,std::string& line);

/** DEFINITIONS */

template<typename T> void parseValues(const std::string& line,std::vector<T>& values) {
	std::string number;
	std::string::size_type i=0, j=0;
	values.clear();
	
	while (i < line.size()) {
		while ((line[i] == ' ' || line[i] == '\t') && i < line.size()) ++i;
		j = i;
		while ((line[i] != ' ' && line[i] != '\t') && i < line.size()) ++i;
		number = line.substr(j,i-j);
		values.push_back(atof(number.c_str()));
	}
}

// Specialized template version of parseValues for strings.
template<> void parseValues<std::string>(const std::string& line,std::vector<std::string>& values) {
	std::string::size_type i=0, j=0;
	values.clear();
	
	while (i < line.size()) {
		while ((line[i] == ' ' || line[i] == '\t') && i < line.size()) ++i;
		j = i;
		while ((line[i] != ' ' && line[i] != '\t') && i < line.size()) ++i;
		values.push_back(line.substr(j,i-j));
	}
}

template<typename T> bool readline(std::iostream& in,std::vector<T>& values) {
	const int n = 4096;
	char line[n];
	char number[n];
	int i = 0,j,chars;
	
	/** Read a new line until end of file is reached, the line is not a comment or empty line. */
	in.getline(line,n);
	chars = in.gcount();
	
	do {
		if (in.good() == false) return false; // Read failed.
		
		while (line[i] == ' ' || line[i] == '\t') {
			++i; // Skip empty chars.
		}
		
		// End of line reached, read a new line.
		if (line[i] == '\n' || line[i] == '\0') {
			in.getline(line,n);
			chars = in.gcount();
			continue;
		}
		
		if (line[i] == '#' || line[i] == '%') { // Check for comment line.
			in.getline(line,n);
			chars = in.gcount();
			continue;
		} else break;
		
	} while (true);
	
	/** At this point line[i] points to the first interesting char. */
	while (i < chars) {
		while (line[i] == ' ' || line[i] == '\t') ++i;
		if (line[i] == '\n' || line[i] == '\0') return true;
		
		j = 0;
		while (line[i] != ' ' && line[i] != '\t') {
			if (line[i] == '\n' || line[i] == '\0') break;
			number[j] = line[i];
			++i;
			++j;
		}
		number[j] = '\0';
		values.push_back(atof(number));
	}	
	return true;
}

bool readline(std::iostream& in,std::string& line) {
	const int n = 4096;
	char cline[n];

	if (in.good() == false) {
		line = "";
		return false;
	} else {
		in.getline(cline,n);
		line = cline;
		return true;
	}
}

// Specialized template version of readline for strings.
template<> bool readline<std::string>(std::iostream& in,std::vector<std::string>& values) {
	const int n = 4096;
	char line[n];
	char number[n];
	int i = 0,j,chars;
	
	/** Read a new line until end of file is reached, the line is not a comment or empty line. */
	in.getline(line,n);
	chars = in.gcount();
	
	do {
		if (in.good() == false) return false; // Read failed.
		
		while (line[i] == ' ' || line[i] == '\t') ++i; // Skip empty chars.
		
		// End of line reached, read a new line.
		if (line[i] == '\n' || line[i] == '\0') {
			in.getline(line,n);
			chars = in.gcount();
			continue;
		}
		
		if (line[i] == '#' || line[i] == '%') { // Check for comment line.
			in.getline(line,n);
			chars = in.gcount();
			continue;
		} else break;
				
	} while (true);
	
	/** At this point line[i] points to the first interesting char. */
	while (i < chars) {
		while (line[i] == ' ' || line[i] == '\t') ++i;
		if (line[i] == '\n' || line[i] == '\0') return true;
		
		j = 0;
		while (line[i] != ' ' && line[i] != '\t') {
			if (line[i] == '\n' || line[i] == '\0') break;
			number[j] = line[i];
			++i;
			++j;
		}
		number[j] = '\0';
		values.push_back(std::string(number));
	}	
	return true;
}

#endif
