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

#ifndef MOVER_H
#define MOVER_H

class Mover {
public:
	Mover();
	~Mover();
	
	void push();
	void initialize(const bool& backw,double& dt,const double& q,const double& m,double* const r,double* const v,double* const E,double* const B);
	void finalize(double& gamma);
private:
	double* r;
	double* v;
	double* E;
	double* B;
	double gamma;
	double CNST;
	
	void ahead();
	void behind();
};

#endif
