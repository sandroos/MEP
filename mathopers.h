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

#ifndef MATHOPERS_H
#define MATHOPERS_H

#include <cmath>

template<typename T> void nullVector(T v[3]) {
	v[0] = 0.0;
	v[1] = 0.0;
	v[2] = 0.0;
}

template<typename T> T dotProduct(const T v1[3],const T v2[3]) {
	    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

template<typename T> T vectorMagnitude(const T v[3]) {
	return sqrt(dotProduct(v,v));
}

template<typename T> T vectorMagnitude2(const T v[3]) {
	return dotProduct(v,v);
}

template<typename T> void unitVector(T v[3]) {
	T mag = vectorMagnitude(v);
	v[0] /= mag;
	v[1] /= mag;
	v[2] /= mag;
}

template<typename T> void unitVector(const T v[3],T result[3]) {
	T mag = vectorMagnitude(v);
	result[0] = v[0] / mag;
	result[1] = v[1] / mag;
	result[2] = v[2] / mag;
}

template<typename T> void crossProduct(const T vec1[3],const T vec2[3],T result[3]) {
	result[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
	result[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
	result[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
}

#endif
