#This file is part of MEP.
#
#    Copyright 2014 Arto Sandroos.
#
#    MEP is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    MEP is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with MEP.  If not, see <http://www.gnu.org/licenses/>.

# Name of the C++ compiler
CMP=g++

# Name of the linker
LNK=g++

# Default compiler options, do not touch
CXXFLAGS=-O3 -Wall

# User-defined compiler options, set these on the command line
FLAGS=

# Header file include(s) for the compiler
INC=

# Library include(s) for the linker
LIB=

# Headers needed for compilation
HDRS = paramreader.h mover.h randomnumbers.h injector.h mathopers.h\
	emfield.h

# Object files needed for compilation
OBJS = paramreader.o mover.o randomnumbers.o main.o injector.o

# Magnetospheric model used
FIELD = mirdipsheet.o

# Make targets

default:
	$(MAKE) all

clean:
	rm -f mercury flines writesphere *.o *~

all:
	$(MAKE) mercury
	$(MAKE) flines
	$(MAKE) writesphere

# Make rules

.cpp.o: $(HDRS)
	$(CMP) $(CXXFLAGS) $(FLAGS) $(INC) -c $<

mercury: $(HDRS) $(OBJS) $(FIELD)
	$(LNK) -o mercury $(OBJS) $(FIELD) $(INC) $(LIB)

flines: flines.cpp paramreader.o $(FIELD)
	$(CMP) $(CXXFLAGS) $(FLAGS) -c flines.cpp $(INC)
	$(LNK) -o flines flines.o paramreader.o $(FIELD) $(LIB)

writesphere: writesphere.cpp
	$(CMP) $(CXXFLAGS) $(FLAGS) -c writesphere.cpp
	$(LNK) -o writesphere writesphere.o $(LIB)

