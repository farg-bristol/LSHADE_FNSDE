
# FORTRAN code for neighbourhood-based speciation niching using L-SHADE

 + Author: Daniel Poole
 + Date: July 2018
 + Version: 1.0

## Summary and reference

FORTRAN code for neighbourhood-based speciation niching using L-SHADE, with feasible selection for constraints. Test suite is newly developed constrained multimodal functions. The reference for this code is:

[D. J. Poole and C. B. Allen. "Constrained Niching Using Differential Evolution". Swarm and Evolutionary Computation, 2018.](https://www.sciencedirect.com/science/article/pii/S2210650217307344?via%3Dihub)

In is kindly requested that any use of the code is referenced to the above paper.

NOTE: the author takes no responsibility for the results provided by the code under the stewardship of another.

Please email bugs to d.j.poole@bristol.ac.uk

## Usage

Code is executed from the RUNDATA directory:

	cd <code directory>/LSHADE_FNSDE
	make
	cd RUNDATA
	../de.exe

Edit DE.conf file to change parameters.  Dat file contains information read by code (DO NOT EDIT UNLESS CHANGING SOURCE CODE).


### Files: 

Language of the code is FORTRAN.

 + de: main optimizer
 + de_config: read configuration file
 + MAIN: main function
 + niching_func_cons: constrained anaytical functions
 + perfanal: post processing subroutine
  
Code has been tested in following environment:

### Linux
 + OS: Ubuntu 14.04.5 LTS 64-bit
 + CPU: core i7 (2.20GHz)
 + RAM: 8GB
 + Compiler: gfortran 4.8.4 with -O2 optimization flag

### Windows
 + OS: Windows 7 64-bit
 + CPU: core i7 (3.80GHz)
 + RAM: 16GB
 + Language: fortran
 + Compiler: GNU Fortran MinGW-w64 -O2 optimization flag


# Copyright 2018 Daniel Poole

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
