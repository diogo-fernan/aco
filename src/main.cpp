/*
@author:	Diogo A. B. Fernandes
@contact:	diogoabfernandes@gmail.com
@license:	see LICENSE
*/

#include <iostream>
#include <cstdlib>

#include <cmath>
#include <limits>
#include <climits>

#include "ACO.h"

#define ITERATIONS		(int) 5

#define NUMBEROFANTS	(int) 4
#define NUMBEROFCITIES	(int) 8

// if (ALPHA == 0) { stochastic search & sub-optimal route }
#define ALPHA			(double) 0.5
// if (BETA  == 0) { sub-optimal route }
#define BETA			(double) 0.8
// Estimation of the suspected best route.
#define Q				(double) 80
// Pheromones evaporation. 
#define RO				(double) 0.2
// Maximum pheromone random number.
#define TAUMAX			(int) 2

#define INITIALCITY		(int) 0

int main() {

	ACO *ANTS = new ACO (NUMBEROFANTS, NUMBEROFCITIES, 
			 			ALPHA, BETA, Q, RO, TAUMAX,
			 			INITIALCITY);

	ANTS -> init();

	ANTS -> connectCITIES (0, 1);
	ANTS -> connectCITIES (0, 2);
	ANTS -> connectCITIES (0, 3);
	ANTS -> connectCITIES (0, 7);
	ANTS -> connectCITIES (1, 3);
	ANTS -> connectCITIES (1, 5);
	ANTS -> connectCITIES (1, 7);
	ANTS -> connectCITIES (2, 4);
	ANTS -> connectCITIES (2, 5);
	ANTS -> connectCITIES (2, 6);
	ANTS -> connectCITIES (4, 3);
	ANTS -> connectCITIES (4, 5);
	ANTS -> connectCITIES (4, 7);
	ANTS -> connectCITIES (6, 7);
	/* ANTS -> connectCITIES(8, 2);
	ANTS -> connectCITIES(8, 6);
	ANTS -> connectCITIES(8, 7); */

	ANTS -> setCITYPOSITION (0,  1,  1);
	ANTS -> setCITYPOSITION (1, 10, 10);
	ANTS -> setCITYPOSITION (2, 20, 10);
	ANTS -> setCITYPOSITION (3, 10, 30);
	ANTS -> setCITYPOSITION (4, 15,  5);
	ANTS -> setCITYPOSITION (5, 10,  1);
	ANTS -> setCITYPOSITION (6, 20, 20);
	ANTS -> setCITYPOSITION (7, 20, 30);
	// ANTS -> setCITYPOSITION(8, 26, 20);

	ANTS -> printGRAPH ();

	ANTS -> printPHEROMONES ();

	ANTS -> optimize (ITERATIONS);

	ANTS -> printRESULTS ();

	return 0;
}
