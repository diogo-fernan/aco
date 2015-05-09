#include "Randoms.cpp"

class ACO {
public:
	ACO (int nAnts, int nCities, 
		double alpha, double beta, double q, double ro, double taumax,
		int initCity);
	virtual ~ACO ();
	
	void init ();
	
	void connectCITIES (int cityi, int cityj);
	void setCITYPOSITION (int city, double x, double y);
	
	void printPHEROMONES ();
	void printGRAPH ();
	void printRESULTS ();
	
	void optimize (int ITERATIONS);

private:
	double distance (int cityi, int cityj);
	bool exists (int cityi, int cityc);
	bool vizited (int antk, int c);
	double PHI (int cityi, int cityj, int antk);
	
	double length (int antk);
	
	int city ();
	void route (int antk);
	int valid (int antk, int iteration);
	
	void updatePHEROMONES ();

	
	int NUMBEROFANTS, NUMBEROFCITIES, INITIALCITY;
	double ALPHA, BETA, Q, RO, TAUMAX;
	
	double BESTLENGTH;
	int *BESTROUTE;

	int **GRAPH, **ROUTES;
	double **CITIES, **PHEROMONES, **DELTAPHEROMONES, **PROBS;

	Randoms *randoms;
};

