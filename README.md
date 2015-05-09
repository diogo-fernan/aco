# aco

*aco* is an ISO C++ Ant Colony Optimization (ACO) algorithm (a metaheuristic optimization technique inspired on ant behavior) for the traveling salesman problem. It releases a number of ants incrementally whilst updating pheromone concentration and calculating the best graph route. In the end, the best route is printed to the command line.

*aco* was developed a few years back for academic and research purposes. It is available for anyone interested in ACO methods. Feel free to explore the algorithm and tweak its parameteres.

# Usage

*aco* was developed with standard libraries and an additional `Randoms.cpp` file that implements generation methods of pseudo-random numbers. Apart from that, the program can be easily compiled and run with the following command at a command line located in `src`:

* `g++ -Wall *.cpp -o aco; ./aco`.

The `main` method is included in the `main.cpp` file which contains a few parameters of the ACO algorithm that can be changed as one sees fit. That file also defines programmatically the connections of the cities that ultimately make up the city graph. The ACO algorithm methods are implemented in the `ACO` class located in the `ACO.h` and `ACO.cpp` files. 