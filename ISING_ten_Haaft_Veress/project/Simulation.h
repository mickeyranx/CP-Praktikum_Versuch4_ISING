#pragma once
#include <vector>
#include <random>
#include <iostream>
using namespace std;
class Simulation {
public:
    
    //tools for randomisation
    static mt19937 gen;
    static uniform_real_distribution<double> distr1;
    static uniform_int_distribution<int> distr2;
    static random_device rand_device;

    //parameters for testing multi-hit
    static unsigned int try_count;
    static unsigned int acceptance_count;
    
    //coldstart
    static vector<int> initializeLatticeCold(int L);

    //hotstart
    static vector<int> initializeLatticeHot(int L);

    //return a small list of 4 indices of Positions of the neighbors (top, right, bottom, left)
    static vector<int> getNeighborPos(int i, int L);

    //change in energy in case of a spin flip (otherwise the energy does not change)
    static double changeInEnergy(vector<int>& config_1, int L, int pos, int s, double h);

    
    //Metropolis sweep with multi-hit
    static void sweepMetropolisMultihit(vector<int>& config_1, int L, double beta, double h, int tries);

    //Heathbath sweep 
    static void sweepHeatbath(vector<int>& config_1, int L, int M, double beta, double h);

    //this is method is for Metropolis multihit, it performs <draw_interval> sweeps
    static void draw(vector<int>& config, int L, double beta, double h, int draw_interval, int multi_hit);

    //this method is for Heathbath, it performs <draw_interval> sweeps
    static void draw(vector<int>& config, int L, double beta, double h, int draw_interval);

    //calcualtes average inner energy per lattice point
    static double averageEnergy(vector<int>& config, int K, int L, double h);

    //calcuates average Magnetisation (!!!not per point!!!)
    static double averageMagnetisation(int M, vector<int>& config);

    
    //prints current state of the lattice the console
    static void printConfig(vector<int>& config);

};

