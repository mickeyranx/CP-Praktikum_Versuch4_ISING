#define _USE_MATH_DEFINES
#include <iostream>
#include <string>
#include <cmath>
#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <set>
#include <tuple>
#include <iomanip>
#include "Simulation.h"
#include <random>
#include <string>
#include <fstream>
#include <numeric>
#include <format> //works only with C++20


using namespace std;


//approximation of PI with MC (A1 a)
static double approximationOfPi(int N, int seed) {
    //initialise random tools
    mt19937 gen = mt19937(seed);
    uniform_real_distribution<double> unif_distr(-1.0, 1.0); //uniform double distribution between -1 and 1

    int square = 0; //counter for point generated insde the unit square
    int circle = 0; //counter for point gerneated inside the uni circle
    for (int i = 0; i < N; i++) //generate samples 
    {
        //roll for random x and y positions
        double x = unif_distr(gen);
        double y = unif_distr(gen);
        if (pow(x, 2) + pow(y, 2) <= 1) {
            circle++;
        }
        square++; //every point will be inside the unit square

    }

    return 4 * (double)circle / ((double)square); //the value of the ratio times 4 gives an approximation for pi depending on the sample size  

}

// RNG with the box müller method, samples a values from the normaldistribution with sigma and mu
static double boxMueller(double mu, double sigma, uniform_real_distribution<double> &unif_distr, mt19937 &gen) {
    //roll 2 numbers from a uniform distribution between 0 and 1 
    double m = unif_distr(gen);
    double n = unif_distr(gen);

    //calcualte new random variable
    double r = sigma * sqrt(-2.0 * log(m));
    return r * cos(2.0 * M_PI * n) + mu;

}

// approximation of a certain gaussian Integral with MC Integration (A1 b)
static double approxGaussianIntegral(int N) {
    //initialise random tools
    mt19937 gen = mt19937(38171);
    uniform_real_distribution<double> unif_distr(0.0001, 1.0); //uniform double distribution between -1 and 1
    double sigma = 1; //stdev
    double sum = 0;

    
    vector<double> x_values;
    //calculate with box-müller
    for (int i = 0; i < N; i++)
    {
        double x = boxMueller(0, sigma, unif_distr, gen); //mu is zero for the gaussian function we want to integrate 
        x_values.push_back(x);
    }
   
    for (double x : x_values)
    {
        double n = cos(x);
       
        sum += n;
    }

    //approximate with integration and summs
    return  sum / (N * std::pow(std::sqrt(2 * M_PI), 3));

}

//backtrack algorithm to invoke all possible configurations of a LxL lattice 
static void backtrack(int L, int counter, vector<int> config, vector<vector<int>>& list_of_configs) {
    if (counter >= pow(L, 2)) { //if the counter reaches the maximum amount of spins possible
        list_of_configs.push_back(config);
        return;
    }
    config.push_back(1);
    counter++;
    backtrack(L, counter, config, list_of_configs); //explore possibilites for 1 at this position until limit is reached
    
    
    counter--;
    config.pop_back(); //remove the previouls inserted 1 at this position
    config.push_back(-1);
    counter++;
    backtrack(L, counter, config, list_of_configs);//explore possibilites for -1 at this position until limit is reached
    return;

}

//calculates the mean energy pp, mean absolute magnetism pp 
//and mean magnetism pp explicitly with the partition function over all possible configurations
static void explicitIsing(int L) {
    //-----------------------------------
    //               setup
    //-----------------------------------
    //list of all possible configurations
    vector<vector<int>> list_of_configs = {};
    //create all configs with backtrack-algorithm
    backtrack(L, 0, {}, list_of_configs);
    cout << list_of_configs.size() << endl;
    //initialize list of betas
    vector<double> betas = {0.1, 0.2, 0.3, 0.35, 0.4, 0.415, 0.44, 0.4406868, 0.445, 0.45, 0.465, 0.48, 0.5, 0.55, 0.6, 0.7, 0.9};
    double K = (double) L * L;
    //-----------------------------------
    //      calculate observables
    //-----------------------------------
    ofstream File("explicit_vals_L=" + to_string(L) + ".txt");
    File << "beta" << "\t" << "<e>" << "\t" << "<m>" << "\t" << "<|m|>" << "\n";
    File << fixed << setprecision(7);
    for (double beta : betas) {
        //1.partition-function
        double Z = 0;
        vector<double> energies = {};
        for (vector<int> config : list_of_configs) {
            double H = Simulation::averageEnergy(config, L*L,L , 0 )*L*L;
            energies.push_back(H);
            Z += exp(-beta * H);
        }

        //2. mean energy pp, mean magnetism pp, mean absolute magnetism pp
        double mean_energy = 0;
        double mean_mag = 0;
        double mean_abs_mag = 0;
        int i = 0;
        for (vector<int> config : list_of_configs) {
            double H_i = energies[i];
            mean_energy += 1.0 / Z * exp(-beta * H_i) * H_i / ((double)L * L);
            double M_i = Simulation::averageMagnetisation(K, config);
            mean_mag += 1.0 / (K*Z) *  M_i * exp(-beta * H_i);
            mean_abs_mag += 1.0 / (K*Z) * abs(M_i) * exp(-beta * H_i);
            i++;
        }

        File << beta << "\t" << mean_energy << "\t" << mean_mag << "\t" << mean_abs_mag << "\n";


    }
    
    File.close();





}

//starts a single measurement (h, beta) with metropolis
static vector<double> startSimulationMetropolis(int L, double beta, double h,int therm_steps, int N, int draw_interval, bool hot_start, int multihit) {
    //------------------------------
    //            setup
    //------------------------------
    vector<int> config = {};
    if (hot_start) {
        config = Simulation::initializeLatticeHot(L);

    }
    else {
        config = Simulation::initializeLatticeCold(L);
    }
    
   
    for (int i = 0; i < therm_steps; i++) //thermalisation
    {
        Simulation::sweepMetropolisMultihit(config, L, beta, h, multihit);
    }
    
    double M = Simulation::averageMagnetisation(L * L, config);
    cout << std::format("thermalisation for beta = {:.3f} done, checking value after thermalisation of magnetisation: {:.5f}", beta, (double) M / (L * L)); //cecking magnetisation after thermalisation

    //initialise Observables, 2 are square values
    vector<double> energies1(N, 0);
    vector<double> energies2(N, 0);
    vector<double> absmag1(N, 0);
    vector<double> absmag2(N, 0);
    vector<double> mag1(N, 0);
    vector<double> mag2(N, 0);
    //------------------------------
    //            sweeps
    //------------------------------
    int K = L * L;
    int counter = 0;
    for (int i = 0; i < N; i++)
    {
        Simulation::draw(config, L ,beta , h , draw_interval, multihit);
        double E = Simulation::averageEnergy(config, K, L, h);
        energies1[i] = E;
        energies2[i] = E * E;
        double M = Simulation::averageMagnetisation(L * L, config);
        double absmag_i = abs(M) / K;
        absmag1[i] = absmag_i;
        absmag2[i] = absmag_i * absmag_i;
        mag1[i] = M / K;
        mag2[i] = M / K * M / K;
        
        cout << counter << endl;

        counter++;
    }
   
    //------------------------------
    //    calculate observables
    //------------------------------
    vector<double> vals = {};

    double e_mean = accumulate(energies1.begin(), energies1.end(), 0.0) / N; //accumulate(start, end) is a sum over all members
    double e2_mean = accumulate(energies2.begin(), energies2.end(), 0.0) / N;
    double m_mean = accumulate(mag1.begin(), mag1.end(), 0.0) / N;
    double m2_mean = accumulate(mag2.begin(), mag2.end(), 0.0) / N;
    double absm_mean = accumulate(absmag1.begin(), absmag1.end(), 0.0) / N;
    double absm2_mean = accumulate(absmag2.begin(), absmag2.end(), 0.0) / N;
    vals.push_back(e_mean);
    vals.push_back(sqrt((e2_mean - e_mean * e_mean) / (N - 1)));
    vals.push_back(m_mean);
    vals.push_back(sqrt((m2_mean - m_mean * m_mean) / (N - 1)));
    vals.push_back(absm_mean);
    vals.push_back(sqrt((absm2_mean - absm_mean * absm_mean) / (N - 1)));
    vals.push_back(beta * beta * (e2_mean - e_mean * e_mean)); //specific heat
    vals.push_back(e2_mean);
    vals.push_back(m2_mean);
    return vals;

   
}


//starts a single measurement (h, beta) with heatbath
static vector<double> startSimulationHeatbath(int L, double beta, double h, int therm_steps, int N, int draw_interval, bool hot_start) {
    cout << std::format("starting sim for beta = {:.5f}, h = {:.5f}", beta, h) << endl;

    //------------------------------
    //            setup
    //------------------------------
    int K = L * L;
    vector<int> config = {};
    if (hot_start) {
        config = Simulation::initializeLatticeHot(L);

    }
    else {
        config = Simulation::initializeLatticeCold(L);
    }


    for (int i = 0; i < therm_steps; i++) //thermalisation
    {
        Simulation::sweepHeatbath(config, L, K, beta, h);
    }

    //cecking magnetisation after thermalisation, 
    double M = Simulation::averageMagnetisation(K, config); 
    cout << std::format("thermalisation for beta = {:.5f} done, checking value after thermalisation of magnetisation: {:.5f}", beta, (double)M / K) << endl; 

    //initliase Observables, 2 are the square values
    vector<double> energies1(N, 0);
    vector<double> energies2(N, 0);
    vector<double> absmag1(N, 0);
    vector<double> absmag2(N, 0);
    vector<double> mag1(N, 0);
    vector<double> mag2(N, 0);
    //------------------------------
    //            sweeps
    //------------------------------
   
    
    for (int i = 0; i < N; i++)
    {
        Simulation::draw(config, L, beta, h, draw_interval);
        double E = Simulation::averageEnergy(config, K, L, h);
        energies1[i] = E;
        energies2[i] = E * E;
        double M = Simulation::averageMagnetisation(L * L, config);
        double absmag_i = abs(M) / K; //absulate magnetisation per lattic point
        absmag1[i] = absmag_i;
        absmag2[i] = absmag_i * absmag_i;
        mag1[i] = M / K;
        mag2[i] = M / K * M / K;
        

    }
   

    //------------------------------
    //    calculate observables
    //------------------------------
    vector<double> vals = {};

    double e_mean = accumulate(energies1.begin(), energies1.end(), 0.0) / N; //accumulate(start, end) is a sum over all members
    double e2_mean = accumulate(energies2.begin(), energies2.end(), 0.0) / N;
    double m_mean = accumulate(mag1.begin(), mag1.end(), 0.0) / N;
    double m2_mean = accumulate(mag2.begin(), mag2.end(), 0.0) / N;
    double absm_mean = accumulate(absmag1.begin(), absmag1.end(), 0.0) / N;
    double absm2_mean = accumulate(absmag2.begin(), absmag2.end(), 0.0) / N;
    vals.push_back(e_mean);
    vals.push_back( sqrt( (e2_mean - e_mean * e_mean) / (N - 1) ));
    vals.push_back(m_mean);
    vals.push_back(sqrt( (m2_mean - m_mean * m_mean) / (N - 1) ));
    vals.push_back(absm_mean);
    vals.push_back(sqrt( (absm2_mean - absm_mean * absm_mean) / (N - 1) ));
    vals.push_back(beta * beta * (e2_mean - e_mean * e_mean)); //specific heat
    vals.push_back(e2_mean);
    vals.push_back(m2_mean);
    return vals;

}

//testing verions of startSimulation(..) for optimal parameters
static void parameterTesting(string output_filename, int L, double beta, double h, int N, int draw_interval, bool hot_start, bool Metropolis, int multi_hit, int block_size, double epsilon, int therm_time ) {
    //------------------------------
    //            setup
    //------------------------------
    vector<int> config = {};
    if (hot_start) { //initialise 
        config = Simulation::initializeLatticeHot(L);

    }
    else {
        config = Simulation::initializeLatticeCold(L);
    }


    for (int i = 0; i < therm_time; i++)//thermalisation steps
    {
        cout << "i" << endl;
        Simulation::draw(config, L, beta, h, draw_interval, multi_hit);
    }


    ofstream File(output_filename);
    File << fixed << setprecision(5);
    File << "beta" << "\t" << "e" << "\t" << "m" << "\t" << "|m|" << "\n";
    //------------------------------
    //            sweeps
    //------------------------------
    int K = L * L; //no. lattice points

    if (Metropolis) { //if metropolis is algorythm of choice
        //initialise block means for energy and absolute value of magnetisation
        double prev_e_mean = 0;
        double prev_abs_m_mean = 0;
        clock_t start_2 = clock();

        for (int i = 0; i < N; i++)// N x block_size is the total amount of sweeps
        {
            
            
            vector<double> e_block(block_size, 0); //initialise block vectors
            vector <double> abs_m_block(block_size, 0);
            

            for (int j = 0; j < block_size; j++) { //block

                Simulation::draw(config, L, beta, h, draw_interval, multi_hit);
                double e = Simulation::averageEnergy(config, K, L, h); 
               
                double M = Simulation::averageMagnetisation(L * L, config);
                double absmag_i = abs(M) / K;
                

                double m = M / K;
                
                e_block[j] = e;
                abs_m_block[j] = absmag_i;

                //writing the important data for testing thermalisation to a file
                File << beta << "\t" << e << "\t" << m << "\t" << absmag_i << "\n";


            }

            double next_e_mean = accumulate(e_block.begin(), e_block.end(), 0.0) / block_size;
            double next_abs_m_mean = accumulate(abs_m_block.begin(), abs_m_block.end(), 0.0) / block_size;
            
           
            if (i > 0) {
                cout << "comparing " << next_abs_m_mean << " with " << prev_abs_m_mean << " : " << abs(next_abs_m_mean - prev_abs_m_mean) << endl;
                clock_t end_2 = clock();
                double elapsed = double(end_2 - start_2) / CLOCKS_PER_SEC;
                cout << elapsed << endl;
                if (abs(next_abs_m_mean - prev_abs_m_mean) < epsilon && abs(next_e_mean - prev_e_mean) < epsilon) {
                    cout << "convergence reached" << to_string(i * block_size * draw_interval) << endl;
                   
                }

            }
            
            prev_abs_m_mean = next_abs_m_mean;
            prev_e_mean = next_e_mean;

        }
       
    }
    else
    {
        cout << "starting heatbath" << endl;
        //initialise block means for energy and absolute value of magnetisation
        double prev_e_mean = 0;
        double prev_abs_m_mean = 0;
        clock_t start_2 = clock();

        for (int i = 0; i < N; i++)// N x block_size is the total amount of sweeps
        {


            vector<double> e_block(block_size, 0); //initialise block vectors
            vector <double> abs_m_block(block_size, 0);


            for (int j = 0; j < block_size; j++) { //block

                Simulation::draw(config, L, beta, h, draw_interval);
                double e = Simulation::averageEnergy(config, K, L, h);

                double M = Simulation::averageMagnetisation(L * L, config);
                double absmag_i = abs(M) / K;


                double m = M / K;

                e_block[j] = e;
                abs_m_block[j] = absmag_i;

                //writing the important data for testing thermalisation to a file
                File << beta << "\t" << e << "\t" << m << "\t" << absmag_i << "\n";


            }

            double next_e_mean = accumulate(e_block.begin(), e_block.end(), 0.0) / block_size;
            double next_abs_m_mean = accumulate(abs_m_block.begin(), abs_m_block.end(), 0.0) / block_size;


            if (i > 0) {
                cout << "comparing " << next_abs_m_mean << " with " << prev_abs_m_mean << " : " << abs(next_abs_m_mean - prev_abs_m_mean) << endl;
                clock_t end_2 = clock();
                double elapsed = double(end_2 - start_2) / CLOCKS_PER_SEC;
                cout << elapsed << endl;
                if (abs(next_abs_m_mean - prev_abs_m_mean) < epsilon && abs(next_e_mean - prev_e_mean) < epsilon) {
                    cout << "convergence reached" << to_string(i * block_size) << endl;

                }

            }

            prev_abs_m_mean = next_abs_m_mean;
            prev_e_mean = next_e_mean;

        }
    }

    File.close();
    return;

}

//simulating values for 3D-Phaseplot
static void excercise4C() {
    //std::vector<double> h_vals = { -0.2, -0.1, -0.05, -0.01, 0 ,0.01,0.05, 0.1, 0.2 };
    //std::vector<double> h_vals = { -0.02, -0.005, -0.002, -0.001, 0.001, 0.002, 0.005, 0.02 };
    //std::vector<double> betas = { 0.1, 0.2, 0.3, 0.4, 0.42, 0.43, 0.44, 0.46, 0.5, 0.6, 0.8 };
    //std::vector<double> betas = {0.442, 0.444, 0.446, 0.448, 0.45, 0.48};
   
    std::vector<double> betas = { 0.8, 0.6};
    std::vector<double> h_vals = { -0.05, -0.01};
    



    //std::vector<int>  therm_steps = { 20000, 20000, 20000, 500, 10000, 10000, 10000, 10000, 500, 500, 500 };
    int therm_steps = 30000;
    

    bool hot = true;
    int draw_interval = 100;
    int N = 50;
    int L = 32;
    //string output_files = "output_files/3DPhasediagramm.txt";
    //string output_files = "output_files/3DPhasediagramm2.txt";
    string output_files = "output_files/3DPhasediagramm3.txt";


    ofstream File(output_files);
    File << "beta" << "\t" << "h" << "\t" << "<m>" << "\t" << "dm" << "\t" << "<|m|>" << "\t" << "d|m|" << "\n";
    File << fixed << setprecision(12);
    //int i = 0;
    for (int i = 0; i < betas.size(); i++)
    {
        mt19937 gen(random_device());
        //vector<double> vals = startSimulationHeatbath(L, beta, h, therm_steps[i], N, draw_interval, hot);
        vector<double> vals = startSimulationHeatbath(L, betas[i], h_vals[i], therm_steps, N, draw_interval, hot);
        File << betas[i] << "\t" << h_vals[i] << "\t" << vals[2] << "\t" << vals[3] << "\t" << vals[4] << "\t" << vals[5] << "\n";
    }

    /*
    for (double beta : betas)
    {
        for (double h : h_vals) {
            mt19937 gen(random_device());
            //vector<double> vals = startSimulationHeatbath(L, beta, h, therm_steps[i], N, draw_interval, hot);
            vector<double> vals = startSimulationHeatbath(L, beta, h, therm_steps, N, draw_interval, hot);
            File << beta << "\t" << h << "\t" << vals[2] << "\t" << vals[3] << "\t" << vals[4] << "\t" << vals[5] << "\n";
        }
        i++;
    }
    */

    File.close();
    //TODO: testit with small N and therm steps

}


//checking hysteris
static void exercise4B() {
    int therm_time = 10000;
    double beta = 0.5; //beta in ferromagnetic phase
    int L = 32;
    int K = L * L;
    double h = -0.5;
    //double h_start = -0.1;
    vector<int> config = {};
    config = Simulation::initializeLatticeCold(L);

    for (int i = 0; i < therm_time; i++)
    {
        Simulation::sweepHeatbath(config, L, K, beta, h);
    }
    double M = Simulation::averageMagnetisation(K, config);
    cout << std::format("thermalisation for beta = {:.5f} done, checking value after thermalisation of magnetisation: {:.5f}", beta, (double)M / K); //cecking magnetisation after thermalisation



    string output_files = "output_files/Hysterese_neg_3.txt";
    ofstream File(output_files);
    File << h << "\t" << "<m>" << "\n";
    File << fixed << setprecision(12);
    M = Simulation::averageMagnetisation(L * L, config);
    double m = M / K;
    File << h << "\t" << m << "\n";

    for (int i = 0; i < 10; i++)
    {
        h += 0.05;
        for (int k = 0; k < 200; k++)
        {
            for (int j = 0; j < 10; j++) {
                Simulation::sweepHeatbath(config, L, K, beta, h);
                
            }
            M = Simulation::averageMagnetisation(L * L, config);
            m = M / K;
            File << h << "\t" << m << "\n";
        }
        
        M = Simulation::averageMagnetisation(L * L, config);
        m = M / K;
        File << h << "\t" << m << "\n";
    }
    File.close();


    //Observables
    int N = 50;
    int draw_interval = 200;
    vector<double> mag1(N, 0);
    vector<double> mag2(N, 0);
    
    //------------------------------
    //            sweeps
    //------------------------------

    int counter = 0;
    for (int i = 0; i < N; i++)
    {
        Simulation::draw(config, L, beta, h, draw_interval);
        double E = Simulation::averageEnergy(config, K, L, h);
        
        double M = Simulation::averageMagnetisation(L * L, config);
        
        mag1[i] = M / K;
        mag2[i] = M / K * M / K;
        //cout << counter << endl;
        counter++;

    }
    double m_mean = accumulate(mag1.begin(), mag1.end(), 0.0) / N;
    double m2_mean = accumulate(mag2.begin(), mag2.end(), 0.0) / N;
    cout << "\n" << m_mean << endl;
    cout << sqrt((m2_mean - m_mean * m_mean) / (N - 1)) << endl;
}

//exercises 3a and 4a to
static void single_measurements() {
    double h = 0.0; //external magnetic field
    
    int therm_steps = 20000; //number of thermalize sweeps
    //number of draws
    int N = 50; //actual number of sweeps is draw_interval * N
    int draw_interval = 200; //sweeps between drawing
    //lattice size
    int L = 128; //actual size is LxL

    //vector<double> betas = {0.4406868};
    vector<double> betas = { 0.8 }; //each beta also needs a seed, if this is not wanted switch below in the for loop
    vector<double> seeds = { 8731 };
    int multi_hit = 1;
    bool hot = true;
    
    int counter = 0;
    //string algorithm = "Metropolis";
    string algorithm = "Heatbath";

    for (double beta : betas)
    {
        Simulation::gen = mt19937(seeds[counter]); //switch to below if not wanting to provide seeds
        //Simulation::gen = mt19937(random_device()); 

        string path = "output_files/heatbath_output/";
        string filename = to_string(L) + "_" + algorithm + "_beta = " + to_string(beta) + "_start = " + to_string(hot) + ".txt";
        //vector<double> vals = startSimulationMetropolis(L, beta, h, therm_steps, N, draw_interval, hot, multi_hit);
        vector<double> vals = startSimulationHeatbath(L, beta, h, therm_steps, N, draw_interval, hot);
        //Writing to file, each simulation gets its own seperate file, basically everything is written there (parameters and observables) expcept hor h
        ofstream File(path + filename);
        File << fixed << setprecision(12);
        File << "beta" << "\t" << "<e>" << "\t" << "de" << "\t" << "<m>" << "\t" << "dm" << "\t";
        File << "<|m|>" << "\t" << "d|m|" << "\t" << "c_v/(L^2)" << "\t" << "<e^2>" << "\t" << "<m^2>" << "\t";
        File << "t_therm" << "\t" << "draw_interval" << "\t" << "draws" << "\t" << "seed" << "\n";
        File << beta << "\t" << vals[0] << "\t" << vals[1] << "\t" << vals[2] << "\t" << vals[3] << "\t" << vals[4] << "\t";
        File << vals[5] << "\t" << vals[6] << "\t" << vals[7] << "\t" << vals[8] << "\t";
        File << therm_steps << "\t" << draw_interval << "\t" << N << "\t" << seeds[counter] << "\n";
        File.close();


        counter++;

    }
    

    



}

//testing thermalisation, multi-hit etc.
static void testing() {
    
    double h = 0.0;//external magnetic field

    int therm_steps = 0;//number of thermalize sweeps

    int draw_interval = 1;//sweeps between drawing

    int L = 128; //lattice dimension -> actual size is LxL

    int block_size = 100; //averages over <block_size> elements to test for convergence of the average
    double epsilon = 0.001; //convergence paramter, a lot times helpful but not used to stop simulation

    int N = 20000; //actual number of sweeps is draw_interval * N * blocksize

    int mh = 1; //multihit parameter
    vector<double> betas = { 0.44 };
    
    //vector<double> betas = { 0.1, 0.2, 0.3,0.35 }; 
    
    //string output_path = "output_files/testing_therm_params/testing_heatbath_L=" + std::to_string(L) + "_";
    string output_path = "output_files/testing_therm_params/testing_metro_L=" + std::to_string(L) + "_";
    
    

    bool hot = false;//hot or coldstart


    for (double beta : betas) {
        parameterTesting(output_path + std::format("beta={:.3f}_MH={}_start={}.txt", beta, mh, hot), L, beta, h, N, draw_interval, hot, true, mh, block_size, epsilon, therm_steps);
    }
    cout << (double)Simulation::acceptance_count / (double)Simulation::try_count << endl; //print acceptance rate for multi-hit, not really used it
}



int main()
{
    clock_t start = clock();
    //--------------------------------------------
    //                 exercise 1
    //--------------------------------------------
    /*
    cout << std::format("approximation of pi with MC : {:.4f}", approximationOfPi(1000, 4242)) << endl;
    cout << std::format("approximation of pi with MC : {:.4f}", approximationOfPi(10000000, 2424)) << endl;
    cout << std::format("approximation of gaussian integral : {:.4f}", approxGaussianIntegral(100000)) << endl;
    */
    //--------------------------------------------
    //                 exercise 2
    //--------------------------------------------
    /*
    explicitIsing(2);
    explicitIsing(3);
    explicitIsing(4);
    */
    //--------------------------------------------
    //          exercises 3 & 4
    //--------------------------------------------
    //--------------------------------------------
    //                 1.testing
    //--------------------------------------------

   
    testing();
    
    //--------------------------------------------
    //               2.measurements
    //--------------------------------------------
    //single_measurements();
    //excercise4C();
    //exercise4B();

    
    


    
    

   
    clock_t end = clock();
    double elapsed = double(end - start) / CLOCKS_PER_SEC;
    printf("execution time: %.3f sec", elapsed);
    
}

