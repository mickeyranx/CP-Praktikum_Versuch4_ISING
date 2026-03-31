#include <iostream>
#include <vector>
#include "Simulation.h"
#include <tuple>
using namespace std;

//initialise tools for randomization
uniform_real_distribution<double> Simulation::distr1 = uniform_real_distribution<double>(0.0, 1.0);
uniform_int_distribution<int> Simulation::distr2 = uniform_int_distribution<int>(0, 1);
random_device Simulation::rand_device;
mt19937 Simulation::gen = mt19937(rand_device());
//default_random_engine Simulation::gen2 = default_random_engine();

//variables for testing optimal multihit parameter
unsigned int Simulation::try_count = 0;
unsigned int Simulation::acceptance_count = 0;




vector<int> Simulation::initializeLatticeCold(int L) {
    vector<int> lattice(pow(L, 2), 1);
    return lattice;
}

vector<int> Simulation::initializeLatticeHot(int L)
{
    vector<int> lattice = {};
    for (int i = 0; i < pow(L, 2); i++)
    {
        int r = distr2(gen); //generate random number between 0 and 1
        //fill latice with random -1 or 1
        lattice.push_back((r == 0) ? -1 : 1);
    }
    return lattice;
}


vector<int> Simulation::getNeighborPos(int i, int L) {
    //determine which row we are in
    double row = floor(double(i) / L);
    vector<int> neigbor_pos;


    //check top neighbor
    if (i - L < 0) {
        //add L^2 - (L - i) 
        neigbor_pos.push_back(pow(L, 2) - (L - i));
    }
    else {
        //add i - L
        neigbor_pos.push_back(i - L);
    }

    //check right neighbor
    if (i + 1 > ((row + 1) * L) - 1) {
        //add i - L - 1

        neigbor_pos.push_back(i - L + 1);
    }
    else {
        //add i+1

        neigbor_pos.push_back(i + 1);
    }

    //check bottom neighbor
    if (i + L > pow(L, 2) - 1) {
        //add i - (L*(L-1))

        neigbor_pos.push_back(i - (L * (L - 1)));
    }
    else {
        //add i+L

        neigbor_pos.push_back(i + L);
    }

    //check left neigbor
    if (i - 1 < ((row + 1) - 1) * L) {
        //i + (L - 1)

        neigbor_pos.push_back(i + (L - 1));
    }
    else {
        //add i-1

        neigbor_pos.push_back(i - 1);
    }
   
    return neigbor_pos;

}



double Simulation::changeInEnergy(vector<int> &config_1,int L ,int pos, int s, double h) {
    vector<int> neighbors = getNeighborPos(pos, L);

    
    double dH = 2.0 * s * (config_1[neighbors[0]] + config_1[neighbors[1]] + config_1[neighbors[2]] + config_1[neighbors[3]] + h);
    return dH;
}



void Simulation::sweepMetropolisMultihit(vector<int> &config_1, int L ,double beta, double h, int tries) {
    for (int i = 0; i < L*L; i++)
    {
        
        
        for (int j = 0; j < tries; j++)
        {
            int s_i = config_1[i];
            try_count++;
           
                       
            int s_i_new = (distr2(gen) == 0) ? -1 : 1; //select new random spin -1 or 1 

           
            double dH = 0;
            //calculate change in energy
            if (s_i != s_i_new) { //in case of spinflip
                dH = changeInEnergy(config_1, L,i, s_i, h);
                //std::cout << dH << std::endl;
                if (dH < 0) {
                    config_1[i] = s_i_new; //accept immediatly if dH < 0
                    Simulation::acceptance_count++;
                    
                }
              
                if (distr1(gen) < exp(-beta * dH)) {
                    
                    config_1[i] = s_i_new;
                    Simulation::acceptance_count++;
                    
                }
            }
            
     
        }


    }


}


void Simulation::sweepHeatbath(vector<int> &config_1, int L,int M ,double beta, double h) {
    for (int i = 0; i < M; i++)
    {
        vector<int> neighbors = getNeighborPos(i, L);
        double k = beta * ((double)(config_1[neighbors[0]] + config_1[neighbors[1]] + config_1[neighbors[2]] + config_1[neighbors[3]]) + h);
        (distr1(gen) < (1.0 / (1.0 + exp(-2.0 * k)))) ? config_1[i] = 1 : config_1[i] = -1; 
        
    }

}

void Simulation::draw(vector<int> &config, int L ,double beta, double h, int draw_interval) {
    int K = L * L;
    for (int i = 0; i < draw_interval; i++)
    {
        sweepHeatbath(config, L , K ,beta, h);
    }
   
}

void Simulation::draw(vector<int>& config, int L, double beta, double h, int draw_interval, int multi_hit) {
    for (int i = 0; i < draw_interval; i++)
    {
        sweepMetropolisMultihit(config, L, beta, h, multi_hit);
        
    }

}


double Simulation::averageEnergy(vector<int>& config, int K, int L ,double h) {
   
    double sum = 0;
    for (int i = 0; i < K; i++) //iterate over all lattice points
    {
        int s = config[i];
        vector<int> neigbors = getNeighborPos(i, L);
        //calculate Hamiltonian
        sum += (double) -s * (config[neigbors[0]] + config[neigbors[1]]) - (double)h * s;//only look at 2 neighbors per orientation to avoid double counting

    }
    return sum/K; //return averaged energy
}


double Simulation::averageMagnetisation(int M, vector<int> &config) {
    double sum = 0;
    for (int i = 0; i < M; i++)
    {
        sum += config[i];

    }
    return  sum;

}




//for debugging
void Simulation::printConfig(vector<int> &config) {
    int row = 0;
    int L_square = config.size();
    for (int i = 0; i < L_square; i++)
    {
        int current_spin = config[i];
        if (floor(double(i) / sqrt(L_square)) > row) {
            row++;
            cout << "\n";
        }

        if(current_spin == -1){
            cout << current_spin << " ";
        }
        else {
            cout << " " << current_spin << " ";
        }

            
   
    }

}








