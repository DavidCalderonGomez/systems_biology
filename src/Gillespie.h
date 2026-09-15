#pragma once
#include <cmath>
#include "Random64.h"

const double kr=1, gammar=0.2, kp=60, gammap=0.03; //creation and degradation rates for mRNA and protein
const double alpha=0.1, beta=1, K=4000, h=2; //parameters for the Hill function
const int steps=250000; //number of steps for the simulation

//Class for the cell with mRNA and protein
class Cell{
    double mRNA, protein, time;//attributes for mRNA and protein at a time t
    public://methods for the cell
    //Constructor
    Cell(double m, double p, double t): mRNA(m), protein(p), time(t) {}
    //Getters and setters
    double get_mRNA(){return mRNA;}
    double get_protein(){return protein;}
    double get_time(){return time;}
    void set_mRNA(double m){mRNA=m;}
    void set_protein(double p){protein=p;}
    void set_time(double t){time=t;}
    //Hill function possible to use for repression or activation
    double hill_function(double protein, bool repression);
    //Gillespie algorithm for the cell
    void Gillespie_step(Crandom &ran64, bool repression);
};

//Hill function possible to use for repression or activation
double Cell::hill_function(double protein, bool repression){ 
    if(repression){
        return alpha + beta*pow(protein,h)/(pow(K,h)+pow(protein,h));
    }
    return alpha + beta*pow(K,h)/(pow(K,h)+pow(protein,h));
}

void Cell::Gillespie_step(Crandom &ran64, bool repression){

        //Calculate the propensities
        double a1=kr*hill_function(protein, repression); //mRNA creation
        double a2=gammar*mRNA; //mRNA degradation
        double a3=kp*mRNA; //protein creation
        double a4=gammap*protein; //protein degradation
        double a0=a1+a2+a3+a4; //total propensity
        //generate a random number uniformly distributed between 0 and 1
        double r=ran64.r();
        //update the time using the exponential distribution
        time+=ran64.exponencial(1/a0); 
        //Determine which reaction occurs
        if(r*a0<a1){
            mRNA+=1; //mRNA creation
        }else if(r*a0<a1+a2){
            mRNA-=1; //mRNA degradation
        }else if(r*a0<a1+a2+a3){
            protein+=1; //protein creation
        }else{
            protein-=1; //protein degradation
        }
}


