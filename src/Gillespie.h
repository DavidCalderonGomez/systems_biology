#pragma once
#include <cmath>
#include "Random64.h"

//Asked Claude typical values for lacZ
const double kr=2, gammar=0.35, kp=10, gammap=0.023; 
const double alpha=0.001, beta=1, K=10, h=2;  

//Time step for the naive method
const double dt=0.1;

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

    //método inocente de tiempos fijos
    void naive_step_core(Crandom &ran64, double a1); //naive method regardless of feedback
    void naive_step(Crandom &ran64); //naive method without feedback
    void naive_step(Crandom &ran64, bool repression); //naive method with feedback
    
    //Gillespie
    void Gillespie_step_core(Crandom &ran64, double a1); //method for Gillespie step regardless of feedback
    void Gillespie_step(Crandom &ran64); //method for Gillespie step without feedback
    void Gillespie_step(Crandom &ran64, bool repression); //method for Gillespie step with feedback

    //protein bursts (do not simulate each protein event, just the arn event that creates a bunch of proteins super fast )
    void Gillespie_burst_core(Crandom &ran64, double a1);   // burst step regardless of feedback
    void Gillespie_burst(Crandom &ran64);                   // no feedback
    void Gillespie_burst(Crandom &ran64, bool repression);  // feedback
};

//Hill function possible to use for repression or activation
double Cell::hill_function(double protein, bool repression){ 
    if(repression){
        return alpha + beta*pow(K,h)/(pow(K,h)+pow(protein,h));
    }
    return alpha + beta*pow(protein,h)/(pow(K,h)+pow(protein,h));
}

void Cell::naive_step_core(Crandom &ran64, double a1){
    //Calculate the propensities
    double a2=gammar*mRNA; //mRNA degradation
    double a3=kp*mRNA; //protein creation
    double a4=gammap*protein; //protein degradation
    double a0=a1+a2+a3+a4; //total propensity
    //generate a random number uniformly distributed between 0 and 1
    if(ran64.r()<a1*dt) mRNA+=1;
    if(ran64.r()<a2*dt) mRNA-=1;
    if(ran64.r()<a3*dt) protein+=1;
    if(ran64.r()<a4*dt) protein-=1;
    time+=dt;
}

void Cell::naive_step(Crandom &ran64){
    naive_step_core(ran64, kr);
}

void Cell::naive_step(Crandom &ran64, bool repression){
    naive_step_core(ran64, kr*hill_function(protein, repression));
}

void Cell::Gillespie_step_core(Crandom &ran64, double a1){

        //Calculate the propensities
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

void Cell::Gillespie_step(Crandom &ran64){
    Gillespie_step_core(ran64, kr);
}

// Sobrecarga: solo calcula a1 con Hill
void Cell::Gillespie_step(Crandom &ran64, bool repression){
    Gillespie_step_core(ran64, kr*hill_function(protein, repression));
}

// en el .cpp
void Cell::Gillespie_burst_core(Crandom &ran64, double a1){
    double b=kp/gammar;          // proteínas por ARN durante su vida
    double a2=gammap*protein;    // degradación de proteína
    double a0=a1+a2;
    time+=ran64.exponencial(1/a0);
    if(ran64.r()*a0<a1) protein+=ran64.poisson(b);//here we sum a poisson variable instead of a constant b each burst. 
    else protein-=1;
}

void Cell::Gillespie_burst(Crandom &ran64){
    Gillespie_burst_core(ran64, kr);
}

void Cell::Gillespie_burst(Crandom &ran64, bool repression){
    Gillespie_burst_core(ran64, kr*hill_function(protein, repression));
}