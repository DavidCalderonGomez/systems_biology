#include <iostream>
#include <cmath>
#include "Random64.h"
#include "Gillespie.h"

using std::cout, std::endl;

int main(){

    Crandom ran64(123456789); //initialize the random number generator with a seed
    Cell cell(0, 0, 0); //initialize the cell with 0 mRNA, 0 protein, and time 0
    bool repression = true; //set the repression flag to true for the Hill function

    //Run the Gillespie algorithm for a number of steps

    for(int k=0; k<50; k++){ //run 50 independent simulations

        //It is better to use a while loop and impose the condition for some final time, for example 1000, instead of a fixed number of steps.
        while(cell.get_time() < 500){ //run until time reaches 500
            cell.Gillespie_step(ran64, repression); //perform a Gillespie step
            //Output the current state of the cell
            //time, mRNA, protein a single line separated by tabs
            cout << cell.get_time() << "\t" << cell.get_mRNA() << "\t" << cell.get_protein() << "\t";    
        }
        cout << endl; //separate each simulation with a new line
        cell.set_mRNA(0); //reset mRNA to 0
        cell.set_protein(0); //reset protein to 0
        cell.set_time(0); //reset time to 0
    }
return 0;
}

