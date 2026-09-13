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
    for(int i=0; i<steps; i++){
        cell.Gillespie_step(ran64, repression); //perform a Gillespie step
        //Output the current state of the cell
        //time, mRNA, protein a single line separated by tabs
        cout << cell.get_time() << "\t" << cell.get_mRNA() << "\t" << cell.get_protein() << "\t";    
    }
return 0;
}

