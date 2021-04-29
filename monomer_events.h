#include "monomer.h"

#ifndef MONO_EVENTS_H
#define MONO_EVENTS_H

void monomer::polymerize(monomer *mono2)
{
     //set the pointers and bools
     //for simulations with active polymerization/depolymerization, rates can be added to object and adjusted in this function
        plus = mono2;//this mono's plus points to monomer2
        plus_is_polymerized = true;
        (*mono2).minus = this;//Mono2's minus points to this mono
        (*mono2).minus_is_polymerized = true;

}//end of polymerize()



int monomer::depolymerize()
{ 
//function for depolymerizing. not used in current version of code (April 29, 2021)  
double rrr = ranf5();//unif_rand();

   if(rrr < depolymerize_prob)
   { /*note that this is okay, because depolymerize_prob = 0 unless polymerization has occured*/
       //Change minus' plus pointer & bool, wait until end to change this' minus ptr 
       (*minus).plus = NULL;
       (*minus).plus_is_polymerized = false;
       minus_is_polymerized = false;
       depolymerize_prob = 0.;

       //Reset minus pointer,bool
       minus = NULL;       
         
       return 1;
   }
   return 0;
}//end of depolymerize


#endif
