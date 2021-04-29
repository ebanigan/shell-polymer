#include <stdlib.h>
#include <sys/types.h>
#include <fcntl.h>
//#include <sys/uio.h>
#include <unistd.h>
#include <iostream>
#include <vector>

#include "common.h"
#include "protos.h"

using namespace std;

#ifndef MONO_H
#define MONO_H


/*The monomer is the basic building block.  It can connect and disconnect from other monomers (polymerize/depolymerize)*/

class monomer{
   protected:
       double depolymerize_prob;//probability to depolymerize 

       long id; // identification for monomers

   public:
      //functions to get values of protected variables
      long get_id(){return id;}      
      double get_depolymerize_prob(){return depolymerize_prob;}
      void set_depolymerize_prob(double prob){depolymerize_prob = prob;}

      monomer *plus;  //plus member of a monomer points to another monomer or NULL
      monomer *minus; //minus member of monomer1 points to monomer2, whose plus points to monomer1

      //boolean flags to indicate states of monomer
      bool plus_is_polymerized;
      bool minus_is_polymerized;
      bool crosslinked;
      bool shell_poly_link;

      bool get_crosslinked(){return crosslinked;}
      void set_crosslinked(bool input){crosslinked = input;}

      bool get_shell_poly_link(){return shell_poly_link;}
      void set_shell_poly_link(bool input){shell_poly_link = input;}

      virtual ~monomer(){};
      //default constructor
      monomer(){}
      //constructor for a free monomer - sets id #, ptrs, flags, and initial rate constants
      monomer(long num)
      { 
        id = num;
        depolymerize_prob = 0.;
       
        plus = NULL;
        minus = NULL;
        plus_is_polymerized = false;
        minus_is_polymerized = false;
      }




      //Returns the monomer id that the plus points to. Returns -1 for pointing to NULL.
      long plus_points_to()
      {
          if(plus_is_polymerized)
           return (*plus).get_id();
          else
           return -1;
      }
      long minus_points_to()
      {
          if(minus_is_polymerized)
           return (*minus).get_id();
          else
           return -1;
      }

      //if monomer is not attached to anything, returns true.
      bool free_monomer()
      {
           return (!(plus_is_polymerized || minus_is_polymerized));
      }
      
      
/*Connects ("polymerizes")monomers (used to create dimers, polymers).
Also makes appropriate changes to pointers and rates.*/
      void polymerize(monomer *mono2);
  
/*If rand# < depol prob, this function is cuts the connection between two monomers 
by setting the appropriate plus/minus pointers to NULL. add adjustments to other rates here if needed.*/
      int depolymerize();



};


#endif
