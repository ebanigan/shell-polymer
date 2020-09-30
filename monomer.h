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


/*The monomer is the basic building block.  It can connect and disconnect from other monomers
(nucleate, branch, denucleate, debranch). It can be capped (prevented from nucleating) and uncapped.
It can be tagged (associated) with arp23 or detagged (dissociated).*/

class monomer{
   protected:
       double nucleate_prob; //probability to nucleate plus -- adjustable (for ex. set = 0 if monomer capped)
       //double polym_prob; //possibly a future feature to separate nucleation and polymerization processes
       double denucleate_prob;//probability to denucleate - in the future may add depolym_prob too
       double new_denucleate_prob;
       double new_debranching_prob;
       double capping_prob; //probability to cap, i.e, prevent nucleation a monomer
       double uncapping_prob; //prob to uncap
       double associate_arp23_prob, dissociate_arp23_prob; //prob to tag or untag with arp23
       double branching_prob; //once mono has arp23, it has prob. to branch after nucleated in a filament.
       double debranching_prob; //probability for branch to fall off

       long id; // identification for monomers

   public:
   //functions to get values of protected variables
      long get_id(){return id;}      
      double get_nucleate_prob(){return nucleate_prob;}
      double get_denucleate_prob(){return denucleate_prob;}
      double get_new_denucleate_prob(){return new_denucleate_prob;}
      double get_new_debranching_prob(){return new_debranching_prob;}
      double get_capping_prob(){return capping_prob;}
      double get_uncapping_prob(){return uncapping_prob;}
      double get_associate_arp23_prob(){return associate_arp23_prob;}
      double get_dissociate_arp23_prob(){return dissociate_arp23_prob;}
      double get_branching_prob(){return branching_prob;}
      double get_debranching_prob(){return debranching_prob;}
      
      void set_associate_arp23_prob(double prob){associate_arp23_prob = prob;}
      void set_dissociate_arp23_prob(double prob){dissociate_arp23_prob = prob;}
      void set_branching_prob(double prob){branching_prob = prob;}
      void set_capping_prob(double prob){capping_prob = prob;}
      void set_uncapping_prob(double prob){capping_prob = prob;}
      void set_nucleate_prob(double prob){nucleate_prob = prob;}
      void set_denucleate_prob(double prob){denucleate_prob = prob;}
      void set_new_denucleate_prob(double prob){new_denucleate_prob = prob;}
      void set_new_debranching_prob(double prob){new_debranching_prob = prob;}
      void update_denucleate_prob(){denucleate_prob = new_denucleate_prob;}
      void update_debranching_prob(){debranching_prob = new_debranching_prob;}
      void set_debranching_prob(double prob){debranching_prob = prob;}

      monomer *plus;  //plus member of a monomer points to another monomer or NULL
      monomer *minus; //minus member of monomer1 points to monomer2, whose plus points to monomer1
                      //also, the first monomer on a new branch has has minus pointing to a monomer on old filament      
      monomer *new_branch; // points to first monomer in the new branch

//boolean flags to indicate states of monomer
      bool plus_is_nucleated;
      bool minus_is_nucleated;
      bool plus_is_capped;
      bool mono_has_arp23;
      bool mono_has_branch;
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
        #if INITIALIZE_DIMERS
        nucleate_prob = 0.;//want nucleate prob = 0.; when initiating dimers.
        #else
        nucleate_prob = NUCLEATE_PARAM;
        #endif
        denucleate_prob = 0.;
	new_denucleate_prob = 0.;
        new_debranching_prob = 0.;
        capping_prob = 0.;
        //capping_prob = CAPPING_PARAM;
        uncapping_prob = 0.;
        branching_prob = 0.;
        debranching_prob = 0.;
        associate_arp23_prob = 0.;
        //associate_arp23_prob = ASSOC_ARP23_PARAM;
        dissociate_arp23_prob = 0.;
       
        plus = NULL;
        minus = NULL;
        new_branch = NULL;
        plus_is_nucleated = false;
        minus_is_nucleated = false;
        plus_is_capped = false;
        mono_has_arp23 = false;
        mono_has_branch = false;
      }




//Returns the monomer id that the plus points to. Returns -1 for pointing to NULL.
      long plus_points_to()
      {
          if(plus_is_nucleated)
           return (*plus).get_id();
          else
           return -1;
      }
      long minus_points_to()
      {
          if(minus_is_nucleated)
           return (*minus).get_id();
          else
           return -1;
      }
      long new_branch_points_to()
      {
          if(mono_has_branch)
           return (*new_branch).get_id();
          else
           return -1;
      }

//if monomer is not attached to anything, returns true.
      bool free_monomer()
      {
           //return ((!plus_is_nucleated)&&(!minus_is_nucleated)&&(!mono_has_branch));
           return (!(plus_is_nucleated || minus_is_nucleated || mono_has_branch));
      }
      
/*If unif_rand# < nucleation prob, this function connects this monomer to second monomer
by making plus point to second monomer and making second monomer's minus
point to this monomer.  Also, this will make appropriate adjustments to nucleation, denucleation,
capping, etc. rates.*/
      int nucleate(monomer *mono2);
      
/*Connects ("polymerizes")monomers (used to create dimers, polymers, and in "nucleate" function).
Also makes appropriate changes to pointers and rates.*/
      void polymerize(monomer *mono2);
  
/*If rand# < denuc prob, this function is cuts the connection between two monomers 
by setting the appropriate plus/minus pointers to NULL.  Also, this will make appropriate
adjustments to nucleation, capping, etc. probabilities.*/
      int denucleate();

/*If rand# < capping prob, monomer will become capped (and nucleation prob set to 0).*/
      int cap();

/*Reverse of cap. If rand# < uncapping prob, monomer loses cap, and for instance can nucleate its plus end again.
Can call remove_cap()*/
      int uncap();
      
      void remove_cap();

/*If rand# < associate prob, monomer becomes tagged with arp23, and can be a base for branches
when it is or becomes part of a filament*/
      int associate_arp23();

/*If rand# < dissoc prob, monomer spontaneously loses arp23.*/
      int dissociate_arp23();
      
      void remove_arp23();

      void construct_branch(monomer *mono2);

/*If rand#< branch prob, new_branch will point to mono2 and minus of minus2 will point to this monomer.
Thus a branch is formed*/
      int branch(monomer *mono2);

/*Cuts off branches the same way denucleate cuts ties between other monomers.*/
      int debranch();      
      void detach_branch();
//debranch calls detach_branch based on random number & probability, detach_branch performs the action

};


#endif
