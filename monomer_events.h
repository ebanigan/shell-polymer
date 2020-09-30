#include "monomer.h"

#ifndef MONO_EVENTS_H
#define MONO_EVENTS_H

void monomer::polymerize(monomer *mono2)
{
     //set the pointers and bools
        plus = mono2;//this mono's plus points to monomer2
        plus_is_nucleated = true;
        (*mono2).minus = this;//Mono2's minus points to this mono
        (*mono2).minus_is_nucleated = true;

     //Change probabilities for this monomer
     //In particular, MUST have the next two lines
        nucleate_prob = 0.;//this mono can't nucleate anymore
        capping_prob = 0.;//plus can't be capped 
        #if (!INITIALIZE_DIMERS)
        if(mono_has_arp23)
         if(!mono_has_branch)
           branching_prob = BRANCHING_PARAM;//However, if it has arp2/3, and no branch, it can branch.      
        #endif//with INITIALIZE_DIMERS = true, this mono can already branch if it has arp2/3
        
     //Change probabilities for mono2
/*
//next 3 lines only necessary if associate_arp23_prob = assoicate_arp23_param for FREE monomers
        if((*mono2).mono_has_arp23)
         if(!((*mono2).mono_has_branch))//this line really isn't necessary...
           (*mono2).branching_prob = BRANCHING_PARAM;//This mono can branch too if it has arp23 and no branch
*/
           
        #if INITIALIZE_DIMERS
        //if((*mono2).!plus_is_nucleated)//shouldn't be necessary since I'm only nucleating individual monos
        (*mono2).nucleate_prob = NUCLEATE_PARAM;//For INITIALIZE_DIMERS case, we only want barbs to polym.
                                                //no if statements needed since only free monos are added to barb
        #endif
        
        (*mono2).capping_prob = CAPPING_PARAM;//Now that it is a barb, mono2 can be capped
        
        
        //If these monomers don't have arp23, now they're allowed to associate with arp23
        
        //unnecessary as of 8/11/08 since arp can associate with free monomers
        //restored 8/26/08
        //if(!((*mono2).mono_has_arp23))
         (*mono2).associate_arp23_prob = ASSOC_ARP23_PARAM;
       // #if (!(INITIALIZE_DIMERS))
        if(!mono_has_arp23) 
         associate_arp23_prob = ASSOC_ARP23_PARAM;
       // #endif
        
        
       
         (*mono2).denucleate_prob = DENUCLEATE_PARAM;
}//end of polymerize()


int monomer::nucleate(monomer *mono2)
{
   if(!((*mono2).minus_is_nucleated))//Can't polymerize to a monomer that's already polymerized... unnecessary (8/11/08)
      if(unif_rand() < nucleate_prob)//check random # against probability
      {
       polymerize(mono2);
       return 1;
      }

  return 0;
}//end of nucleate


int monomer::denucleate()
{  
double rrr = ranf5();//unif_rand();

   if(rrr < denucleate_prob)
   { /*note that this is okay, because denucleate_prob = 0 unless nucleation has occured!*/
       //Change minus' plus pointer & bool, wait until end to change this' minus ptr 
       (*minus).plus = NULL;
       (*minus).plus_is_nucleated = false;
       minus_is_nucleated = false;
       denucleate_prob = 0.;


       //Reset minus pointer,bool
       minus = NULL;       
         
       return 1;
   }
   return 0;
}//end of denucleate


int monomer::cap()
{
    if(unif_rand() < capping_prob)
    {
         plus_is_capped = true;
         nucleate_prob = 0.;
         capping_prob = 0.;
         uncapping_prob = UNCAPPING_PARAM;
         //should this also kill branching prob?
         return 1;
    }
    return 0;
}//end of cap


int monomer::uncap()
{
    if(unif_rand() < uncapping_prob)
    {
        remove_cap();
        return 1;
    }
    return 0;
}//end of uncap


void monomer::remove_cap()
{
  plus_is_capped = false;
  uncapping_prob = 0.;
  if(minus_is_nucleated || mono_has_branch)
  {
    nucleate_prob = NUCLEATE_PARAM;
    capping_prob = CAPPING_PARAM;
  }
}
   


int monomer::associate_arp23(/*arp23_particle *arp23*/)
{
    if(unif_rand() < associate_arp23_prob)
    {
      mono_has_arp23 = true;
//next lines not necessary since arp23 only sticks to f-actin
/*
      if(plus_is_nucleated || minus_is_nucleated)
      { */
//       if(!mono_has_branch)//this line shouldn't be necessary
        branching_prob = BRANCHING_PARAM; //should I just let branching prob be non-zero whenever there is arp? (No. 7/22/08)
//      }
      associate_arp23_prob = 0.;
      dissociate_arp23_prob = DISSOC_ARP23_PARAM;
      //arp23 also protects point from denucleation. may include in future.
      return 1;
    }
    return 0;
}//end of associate_arp23


int monomer::dissociate_arp23(/*arp23_particle *arp23*/)
{
    if(unif_rand() < dissociate_arp23_prob)
    {
      remove_arp23();
      return 1;
    }
    return 0;
}//end of dissociate_arp23

void monomer::remove_arp23()
{
   mono_has_arp23 = false;
   branching_prob = 0.;
   if(plus_is_nucleated || minus_is_nucleated)
     associate_arp23_prob = ASSOC_ARP23_PARAM;
   else
     associate_arp23_prob = 0.;
   dissociate_arp23_prob = 0.;
}

void monomer::construct_branch(monomer *mono2)
{
  //Set bools and pointers
  new_branch = mono2;
  mono_has_branch = true;
  (*mono2).minus = this;
  (*mono2).minus_is_nucleated = true;

  branching_prob = 0.;
  dissociate_arp23_prob = 0.; //Necessary for now. If this line is eliminated, need to add an extra check into branch()
                              //or modify associate_arp23()-- i.e., not let it change branching_param to nonzero.

  (*mono2).debranching_prob = DEBRANCHING_PARAM;
  #if INITIALIZE_DIMERS
  //if(!(*mono2).plus_is_nucleated)//unnecessary...
  (*mono2).nucleate_prob = NUCLEATE_PARAM;
  #endif
  
  //if(!(*mono2).mono_has_branch)//unnecessary since arp only sticks to f-actin
/*
  if((*mono2).mono_has_arp23)
    (*mono2).branching_prob = BRANCHING_PARAM;*/
  
  (*mono2).capping_prob = CAPPING_PARAM;
  (*mono2).associate_arp23_prob = ASSOC_ARP23_PARAM;
}


int monomer::branch(monomer *mono2)
{
  if((*mono2).get_id() != id)
  {
   if(!((*mono2).minus_is_nucleated))
   {
    if(unif_rand() < branching_prob)
    {
       construct_branch(mono2);
       return 1;
    }
   }//when branching, leave denucleate prob. for mono2 at zero (we have debranching for this reason)
  }
  return 0;
}//end of branch


int monomer::debranch()
{   
   if(unif_rand() < debranching_prob)
   {   
       detach_branch();
       return 1;
   }
   return 0;
}//end of debranch     


void monomer::detach_branch()
{
       //ptr and bools for branch base    
       (*minus).new_branch = NULL;
       (*minus).mono_has_branch = false;
       
       //pointer changed to NULL at end
       minus_is_nucleated = false;
       
       
//in the future, may remove arp23 during debranching instead of returning to state with mono and arp23 together
//because there is evidence that this happens
//in that case should I only allow arp23 to associate to nucleated monomers?
//For now, allow arp23 to attach to individual monomers and leave it attached after debranching.
       
       debranching_prob = 0.;

 
       //ptr and bool for this monomer
       minus = NULL;
}//end detach_branch


#endif
