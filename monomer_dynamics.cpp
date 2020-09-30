/*
This function is responsible for giving monomers a random displacement at each
time step and for executing any other movements due to systematic forces.
In addition, this subroutine carries out rotational Brownian motion.

For ideas to improve this subroutine and the see 6/17/08 in notes.txt.

Algorithm: loop over all cells(loop over all monomers in cell(check for and execute interactions))
*/

void monomer_dynamics(bool init, unsigned long long step)
{ 
int icounter, jcounter, kcounter;
int ix, jy, kz;
int ix2, jy2, kz2;
int ii, jj, kk;
int dd;
int mm, nn;
int ll;//for bind site
bool avoid, same_cell;


//loop over all cells
for(ii = 0; ii < N[0]; ii++)
{
 for(jj = 0; jj < N[1]; jj++)
 {
   for(kk = 0; kk < N[2]; kk++)
   {
//start the loop over all monomers in cell
nn = get_firstmonoincell(ii,jj,kk);

while(nn != EMPTY)
{
//loop over nearby cells
icounter = 0;
for(ix2 = (ii - 1 + N[0])%N[0]; icounter < 3; ix2 = (ix2+1)%N[0], icounter++)
{
/*//-1 problem fixed. 151028
  if(ix2 != -1)
    ix = ix2;
  else
    ix = N[0]-1;
*/
 ix = ix2;

 jcounter = 0;
 for(jy2 = (jj - 1 + N[1])%N[1]; jcounter < 3; jy2 = (jy2+1)%N[1], jcounter++)
 {

/*
  if(jy2 != -1)
    jy = jy2;
  else
    jy = N[1]-1;
*/
  jy = jy2;

  kcounter = 0;
  for(kz2 = kk; kcounter < 2; kz2 = (kz2+1)%N[2], kcounter++)
  {
      same_cell = false;

      kz = kz2;
      if(kz == kk)
      {
      if(jy2 == (jj-1+N[1])%N[1] || (jy == jj && ix2 == (ii-1+N[0])%N[0]))
      {
         continue;
      }

        if(jy == jj)
        if(ix == ii)
           same_cell = true;

      }

//compare all monomers mm in particular nearby cell to monomer nn
    mm = get_firstmonoincell(ix,jy,kz);

    while(mm != LAST)
    {
        avoid = false;
        if(same_cell)
        if(mm <= nn)
          avoid = true;

	if(MODIFIED_SHELL_EXC_VOL)
	{//in this case, don't execute excl vol interactions between shell monos
	   if((mm >= NUMBER_IN_POLYMER) && (nn >= NUMBER_IN_POLYMER))
		avoid = true;
	}

//apply systematic forces and torques first, if applicable
      if(!avoid)
      {     
//if(mono_list[nn].calculate_distance_sq(&mono_list[mm],PREVIOUS) < 1.)
//fprintf(stderr, "%i %i m1 %g %g %g m2 %g %g %g\n", nn, mm, mono_list[nn].get_pos(0),mono_list[nn].get_pos(1), mono_list[nn].get_pos(2), mono_list[mm].get_pos(0), mono_list[mm].get_pos(1), mono_list[mm].get_pos(2));
       (mono_list[nn]).sys_force(&mono_list[mm]);
//if(mono_list[nn].calculate_distance_sq(&mono_list[mm],PREVIOUS) < 1.)
//fprintf(stderr, "%i %i m1 %g %g %g m2 %g %g %g\n", nn, mm, mono_list[nn].get_pos(0),mono_list[nn].get_pos(1), mono_list[nn].get_pos(2), mono_list[mm].get_pos(0), mono_list[mm].get_pos(1), mono_list[mm].get_pos(2));


//no binding//
//	if(!init)
//150614: take out binding entirely for now.
//170510: binding between polymer monos added back in (previous was between polymer and shell) and only for "binding monos"
// 	if((mono_list[nn].get_binding_mono())||(mono_list[mm].get_binding_mono()))
//	if(((nn >= NUMBER_IN_POLYMER)&&(mm<NUMBER_IN_POLYMER))||((mm >= NUMBER_IN_POLYMER)&&(nn<NUMBER_IN_POLYMER)))
	if((nn < NUMBER_IN_POLYMER)&&(mm<NUMBER_IN_POLYMER))
		binding_func(nn, mm);
 
       } //avoid
       mm = get_monolinklist(mm);
    }//bracket closes while(mm !=-1)
  }//kz loop
 }//jy loop
}//closes for loop over ix
// Brownian motion


	if(!init)
	{
	 if(THERMAL)
           (mono_list[nn]).translational_brownian_move(step);

	 if(COMPRESS)
	 {
		compress(nn);
         }

	 if(nn >= NUMBER_IN_POLYMER)
	 {
//excess feature commented out         
//        (mono_list[nn]).outward_pressure();

		if(!COMPRESS)
		if(mono_list[nn].get_load_mono())
		{
		  (mono_list[nn]).load();
		}
	 }// if(nn >= NUMBER_IN_POLYMER)
	}//!init

nn = get_monolinklist(nn);

}//bracket closes while(nn != -1)
}//bracket closes for(kk = 0; kk < N[2]; kk++)
}//bracket closes for(jj = 0; jj < N[1]; jj++)
}//bracket closes for(ii = 0; ii < N[0]; ii++)
 

}


