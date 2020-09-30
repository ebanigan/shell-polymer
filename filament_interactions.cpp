/*Cycles through all mononmer pairwise connections and modfies position based on
bond and bending (and torsion, if applicable) potentials*/

void filament_interactions()
{
     long ii,jj;

       for(ii = 0; ii < pairs.size(); ii++)
       {
        (pairs[ii]).bond_potential();
	bending_func(ii);
       }//for loop over pairs



	for(ii = 0; ii < crosslinkpairs.size(); ii++)
	{
                //if(ii==0)
                     //fprintf(stderr, "cl %g\n", crosslinkpairs[ii].get_bond_strength());
		crosslinkpairs[ii].bond_potential();
	}//for loop over crosslink pairs
}


void bending(int ii)
{
        int jj = (*(pairs[ii].first)).get_id();
        if(jj < NUMBER_IN_POLYMER)//shell has no bend pot., so only compute for polymer.
          jj = (*(pairs[ii].first)).get_pair_num2();
        if(jj > -1)
        {
          (pairs[ii]).bending_potential(&pairs[jj]);
        }
}


void bending_null(int ii){}
