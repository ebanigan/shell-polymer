

bool release_pipette(unsigned long long tstep)
{

if(!pipette_released)// if it is released, don't do checks again
{

 if(THRESHOLD_EXT > 0.)
 {
   vector<double> exts=compute_extension();
   if(exts[0] > THRESHOLD_EXT)
   {
       pipette_released=true;
   }
 }
 else if(THRESHOLD_FORCE > 0.)
 {
   if(compute_applied_force() > THRESHOLD_FORCE)
   {
       pipette_released=true;
   }
 }
 else // threshold time
 {
   if(tstep > RELEASE_TIME)
   {
       pipette_released=true;
   }
 }
}

return pipette_released;

}



bool end_condition(unsigned long long tstep)
{
bool end_sim = false;

//For now let's just have zero force or a release relaxation time as the possible end conditions
//because 0 extension we would need good measurement at beginning of sim
if(END_ZERO_FORCE)
{
  if(compute_applied_force() < 0.)
    end_sim = true;
}
else if(tstep > RELEASE_RELAX_TIME)
{
    end_sim = true;
}

return end_sim;
}
