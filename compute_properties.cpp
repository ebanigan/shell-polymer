

double compute_applied_force()
{
double tot_spring_extension=0.;

for(int ii = NUMBER_IN_POLYMER; ii < NUMBER_OF_MONOMERS; ii++)
{

         //compute force
        if(mono_list[ii].get_load_mono())
        {  
           tot_spring_extension += mono_list[ii].get_sign_load()*(0.5*L[0] + mono_list[ii].get_sign_load()*(pipette_dist_from_center - mono_list[ii].get_load_rest_length()) - mono_list[ii].get_prev_pos(0));
        }
}

return PIPETTE_STIFFNESS*tot_spring_extension;

}




vector<double> compute_extension()
{
double min_shell = 9.*L[0];
double max_shell = -9.*L[0];
double max_y = -9.*L[1];
double max_z = -9.*L[2];
double min_y = 9.*L[1];
double min_z = 9.*L[2];
double monox, monoy, monoz;
vector<double> exts;

for(int ii = NUMBER_IN_POLYMER; ii < NUMBER_OF_MONOMERS; ii++)
{
        monox = mono_list[ii].get_prev_pos(0);
        monoy = mono_list[ii].get_prev_pos(1);
        monoz = mono_list[ii].get_prev_pos(2);


        if(monox < min_shell)
                min_shell = monox;
        else if(monox > max_shell)
                max_shell = monox;

        if(monoy < min_y)
                min_y = monoy;
        else if(monoy > max_y)
                max_y = monoy;

        if(monoz < min_z)
                min_z = monoz;
        else if(monoz > max_z)
                max_z = monoz;
}

exts.push_back(max_shell-min_shell);
exts.push_back(max_y-min_y);
exts.push_back(max_z-min_z);

return exts; 
}
