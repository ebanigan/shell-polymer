void compress(int ii)
{

double dist_bottom = BOTTOM_WALLPOS + WALLTHICKNESS + mono_list[ii].get_radius() - mono_list[ii].get_prev_pos(2);
double dist_top = TOP_WALLPOS - WALLTHICKNESS - mono_list[ii].get_radius() - mono_list[ii].get_prev_pos(2); 
double compressive_force;

if(dist_bottom > 0.)
{
	compressive_force = WALLSPRING * dist_bottom * dt;
	mono_list[ii].move_1d(2, compressive_force*mono_list[ii].get_tdiffusion_coeff());
	BOTTOM_WALLFORCE -= compressive_force * WALLDIFF;
}
else if(dist_top < 0.)
{
	compressive_force = WALLSPRING * dist_top * dt;
	mono_list[ii].move_1d(2, compressive_force*mono_list[ii].get_tdiffusion_coeff());
	TOP_WALLFORCE -= compressive_force * WALLDIFF;
}


}

