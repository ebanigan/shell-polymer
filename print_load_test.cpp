void print_load_test(long step)
{

int ii;
int mono1, mono2, tempm;
bool crossing_left_cap, crossing_left_side, crossing_middle, crossing_right_side, crossing_right_cap;
double x1, x2, tempx;

FILE *loadtestfile;
char loadtestname[96];
sprintf(loadtestname, "output/loadtest%6.6i", TRIALNUMBER);
if(step == 0)
	loadtestfile = fopen(loadtestname, "w");
else
	loadtestfile = fopen(loadtestname, "a");


fprintf(loadtestfile, "%li %g %g %g %g %g %g %g\n", step, mono_list[spring_mono_id].get_prev_pos(0) - spring_att_point[0], mono_list[310].get_prev_pos(0) - spring_att_point[0], mono_list[221].get_prev_pos(0) - spring_att_point[0], mono_list[287].get_prev_pos(0) - spring_att_point[0], mono_list[311].get_prev_pos(0) - spring_att_point[0], mono_list[363].get_prev_pos(0) - spring_att_point[0], mono_list[97].get_prev_pos(0) - spring_att_point[0]);

fflush(loadtestfile);
fclose(loadtestfile);

//commented out below: code identifying pairs that traverse 5 cuts of the shell
/*
for(ii = 0; ii < pairs.size(); ii++)
{
	crossing_left_cap = false;
	crossing_left_side = false;
	crossing_middle = false;
	crossing_right_side = false;
	crossing_right_cap = false;
	mono1 = (*(pairs[ii].first)).get_id();
	mono2 = (*(pairs[ii].second)).get_id();
	x1 = mono_list[mono1].get_prev_pos(0);
	x2 = mono_list[mono2].get_prev_pos(0);
	if(x1 > x2)
	{
	  tempx = x1;
	  x1 = x2;
	  x2 = tempx;
	  tempm = mono1;
	  mono1 = mono2;
	  mono2 = tempm;
	}

	if(x1 < shell_cm[0] - SHELL_RADIUS*(1.0-cl_load_frac))
	{
		if(x2 > shell_cm[0] - SHELL_RADIUS*(1.0-cl_load_frac))
		{
		  crossing_left_cap = true;
		  if(x2 > shell_cm[0] - 0.5*SHELL_RADIUS*(1.0-2.0*cl_load_frac))
		  {
		    crossing_left_side = true;
                    if(x2 > shell_cm[0])
                    {
		      crossing_middle = true;
		      if(x2 > shell_cm[0] + 0.5*SHELL_RADIUS*(1.0-2.0*cl_load_frac))
		      {
		        crossing_right_side = true;
			if(x2 > shell_cm[0] + SHELL_RADIUS*(1.0-cl_load_frac))
		        {
		          crossing_right_cap = true;
		        }//crossing_right_cap
		      }//crossing_right_side
		    }//crossing_middle
		  }//crossing_left_side			
		}//crossing_left_cap
	}//x1 < left cap line
	else if(x1 < shell_cm[0] - 0.5*SHELL_RADIUS*(1.0-2.0*cl_load_frac))
	{
                if(x2 > shell_cm[0] - 0.5*SHELL_RADIUS*(1.0-2.0*cl_load_frac))
                {
                  crossing_left_side = true;
                  if(x2 > shell_cm[0])
                  {
                    crossing_middle = true;
                    if(x2 > shell_cm[0] + 0.5*SHELL_RADIUS*(1.0-2.0*cl_load_frac))
                    {
                      crossing_right_side = true;
                      if(x2 > shell_cm[0] + SHELL_RADIUS*(1.0-cl_load_frac))
                      {
                        crossing_right_cap = true;
                      }//crossing_right_cap
                    }//crossing_right_side
                  }//crossing_middle
                }//crossing_left_side 
	}//x1 < left side line
	else if(x1 < shell_cm[0])
	{
                if(x2 > shell_cm[0])
                {
                  crossing_middle = true;
                  if(x2 > shell_cm[0] + 0.5*SHELL_RADIUS*(1.0-2.0*cl_load_frac))
                  {
                    crossing_right_side = true;
                    if(x2 > shell_cm[0] + SHELL_RADIUS*(1.0-cl_load_frac))
                    {
                      crossing_right_cap = true;
                    }//crossing_right_cap
                  }//crossing_right_side
                }//crossing_middle
	}//x1 < middle line
	else if(x1 < shell_cm[0] + 0.5*SHELL_RADIUS*(1.0-2.0*cl_load_frac))
	{
                if(x2 > shell_cm[0] + 0.5*SHELL_RADIUS*(1.0-2.0*cl_load_frac))
                {
                  crossing_right_side = true;
                  if(x2 > shell_cm[0] + SHELL_RADIUS*(1.0-cl_load_frac))
                  {
                    crossing_right_cap = true;
                  }//crossing_right_cap
                }//crossing_right_side
	}//x1 < right side line
	else if(x1 < shell_cm[0] + SHELL_RADIUS*(1.0-cl_load_frac))
	{
                if(x2 > shell_cm[0] + SHELL_RADIUS*(1.0-cl_load_frac))
                {
                  crossing_right_cap = true;
                }//crossing_right_cap
	}//x1 < right cap line
}//for(ii...
*/


}

