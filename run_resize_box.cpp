void run_resize_box()
{
                double min_shell = 9.*L[0];
                double max_shell = -9.*L[0];
                double max_y = -9.*L[1];
                double max_z = -9.*L[2];
                double min_y = 9.*L[1];
                double min_z = 9.*L[2];
                double monox, monoy, monoz;

            for(int ii2 = NUMBER_IN_POLYMER; ii2 < NUMBER_OF_MONOMERS; ii2++)
            {
                monox = mono_list[ii2].get_prev_pos(0);
                monoy = mono_list[ii2].get_prev_pos(1);
                monoz = mono_list[ii2].get_prev_pos(2);

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

            resize_box(max_shell-min_shell, max_y-min_y, max_z-min_z);
}



void resize_box(double shell_width, double transverse1, double transverse2)
{
  double max_predicted_size = 1.25*(shell_width+0.75*SHELL_RADIUS);
  bool box_changed = false;
  
  if(max_predicted_size > LX)
  {
   fprintf(stderr, "extending LX\n");
   cl_lx = 1.14*max_predicted_size;
   box_changed = true;
  }
  else if(max_predicted_size < 0.73*LX)
  {
   fprintf(stderr, "shrinking LX\n");
   cl_lx = 0.9*LX;
   box_changed = true;
  }
  if(transverse1 < 0.6*LY)
  {
   if(cl_ly > 6.*MONO_DIAM)
   if(transverse2 < 0.6*LZ)
   {
     cl_ly = 0.9*LY;
     box_changed = true;
     fprintf(stderr, "shrinking LY, LZ\n");
   }
  }
  else if(1.38*transverse1 > LY)
  {
    cl_ly = 1.09*LY;
    box_changed = true;
    fprintf(stderr, "extending LY, LZ\n");
  }
  else if(1.38*transverse2 > LY)
  {
    cl_ly = 1.09*LY;
    box_changed = true;
    fprintf(stderr, "extending LY, LZ\n");
  }

  if(box_changed)
  {
   initialize_box();
   recenter();
   initialize_neighbor_list();
   update_mono_list();
  }

}

