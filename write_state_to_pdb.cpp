#define HIDE_G_ACTIN false
#define PDB_SKIP (4*NUMSKIP)//NUMSKIP//(4*NUMSKIP)//changed at 12:13pm 6/19/15 //(2*NUMSKIP)//10
//for movies:
//#define PDB_SKIP (2*NUMSKIP)
#define HIGHLIGHT_CENTER true

void write_state_to_pdb(unsigned long long step)
{
//too much data. skip every 5 pdb files
if(step%(5*PDB_SKIP) == 0)
{
  int dummy;
  double scale_factor = L[0]/100.;//(L[DIMENSION-1] / 100.);
  FILE *outfile;
  char outname[72];
  char dirname[72];
  
  sprintf(dirname, "output/pdb_files%6.6i", TRIALNUMBER);
//  check_directory(dirname);
  
  
  if(step == 0)
  {
    char command[64];
    sprintf(command, "mkdir output/pdb_files%6.6i", TRIALNUMBER);
    dummy= system(command);
  }
  
  
  sprintf(outname, "output/pdb_files%6.6i/pdb_file%6.6llu.pdb", TRIALNUMBER, step/PDB_SKIP);
  outfile = fopen(outname, "w");
  
  char color;
  int m;
  for(m = 0; m < NUMBER_OF_MONOMERS; m++)
  {
    if(m < NUMBER_IN_POLYMER)
    {
     if(!(mono_list[m].get_crosslinked() || mono_list[m].get_shell_poly_link()))
     {
        color = 'E';
     }
     else if(mono_list[m].get_crosslinked() && mono_list[m].get_shell_poly_link())
     {
       color = 'C';
     }
     else if(mono_list[m].get_crosslinked())
     {
        color = 'B';
     }
     else
     {
	color = 'F';
     }
    }
    else
    {
      if(!mono_list[m].get_shell_poly_link())
      {
	if(mono_list[m].get_load_mono())
          color = 'S';
	else
          color = 'T';
      }
      else
	color = 'G';
    }
//use these statements to not print part of shell.  
//if(!((color == 'T') && (mono_list[m].get_pos(0)-0.5*LX > 0.) && (mono_list[m].get_pos(1)-0.5*LY > 0.) && (mono_list[m].get_pos(2)-0.5*LZ > 0.)) )
//if(!((color == 'T') && (mono_list[m].get_pos(0)-0.5*LX > 0.) && (mono_list[m].get_pos(1)-0.5*LY < 0.) && (mono_list[m].get_pos(2)-0.5*LZ > 0.)) )
      fprintf(outfile, "ATOM  %5i  O   ASP %c%4i      %6.3f  %6.3f  %6.3f  1.00  0.00\n", m+1, color, 1, mono_list[m].get_pos(0)/scale_factor, mono_list[m].get_pos(1)/scale_factor, mono_list[m].get_pos(2)/scale_factor);
  }

if(SOLID_INTERIOR)
{
  color = 'C';
  fprintf(outfile, "ATOM  %5i  O   ASP %c%4i      %6.3f  %6.3f  %6.3f  1.00  0.00\n", m+1, color, 1, central_mono.get_pos(0)/scale_factor, central_mono.get_pos(1)/scale_factor, central_mono.get_pos(2)/scale_factor);
}
  
  fflush(outfile);
  fclose(outfile);
}//skip every 5 pdb files...

}

