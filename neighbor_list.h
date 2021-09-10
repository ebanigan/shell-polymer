#include "monomer_pair.h"
//#include "dynamic_monomer_functions.h"


 extern vector<int> monolinklist;
 //with vector<long> monolinklist not stored in the function, it no longer matters if I initialize right away. 
 extern vector<vector<vector<int> > > firstmonoincell;




/************/ 
//Functions for neighborlist
/***********/
        void initialize_neighbor_list()
        {
          int ii, jj, kk;

	monolinklist.clear();
	firstmonoincell.clear();

        for(ii = 0; ii < NUMBER_OF_MONOMERS; ii++)
         monolinklist.push_back(LAST);

        vector<int> temp_vector(N[2],EMPTY);
	vector<vector<int> > temp_2dvector;
        for(ii = 0; ii < N[1]; ii++)
	  temp_2dvector.push_back(temp_vector);
	for(jj = 0; jj < N[0]; jj++)
	  firstmonoincell.push_back(temp_2dvector);
        }//initialize_neighbor_list

        int get_firstmonoincell(int xx, int yy, int zz){return firstmonoincell[xx][yy][zz];}
        int get_monolinklist(int ii){return monolinklist[ii];}

	void update_mono_list()
	{
          int nn;
          int ii, jj, kk;
          int dd;
          int temp;
          int cellnum[DIMENSION];

          for(ii = 0; ii < N[0]; ii++)
           for(jj = 0; jj < N[1]; jj++)
            for(kk = 0; kk < N[2]; kk++)
             firstmonoincell[ii][jj][kk] = EMPTY;

          for(nn = 0; nn < NUMBER_OF_MONOMERS; nn++)
          {
             for(dd = 0; dd < DIMENSION; dd++)
             {
               cellnum[dd] = (int)((mono_list[nn]).get_pos(dd) * inv_cell_length[dd]);
	     }

             monolinklist[nn] = firstmonoincell[cellnum[0]][cellnum[1]][cellnum[2]];
             firstmonoincell[cellnum[0]][cellnum[1]][cellnum[2]] = nn;
	  }//nn

	}//end of update

