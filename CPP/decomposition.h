#ifndef DECOMPOSITION_H
#define DECOMPOSITION_H

/** decomposes a given rectangular mesh into rectangular submeshes
 *  where nb_Procs is the total number of processors that are distributed 
 *  as follows: nb_Procs_x in x direction, nb_Procs_y in y direction
 *  the user needs to make sure that this distribution is possible
 */
class decomposition{
public:
	decomposition();
	decomposition(int myrank, int nb_procs, int nb_procs_x, int nx, int ny);
	~decomposition();
  int *get_index_x();
	int *get_index_y();
	int *get_index_global();
	int *get_index_global_top();
	int *get_index_global_bottom();
	int *get_index_global_left();
	int *get_index_global_right();
	int get_myNx();
	int get_myNy();
	int get_myN();
private:
	int *index_x=NULL, *index_y=NULL, *index_global=NULL,
			*index_global_top=NULL, *index_global_bottom=NULL, 
			*index_global_left=NULL, *index_global_right=NULL;
	int myRank, myRank_x, myRank_y,
			N_procs, N_procs_x, N_procs_y,
			N, Nx, Ny,
			myN, myNx, myNy;
// functions
	bool is_admissable();
	void decompose();
	void decompose_x();
	void decompose_y();
	void decompose_global();
	void accumulate_global_borders();
	void accumulate_global_top();
	void accumulate_global_bottom();
	void accumulate_global_left();
	void accumulate_global_right();
};

#endif