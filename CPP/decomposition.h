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
	decomposition(int myrank, int nb_procs, int nb_procs_x, int nx, int ny, bool iscovering);
	~decomposition();
  int *get_index_x();
	int *get_index_y();
	int *get_index_global();
private:
	int *index_x, *index_y, *index_global;
	int myRank, myRank_x, myRank_y,
			N_procs, N_procs_x, N_procs_y,
			N, Nx, Ny,
			myN, myNx, myNy;
	bool isCovering;
// functions
	bool is_admissable();
	void decompose_x();
	void decompose_y();
	void decompose_global();
};

#endif