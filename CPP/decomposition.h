#ifndef DECOMPOSITION_H
#define DECOMPOSITION_H

/** decomposes a given rectangular mesh into rectangular submeshes
 */
class decomposition{
	decomposition();
	decomposition(int myrank, int nb_procs, int nb_procs_x, int nx, int ny);
	~decomposition();
public:
	*int get_index_x();
	*int get_index_y();
	*int get_index_global();
private:
	int *index_x, *index_y, *index_global;
	int myRank, N_procs, N_procs_x, N_procs_y,
			N, Nx, Ny;
// functions
	bool is_admissable();
	void decompose_x();
	void decompose_y();
	void decompose_global();
};

#endif