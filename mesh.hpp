#ifndef MESH_HPP
#define MESH_HPP

class mesh
{
public:

	mesh(const double& S, const double& T, const int& n_x, const int& n_t, const double& sigma);
	
	void print_mesh();
	
	double get_dt();
	double get_dx();
	double get_Smax();
	double get_Smin();
	double get_S();
	double get_mat();
	double get_sigma();
	
	int get_nx();
	int get_nt();
	
private:
	
	int m_nx;
	int m_nt;
	
	double m_Smax;
	double m_Smin;
	double m_dt;
	double m_dx;
	double m_S;
	double m_sigma;
};
	

#endif
