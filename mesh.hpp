#ifndef MESH_HPP
#define MESH_HPP

class mesh
{
public:

	mesh(const double& S, const double& T, const int& n_x, const int& n_t, const double& sigma);
	
	void print_mesh();
	
private:
	int m_nx;
	int m_nt;
	
	double m_Smax;
	double m_Smin;
	double m_dt;
	double m_dx;
};
	

#endif