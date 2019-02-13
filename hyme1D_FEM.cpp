#include <cstddef>
#include <cstdlib>

#include <fstream>
#include <iostream>

#include "hyme1D_FEM.h"

/*
1. Linear shape function for both displacement and pore pressure
2. Complete gauss integration for both
*/

namespace hyme1D_FEM
{

// Model Data
struct GaussPoint
{
	double n;
	double density_avg;

	double s11;
	double e11;

	double xi;
	double dx_dxi;
	double dxi_dx;

	double N1, N2;
	double dN1_dx, dN2_dx;
};

struct Element
{
	// topology
	size_t index_x;

	double density_s, density_f;

	double E;
	double k; // permeability
	double miu; // dynamic viscosity
	double k_div_miu; // k divided by miu

	// All physical parameters needed to be updated
	// are contained in gauss points.
	GaussPoint gp1, gp2;
};

struct Node
{
	double x;

	double pa, pv, p;
	double a, v, u, w;
	double du; // for stress increment
};

struct BodyForce
{
	double bf_x;
	size_t elem_id;
};

struct SurfaceForce
{
	double sf_x;
	size_t elem_id;
	double x;
};

// Acceleration boundary conditions (Dirchlet boundary conditions)
struct AccelerationBC
{
	double a;
	size_t node_id;
};

// init according to pore pressure BC
// Assume that dp_dt = 0 where p is specified.
struct PAccelerationBC
{
	double pa;
	size_t node_id;
};

// Velocity initial boundary conditions
struct VelocityBC
{
	double v;
	size_t node_id;
};

// Displacement initial boundary conditions
struct DisplacementBC
{
	double d;
	size_t node_id;
};

// pore pressure initial boundary conditions
struct PorePressureBC
{
	double p;
	size_t node_id;
};

size_t node_num;
Node *nodes;

size_t elem_num;
Element *elems;

// Dirichlet boundary conditions
// We assume that the same body force is applied to both solid and fluid phase.
size_t body_force_num;
BodyForce *body_force;

size_t surface_force_num;
SurfaceForce *surface_force;

// Neumann boundary conditions
size_t acceleration_num;
AccelerationBC *acceleration_bc;

// assume that dp_dt = 0 where p is specified as bc.
size_t p_acceleration_num;
PAccelerationBC *p_acceleration_bc;

// Initial boundary conditions
size_t displacement_num;
DisplacementBC *displacement_bc;

size_t velocity_num;
VelocityBC *velocity_bc;

size_t pore_pressure_num;
PorePressureBC *pore_pressure_bc;

// Step time length
double total_t, cur_t;
// Substep time increment
double dt;

// Total degrees of freedom
size_t total_dof;
// Temporary calculation variables
// v + 0.5 * a * dt
double *tmp_coef3;
// v * dt + 0.25 * a * dt * dt
double *tmp_coef4;

inline void calNaturalCoord(Element *elem, double *coord, double *na_coord)
{
	size_t id_x = elem->index_x;
	double n1_x = nodes[id_x].x;
	double n2_x = nodes[id_x + 1].x;
	na_coord[0] = (2.0 * coord[0] - (n2_x + n1_x)) / (n2_x - n1_x);
}

inline double dx_dxi(Element *elem, double xi)
{
	size_t id_x = elem->index_x;
	return (nodes[id_x + 1].x - nodes[id_x].x) / 2.0;
}

#define _N1(xi) ((1.0 - (xi))/2.0)
#define _N2(xi) ((1.0 + (xi))/2.0)
#define _dN1_dxi(xi) (-0.5)
#define _dN2_dxi(xi) (0.5)

// initialize calculation
int init(void)
{
	cur_t = 0.0;

	// Generate Dirichlet boundary conditions for pore pressure
	p_acceleration_num = pore_pressure_num;
	if (p_acceleration_num)
		p_acceleration_bc = (PAccelerationBC *)malloc(sizeof(PAccelerationBC) * p_acceleration_num);
	else
		p_acceleration_bc = NULL;
	for (size_t i = 0; i < p_acceleration_num; i++)
	{
		// Assume that dp_dt = 0 where p is specified.
		p_acceleration_bc[i].pa = 0.0;
		p_acceleration_bc[i].node_id = pore_pressure_bc[i].node_id;
	}

	// Apply initial boundary conditions
	// Acceleration
	for (size_t i = 0; i < acceleration_num; i++)
	{
		nodes[acceleration_bc[i].node_id].a = acceleration_bc[i].a;
	}
	// Velocity 
	for (size_t i = 0; i < velocity_num; i++)
	{
		nodes[velocity_bc[i].node_id].v = velocity_bc[i].v;
	}
	// Displacement 
	for (size_t i = 0; i < displacement_num; i++)
	{
		nodes[displacement_bc[i].node_id].u = displacement_bc[i].d;
		// update stress
	}
	// Pore pressure
	for (size_t i = 0; i < pore_pressure_num; i++)
	{
		nodes[pore_pressure_bc[i].node_id].p = pore_pressure_bc[i].p;
	}

	return 0;
};

// Calculate each substeps
int substep(void);
int output(std::ofstream &out, double cur_time, Node *nodes, size_t node_num);

// calculate the whole step
// simplify ver.
int cal(std::ofstream &out)
{
	while (cur_t < total_t)
	{
		substep();
		cur_t += dt;
		output(out, cur_t, nodes, node_num);
	}

	return 0;
}


int hyme1D_FEM(void)
{
	/* init mesh data */
	node_num = 21;
	elem_num = node_num - 1;
	nodes = (Node *)malloc(sizeof(Node) * node_num);
	elems = (Element *)malloc(sizeof(Element) * elem_num);

	for (size_t i = 0; i < node_num; i++)
	{
		nodes[i].x = (double)i;

		// default initialization
		nodes[i].pa = 0.0;
		nodes[i].pv = 0.0;
		nodes[i].p = 0.0;
		nodes[i].a = 0.0;
		nodes[i].v = 0.0;
		nodes[i].u = 0.0;
		nodes[i].w = 0.0;
	}

	for (size_t i = 0; i < elem_num; i++)
	{
		elems[i].index_x = i;

		elems[i].density_s = 3000.0;
		elems[i].density_f = 1000.0;

		elems[i].E = 40.0e6;
		elems[i].k = 1.0e-3;
		elems[i].miu = 1.0;
		
		// init gauss point
		elems[i].gp1.s11 = 0.0;
		elems[i].gp1.e11 = 0.0;
		elems[i].gp1.n = 0.5;
		elems[i].gp2.s11 = 0.0;
		elems[i].gp2.e11 = 0.0;
		elems[i].gp2.n = 0.5;

		// default initialization
		elems[i].k_div_miu = elems[i].k / elems[i].miu;
	}

	/* init BCs */
	// Body force
	body_force_num = elem_num;
	body_force = (BodyForce *)malloc(sizeof(BodyForce) * body_force_num);
	for (size_t i = 0; i < body_force_num; i++)
	{
		body_force[i].bf_x = 0.0;
		body_force[i].elem_id = i;
	}
	// Surface force
	surface_force_num = 1;
	surface_force = (SurfaceForce *)malloc(sizeof(SurfaceForce) * surface_force_num);
	surface_force[0].sf_x = -40.0e3;
	surface_force[0].elem_id = elem_num - 1;
	surface_force[0].x = (double)(node_num - 1);

	// Acceleration BC
	acceleration_num = 1;
	acceleration_bc = (AccelerationBC *)malloc(sizeof(AccelerationBC) * acceleration_num);
	acceleration_bc[0].a = 0.0;
	acceleration_bc[0].node_id = 0;

	// Pore pressure BC
	pore_pressure_num = 1;
	pore_pressure_bc = (PorePressureBC *)malloc(sizeof(PorePressureBC) * pore_pressure_num);
	pore_pressure_bc[0].p = 0.0;
	pore_pressure_bc[0].node_id = node_num - 1;
	// Velocity BC 
	velocity_num = 0;
	velocity_bc = NULL;
	// Displacement BC
	displacement_num = 0;
	displacement_bc = NULL;

	std::ofstream out("res.txt");

	total_dof = node_num + node_num;
	// v + 0.5 * a * dt
	tmp_coef3 = (double *)malloc(sizeof(double) * node_num);
	// v * dt + 0.25 * a * dt * dt
	tmp_coef4 = (double *)malloc(sizeof(double) * total_dof);

	// Initialize calculation
	init();

	total_t = 1.0;
	dt = 0.01;
	// Calculate
	cal(out);

	// Delete data
	if (nodes) free(nodes);
	if (elems) free(elems);
	// Dirichlet boundary conditions
	if (acceleration_bc) free(acceleration_bc);
	if (p_acceleration_bc) free(p_acceleration_bc);
	// Neumann boundary conditions
	if (body_force) free(body_force);
	if (surface_force) free(surface_force);
	// Initial boundary conditions
	if (velocity_bc) free(velocity_bc);
	if (displacement_bc) free(displacement_bc);
	if (pore_pressure_bc) free(pore_pressure_bc);
	
	// Output file
	out.close();
	
	if (tmp_coef3) free(tmp_coef3);
	if (tmp_coef4) free(tmp_coef4);

	return 0;
}


int substep(void)
{
	// Global stiffness matrix
	Eigen::SparseMatrix<double> g_st_mat(total_dof, total_dof);
	std::list<Eigen::Triplet<double> > g_st_mat_coef;
	// Global force vector
	Eigen::VectorXd g_f_vec(total_dof);
	// Global acceleration vector
	Eigen::VectorXd g_a_vec(total_dof);
	
	// Elemental stiffness matrix
	Eigen::Matrix<double, 4, 4> elem_st_mat;
	// Elemental force vector
	Eigen::VectorXd elem_f_vec(4);
	Element *cur_elem;
	GaussPoint *cur_gp1, *cur_gp2;
	Node *cur_node1, *cur_node2;
	/* Index mapping array
		used to map coefficient of elemental stiffness
		matrix into glabal stiffness matrix.
	*/
	size_t cur_g_id[4];

	// Temporary calculation param
	const static double dt_div_2 = 0.5 * dt;
	const static double dt_square_div_4 = 0.25 * dt * dt;
	for (size_t i = 0; i < node_num; i++)
	{
		// v + 0.5 * a * dt
		tmp_coef3[i] = nodes[i].v + dt_div_2 * nodes[i].a;
		// v * dt + 0.25 * a * dt * dt
		tmp_coef4[i] = dt * nodes[i].v + dt_square_div_4 * nodes[i].a;
		tmp_coef4[node_num + i] = dt * nodes[i].pv + dt_square_div_4 * nodes[i].pa;
	}

	g_f_vec.setZero();

	// Calculate elemental stiffness matrix and force vector
	for (size_t elem_id = 0; elem_id < elem_num; elem_id++)
	{
		cur_elem = elems + elem_id;
		cur_node1 = nodes + cur_elem->index_x;
		cur_node2 = nodes + cur_elem->index_x + 1;
		cur_gp1 = &(cur_elem->gp1);
		cur_gp2 = &(cur_elem->gp2);

		// init index mapping array
		cur_g_id[0] = cur_elem->index_x;
		cur_g_id[1] = cur_elem->index_x + 1;
		cur_g_id[2] = node_num + cur_elem->index_x;
		cur_g_id[3] = node_num + cur_elem->index_x + 1;

		// init integration point
#define INIT_INTEGRATION_POINT(elem, gp, gp_coord)           \
		{                                                    \
			(gp).xi = gp_coord;                              \
			(gp).N1 = _N1((gp).xi);                          \
			(gp).N2 = _N2((gp).xi);                          \
			(gp).dx_dxi = dx_dxi(&(elem), (gp).xi);          \
			(gp).dxi_dx = 1.0 / (gp).dx_dxi;                 \
			(gp).dN1_dx = _dN1_dxi(gp_coord) * (gp).dxi_dx;  \
			(gp).dN2_dx = _dN2_dxi(gp_coord) * (gp).dxi_dx;  \
		} while (false)
		INIT_INTEGRATION_POINT(*cur_elem, *cur_gp1, -0.57735);
		INIT_INTEGRATION_POINT(*cur_elem, *cur_gp2,  0.57735);

		// Init elemental stiffness matrix
		elem_st_mat.setZero();
		// M
		cur_gp1->density_avg = (1.0 - cur_gp1->n) * cur_elem->density_s + cur_gp1->n * cur_elem->density_f;
		cur_gp2->density_avg = (1.0 - cur_gp2->n) * cur_elem->density_s + cur_gp2->n * cur_elem->density_f;
		elem_st_mat(0, 0) = cur_gp1->density_avg * cur_gp1->N1 * cur_gp1->N1 * cur_gp1->dx_dxi + cur_gp2->density_avg * cur_gp2->N1 * cur_gp2->N1 * cur_gp2->dx_dxi;
		elem_st_mat(0, 1) = cur_gp1->density_avg * cur_gp1->N1 * cur_gp1->N2 * cur_gp1->dx_dxi + cur_gp2->density_avg * cur_gp2->N1 * cur_gp2->N2 * cur_gp2->dx_dxi;
		elem_st_mat(1, 0) = cur_gp1->density_avg * cur_gp1->N2 * cur_gp1->N1 * cur_gp1->dx_dxi + cur_gp2->density_avg * cur_gp2->N2 * cur_gp2->N1 * cur_gp2->dx_dxi;
		elem_st_mat(1, 1) = cur_gp1->density_avg * cur_gp1->N2 * cur_gp1->N2 * cur_gp1->dx_dxi + cur_gp2->density_avg * cur_gp2->N2 * cur_gp2->N2 * cur_gp2->dx_dxi;
		elem_st_mat(2, 0) = cur_elem->density_f * cur_elem->k_div_miu * (cur_gp1->dN1_dx * cur_gp1->N1 * cur_gp1->dx_dxi + cur_gp2->dN1_dx * cur_gp2->N1 * cur_gp2->dx_dxi);
		elem_st_mat(2, 1) = cur_elem->density_f * cur_elem->k_div_miu * (cur_gp1->dN1_dx * cur_gp1->N2 * cur_gp1->dx_dxi + cur_gp2->dN1_dx * cur_gp2->N2 * cur_gp2->dx_dxi);
		elem_st_mat(3, 0) = cur_elem->density_f * cur_elem->k_div_miu * (cur_gp1->dN2_dx * cur_gp1->N1 * cur_gp1->dx_dxi + cur_gp2->dN2_dx * cur_gp2->N1 * cur_gp2->dx_dxi);
		elem_st_mat(3, 1) = cur_elem->density_f * cur_elem->k_div_miu * (cur_gp1->dN2_dx * cur_gp1->N2 * cur_gp1->dx_dxi + cur_gp2->dN2_dx * cur_gp2->N2 * cur_gp2->dx_dxi);
		//std::cout << "elemental stiffness matrix: " << std::endl << elem_st_mat << std::endl;
		// + 0.5 * dt * C
		elem_st_mat(2, 0) += dt_div_2 * (cur_gp1->N1 * cur_gp1->dN1_dx * cur_gp1->dx_dxi + cur_gp2->N1 * cur_gp2->dN1_dx * cur_gp2->dx_dxi);
		elem_st_mat(2, 1) += dt_div_2 * (cur_gp1->N1 * cur_gp1->dN2_dx * cur_gp1->dx_dxi + cur_gp2->N1 * cur_gp2->dN2_dx * cur_gp2->dx_dxi);
		elem_st_mat(3, 0) += dt_div_2 * (cur_gp1->N2 * cur_gp1->dN1_dx * cur_gp1->dx_dxi + cur_gp2->N2 * cur_gp2->dN1_dx * cur_gp2->dx_dxi);
		elem_st_mat(3, 1) += dt_div_2 * (cur_gp1->N2 * cur_gp1->dN2_dx * cur_gp1->dx_dxi + cur_gp2->N2 * cur_gp2->dN2_dx * cur_gp2->dx_dxi);
		//std::cout << "elemental stiffness matrix: " << std::endl << elem_st_mat << std::endl;
		// + 0.25 * dt * dt * K
		elem_st_mat(0, 0) += dt_square_div_4 * (cur_gp1->dN1_dx * cur_elem->E * cur_gp1->dN1_dx * cur_gp1->dx_dxi
			                                  + cur_gp2->dN1_dx * cur_elem->E * cur_gp2->dN1_dx * cur_gp2->dx_dxi);
		elem_st_mat(0, 1) += dt_square_div_4 * (cur_gp1->dN1_dx * cur_elem->E * cur_gp1->dN2_dx * cur_gp1->dx_dxi
			                                  + cur_gp2->dN1_dx * cur_elem->E * cur_gp2->dN2_dx * cur_gp2->dx_dxi);
		elem_st_mat(0, 2) += dt_square_div_4 * -(cur_gp1->dN1_dx * cur_gp1->N1 * cur_gp1->dx_dxi + cur_gp2->dN1_dx * cur_gp2->N1 * cur_gp2->dx_dxi);
		elem_st_mat(0, 3) += dt_square_div_4 * -(cur_gp1->dN1_dx * cur_gp1->N2 * cur_gp1->dx_dxi + cur_gp2->dN1_dx * cur_gp2->N2 * cur_gp2->dx_dxi);
		elem_st_mat(1, 0) += dt_square_div_4 * (cur_gp1->dN2_dx * cur_elem->E * cur_gp1->dN1_dx * cur_gp1->dx_dxi
			                                  + cur_gp2->dN2_dx * cur_elem->E * cur_gp2->dN1_dx * cur_gp2->dx_dxi);
		elem_st_mat(1, 1) += dt_square_div_4 * (cur_gp1->dN2_dx * cur_elem->E * cur_gp1->dN2_dx * cur_gp1->dx_dxi
			                                  + cur_gp2->dN2_dx * cur_elem->E * cur_gp2->dN2_dx * cur_gp2->dx_dxi);
		elem_st_mat(1, 2) += dt_square_div_4 * -(cur_gp1->dN2_dx * cur_gp1->N1 * cur_gp1->dx_dxi + cur_gp2->dN2_dx * cur_gp2->N1 * cur_gp2->dx_dxi);
		elem_st_mat(1, 3) += dt_square_div_4 * -(cur_gp1->dN2_dx * cur_gp1->N2 * cur_gp1->dx_dxi + cur_gp2->dN2_dx * cur_gp2->N2 * cur_gp2->dx_dxi);
		elem_st_mat(2, 2) += dt_square_div_4 * cur_elem->k_div_miu * (cur_gp1->dN1_dx * cur_gp1->dN1_dx * cur_gp1->dx_dxi + cur_gp2->dN1_dx * cur_gp2->dN1_dx * cur_gp2->dx_dxi);
		elem_st_mat(2, 3) += dt_square_div_4 * cur_elem->k_div_miu * (cur_gp1->dN1_dx * cur_gp1->dN2_dx * cur_gp1->dx_dxi + cur_gp2->dN1_dx * cur_gp2->dN2_dx * cur_gp2->dx_dxi);
		elem_st_mat(3, 2) += dt_square_div_4 * cur_elem->k_div_miu * (cur_gp1->dN2_dx * cur_gp1->dN1_dx * cur_gp1->dx_dxi + cur_gp2->dN2_dx * cur_gp2->dN1_dx * cur_gp2->dx_dxi);
		elem_st_mat(3, 3) += dt_square_div_4 * cur_elem->k_div_miu * (cur_gp1->dN2_dx * cur_gp1->dN2_dx * cur_gp1->dx_dxi + cur_gp2->dN2_dx * cur_gp2->dN2_dx * cur_gp2->dx_dxi);
		//std::cout << "elemental stiffness matrix: " << std::endl << elem_st_mat << std::endl;

		// Init elemental force vector
		elem_f_vec.setZero();
		// - C(v+ 0.5 * dt * a)
		elem_f_vec[2] -= (cur_gp1->N1 * cur_gp1->dN1_dx * cur_gp1->dx_dxi + cur_gp2->N1 * cur_gp2->dN1_dx * cur_gp2->dx_dxi) * tmp_coef3[cur_g_id[0]]  // node1
			           + (cur_gp1->N1 * cur_gp1->dN2_dx * cur_gp1->dx_dxi + cur_gp2->N1 * cur_gp2->dN2_dx * cur_gp2->dx_dxi) * tmp_coef3[cur_g_id[1]]; // node2
		elem_f_vec[3] -= (cur_gp1->N2 * cur_gp1->dN1_dx * cur_gp1->dx_dxi + cur_gp2->N2 * cur_gp2->dN1_dx * cur_gp2->dx_dxi) * tmp_coef3[cur_g_id[0]]  // node1
			           + (cur_gp1->N2 * cur_gp1->dN2_dx * cur_gp1->dx_dxi + cur_gp2->N2 * cur_gp2->dN2_dx * cur_gp2->dx_dxi) * tmp_coef3[cur_g_id[1]]; // node2
		//std::cout << "force vector matrix: " << std::endl << elem_f_vec << std::endl;
		// -N (dn)
		elem_f_vec[0] -= cur_gp1->dN1_dx * (cur_gp1->s11 - cur_gp1->N1 * cur_node1->p - cur_gp1->N2 * cur_node2->p) * cur_gp1->dx_dxi
			           + cur_gp2->dN1_dx * (cur_gp2->s11 - cur_gp2->N1 * cur_node1->p - cur_gp2->N2 * cur_node2->p) * cur_gp2->dx_dxi;
		elem_f_vec[1] -= cur_gp1->dN2_dx * (cur_gp1->s11 - cur_gp1->N1 * cur_node1->p - cur_gp1->N2 * cur_node2->p) * cur_gp1->dx_dxi
			           + cur_gp2->dN2_dx * (cur_gp2->s11 - cur_gp2->N1 * cur_node1->p - cur_gp2->N2 * cur_node2->p) * cur_gp2->dx_dxi;
		elem_f_vec[2] -= cur_elem->k_div_miu * (cur_gp1->dN1_dx * (cur_gp1->dN1_dx * cur_node1->p + cur_gp1->dN2_dx * cur_node2->p) * cur_gp1->dx_dxi
			                                  + cur_gp2->dN1_dx * (cur_gp2->dN1_dx * cur_node1->p + cur_gp2->dN2_dx * cur_node2->p) * cur_gp2->dx_dxi);
		elem_f_vec[3] -= cur_elem->k_div_miu * (cur_gp1->dN2_dx * (cur_gp1->dN1_dx * cur_node1->p + cur_gp1->dN2_dx * cur_node2->p) * cur_gp1->dx_dxi
			                                  + cur_gp2->dN2_dx * (cur_gp2->dN1_dx * cur_node1->p + cur_gp2->dN2_dx * cur_node2->p) * cur_gp2->dx_dxi);
		//std::cout << "force vector matrix: " << std::endl << elem_f_vec << std::endl;
		// - dN(v * dt + 0.25 * dt * dt * a)
		double cur_gp1_ds11 = cur_elem->E * (cur_gp1->dN1_dx * tmp_coef4[cur_g_id[0]] + cur_gp1->dN2_dx * tmp_coef4[cur_g_id[1]]);
		double cur_gp2_ds11 = cur_elem->E * (cur_gp2->dN1_dx * tmp_coef4[cur_g_id[0]] + cur_gp2->dN2_dx * tmp_coef4[cur_g_id[1]]);
		double cur_node1_dp = tmp_coef4[cur_g_id[2]];
		double cur_node2_dp = tmp_coef4[cur_g_id[3]];
		elem_f_vec[0] -= cur_gp1->dN1_dx * (cur_gp1_ds11 - cur_gp1->N1 * cur_node1_dp - cur_gp1->N2 * cur_node2_dp) * cur_gp1->dx_dxi
			           + cur_gp2->dN1_dx * (cur_gp2_ds11 - cur_gp2->N1 * cur_node1_dp - cur_gp2->N2 * cur_node2_dp) * cur_gp2->dx_dxi;
		elem_f_vec[1] -= cur_gp1->dN2_dx * (cur_gp1_ds11 - cur_gp1->N1 * cur_node1_dp - cur_gp1->N2 * cur_node2_dp) * cur_gp1->dx_dxi
			           + cur_gp2->dN2_dx * (cur_gp2_ds11 - cur_gp2->N1 * cur_node1_dp - cur_gp2->N2 * cur_node2_dp) * cur_gp2->dx_dxi;
		elem_f_vec[2] -= cur_elem->k_div_miu * (cur_gp1->dN1_dx * (cur_gp1->dN1_dx * cur_node1_dp + cur_gp1->dN2_dx * cur_node2_dp) * cur_gp1->dx_dxi
			                                  + cur_gp2->dN1_dx * (cur_gp2->dN1_dx * cur_node1_dp + cur_gp2->dN2_dx * cur_node2_dp) * cur_gp2->dx_dxi);
		elem_f_vec[3] -= cur_elem->k_div_miu * (cur_gp1->dN2_dx * (cur_gp1->dN1_dx * cur_node1_dp + cur_gp1->dN2_dx * cur_node2_dp) * cur_gp1->dx_dxi
			                                  + cur_gp2->dN2_dx * (cur_gp2->dN1_dx * cur_node1_dp + cur_gp2->dN2_dx * cur_node2_dp) * cur_gp2->dx_dxi);
		//std::cout << "force vector matrix: " << std::endl << elem_f_vec << std::endl;

		// Map elemental stiffness matrix to global stiffness matrix.
#define ElementalStiffnessMatrixToGlobalStiffnessMatrix(id_x, id_y) \
		g_st_mat_coef.push_back(Eigen::Triplet<double>(cur_g_id[id_x], cur_g_id[id_y], elem_st_mat(id_x, id_y)))
		ElementalStiffnessMatrixToGlobalStiffnessMatrix(0, 0);
		ElementalStiffnessMatrixToGlobalStiffnessMatrix(0, 1);
		ElementalStiffnessMatrixToGlobalStiffnessMatrix(0, 2);
		ElementalStiffnessMatrixToGlobalStiffnessMatrix(0, 3);
		ElementalStiffnessMatrixToGlobalStiffnessMatrix(1, 0);
		ElementalStiffnessMatrixToGlobalStiffnessMatrix(1, 1);
		ElementalStiffnessMatrixToGlobalStiffnessMatrix(1, 2);
		ElementalStiffnessMatrixToGlobalStiffnessMatrix(1, 3);
		ElementalStiffnessMatrixToGlobalStiffnessMatrix(2, 0);
		ElementalStiffnessMatrixToGlobalStiffnessMatrix(2, 1);
		ElementalStiffnessMatrixToGlobalStiffnessMatrix(2, 2);
		ElementalStiffnessMatrixToGlobalStiffnessMatrix(2, 3);
		ElementalStiffnessMatrixToGlobalStiffnessMatrix(3, 0);
		ElementalStiffnessMatrixToGlobalStiffnessMatrix(3, 1);
		ElementalStiffnessMatrixToGlobalStiffnessMatrix(3, 2);
		ElementalStiffnessMatrixToGlobalStiffnessMatrix(3, 3);

		// Map elemental force vector to global force vector.
#define ElementalForceVectorToGlobalForceVector(id) g_f_vec[cur_g_id[id]] += elem_f_vec[id]
		ElementalForceVectorToGlobalForceVector(0);
		ElementalForceVectorToGlobalForceVector(1);
		ElementalForceVectorToGlobalForceVector(2);
		ElementalForceVectorToGlobalForceVector(3);
	}

	//g_st_mat.setFromTriplets(g_st_mat_coef.begin(), g_st_mat_coef.end());
	//g_st_mat.makeCompressed();
	//std::cout << "global stiffness matrix" << std::endl << g_st_mat << std::endl;
	//std::cout << "global force vector" << std::endl << g_f_vec << std::endl;

	// Apply Neumann BC (directly modify the total force vector)
	double body_force_tmp;
	for (size_t i = 0; i < body_force_num; i++)
	{
		body_force_tmp = body_force[i].bf_x;
		cur_elem = elems + body_force[i].elem_id;;
		cur_node1 = nodes + cur_elem->index_x;
		cur_node2 = nodes + cur_elem->index_x + 1;
		cur_gp1 = &(cur_elem->gp1);
		cur_gp2 = &(cur_elem->gp2);

		// Init index mapping array
		cur_g_id[0] = cur_elem->index_x;
		cur_g_id[1] = cur_elem->index_x + 1;
		cur_g_id[2] = node_num + cur_elem->index_x;
		cur_g_id[3] = node_num + cur_elem->index_x + 1;
		
		// Form external force vector for this element
		elem_f_vec.setZero();
		elem_f_vec[0] = cur_gp1->N1 * body_force_tmp * cur_gp1->dx_dxi + cur_gp2->N2 * body_force_tmp * cur_gp2->dx_dxi;
		elem_f_vec[1] = cur_gp1->N2 * body_force_tmp * cur_gp1->dx_dxi + cur_gp2->N2 * body_force_tmp * cur_gp2->dx_dxi;
		elem_f_vec[2] = cur_elem->density_f * cur_elem->k_div_miu * (cur_gp1->dN1_dx * body_force_tmp * cur_gp1->dx_dxi + cur_gp2->dN1_dx * body_force_tmp * cur_gp2->dx_dxi);
		elem_f_vec[3] = cur_elem->density_f * cur_elem->k_div_miu * (cur_gp1->dN2_dx * body_force_tmp * cur_gp1->dx_dxi + cur_gp2->dN2_dx * body_force_tmp * cur_gp2->dx_dxi);

		// Map to global force vector
		ElementalForceVectorToGlobalForceVector(0);
		ElementalForceVectorToGlobalForceVector(1);
		ElementalForceVectorToGlobalForceVector(2);
		ElementalForceVectorToGlobalForceVector(3);
		//std::cout << "elemental force vector matrix: " << std::endl << elem_f_vec << std::endl;
	}
	//std::cout << "global force vector matrix: " << std::endl << g_f_vec << std::endl;

	double xi_tmp;
	for (size_t i = 0; i < surface_force_num; i++)
	{
		cur_elem = elems + surface_force[i].elem_id;;

		// Init index mapping array
		cur_g_id[0] = cur_elem->index_x;
		cur_g_id[1] = cur_elem->index_x + 1;
		cur_g_id[2] = node_num + cur_elem->index_x;
		cur_g_id[3] = node_num + cur_elem->index_x + 1;

		// Form external force vector for this element
		elem_f_vec.setZero();
		calNaturalCoord(cur_elem, &surface_force[i].x, &xi_tmp);
		elem_f_vec[0] = _N1(xi_tmp) * surface_force[i].sf_x;
		elem_f_vec[1] = _N2(xi_tmp) * surface_force[i].sf_x;

		// Map to global force vector
		ElementalForceVectorToGlobalForceVector(0);
		ElementalForceVectorToGlobalForceVector(1);
		ElementalForceVectorToGlobalForceVector(2);
		ElementalForceVectorToGlobalForceVector(3);
		//std::cout << "elemental force vector matrix: " << std::endl << elem_f_vec << std::endl;
	}
	// currently, the s * w type BC is ignored.
	//std::cout << "global force vector matrix: " << std::endl << g_f_vec << std::endl;

	/** Apply Dirichlet BC into force vector and stiffness matrix. **/
	// (The performance of this part can be enhanced.)
	// Acceleration
	size_t g_id_tmp;
	double Knn_tmp;
	for (size_t i = 0; i < acceleration_num; i++)
	{
		g_id_tmp = acceleration_bc[i].node_id;
		// modify stiffness matrix
		for (std::list<Eigen::Triplet<double> >::iterator iter = g_st_mat_coef.begin();
			iter != g_st_mat_coef.end();)
		{
			if (g_id_tmp == iter->row() || g_id_tmp == iter->col())
			{
				if (iter->row() == iter->col())
				{
					// diagonal term
					Knn_tmp = iter->value();
				}
				else
				{
					// set all non-diagonal term to zero
					iter = g_st_mat_coef.erase(iter);
					continue;
				}
			}
			iter++;
		}
		// modify force vector
		g_f_vec[g_id_tmp] = Knn_tmp * acceleration_bc[i].a;
	}
	// Pore pressure
	for (size_t i = 0; i < p_acceleration_num; i++)
	{
		g_id_tmp = node_num + p_acceleration_bc[i].node_id;
		// modify stiffness matrix
		for (std::list<Eigen::Triplet<double> >::iterator iter = g_st_mat_coef.begin();
			iter != g_st_mat_coef.end();)
		{
			if (g_id_tmp == iter->row() || g_id_tmp == iter->col())
			{
				if (iter->row() == iter->col())
				{
					// diagonal term
					Knn_tmp = iter->value();
				}
				else
				{
					// set all non-diagonal term to zero
					iter = g_st_mat_coef.erase(iter);
					continue;
				}
			}
			iter++;
		}
		// modify force vector
		g_f_vec[g_id_tmp] = Knn_tmp * p_acceleration_bc[i].pa;
	}

	g_st_mat.setFromTriplets(g_st_mat_coef.begin(), g_st_mat_coef.end());
	g_st_mat.makeCompressed();
	//std::cout << "global stiffness matrix" << std::endl << g_st_mat << std::endl;

	//std::cout << "global force vector" << std::endl << g_f_vec << std::endl;

	// Solve for acceleration
	Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<Eigen::Index> > solver;
	solver.analyzePattern(g_st_mat);
	solver.factorize(g_st_mat);
	g_a_vec = solver.solve(g_f_vec);

	//std::cout << "global acceleration vector" << std::endl << g_a_vec << std::endl;

	/** Reapply Dirichlet boundary conditions (may not be needed?) **/
	/*
	for (size_t i = 0; i< acceleration_num; i++)
	{
		g_id_tmp = acceleration_bc[i].node_id;
		g_a_vec[g_id_tmp] = acceleration_bc[i].a;
	}
	for (size_t i = 0; i < p_acceleration_num; i++)
	{
		g_id_tmp = node_num + p_acceleration_bc[i].node_id;
		g_a_vec[g_id_tmp] = p_acceleration_bc[i].pa;
	}
	*/

	// Update velocity and displacement
	for (size_t i = 0; i < node_num; i++)
	{
		g_id_tmp = i;
		nodes[i].du = dt * nodes[i].v + dt_square_div_4 * (nodes[i].a + g_a_vec[g_id_tmp]);
		nodes[i].u += nodes[i].du;
		nodes[i].v += dt_div_2 * (nodes[i].a + g_a_vec[g_id_tmp]);
		nodes[i].a = g_a_vec[g_id_tmp];
		
		g_id_tmp = node_num + i;
		nodes[i].p += dt * nodes[i].pv + dt_square_div_4 * (nodes[i].pa + g_a_vec[g_id_tmp]);
		nodes[i].pv += dt_div_2 * (nodes[i].pa + g_a_vec[g_id_tmp]);
		nodes[i].pa = g_a_vec[g_id_tmp];
	}
	// output nodal data
	/*
	std::cout << "nodal displacement: " << std::endl;
	for (size_t i = 0; i < node_num; i++)
		std::cout << nodes[i].u << ",";
	std::cout << std::endl;
	std::cout << "nodal pore pressure: " << std::endl;
	for (size_t i = 0; i < node_num; i++)
		std::cout << nodes[i].p << ",";
	std::cout << std::endl;
	*/

	// Update variables on Gauss points
	double de11;
	double de_vol_tmp;
	for (size_t i = 0; i < elem_num; i++)
	{
		cur_elem = elems + i;
		cur_node1 = nodes + cur_elem->index_x;
		cur_node2 = nodes + cur_elem->index_x + 1;
		
		cur_gp1 = &elems[i].gp1;
		de11 = cur_node1->du * cur_gp1->dN1_dx + cur_node2->du * cur_gp1->dN2_dx;
		cur_gp1->e11 += de11;
		cur_gp1->s11 += de11 * elems[i].E;
		de_vol_tmp = de11;
		cur_gp1->n = (de_vol_tmp + cur_gp1->n) / (1.0 + de_vol_tmp);
		//std::cout << "n: " << cur_gp1->n << " de11: " << de11
		//	      << " s11: " << cur_gp1->s11 << std::endl;

		cur_gp2 = &elems[i].gp2;
		de11 = cur_node1->du * cur_gp2->dN1_dx + cur_node2->du * cur_gp2->dN2_dx;;
		cur_gp2->e11 += de11;
		cur_gp2->s11 += de11 * elems[i].E;
		de_vol_tmp = de11;
		cur_gp2->n = (de_vol_tmp + cur_gp2->n) / (1.0 + de_vol_tmp);
	
		// need to calculate w in the future
	}
	
	return 0;
}


int output(std::ofstream &out, double cur_time, Node *nodes, size_t node_num)
{
	out << cur_time << ",";
	// output displacement
	for (size_t i = 0; i < node_num; i++)
	{
		out << nodes[i].u << ",";
	}
	for (size_t i = 0; i < node_num; i++)
	{
		out << nodes[i].p << ",";
	}
	out << std::endl;

	return 0;
}

}