#include <cstddef>
#include <cstdlib>

#include <fstream>
#include <iostream>

#include "hyme1D_FEM2.h"

/*
Note: 
1. Quadratic shape function for displacement;
2. Linear shape function for pore pressure;
3. Reduced gauss integration for both fields.
*/

/*
Parameters of the Newmark - beta method
1. gamma: v;
2. beta: d;
*/
//#define gamma 0.5
//#define beta 0.25
#define gamma 0.6
#define beta 0.3025

namespace hyme1D_FEM2
{
	size_t node_s_num;
	Node_s *node_ss;
	size_t node_f_num;
	Node_f *node_fs;

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

	inline void naturalCoord_s(Element *elem, double *coord, double *na_coord)
	{
		size_t id_x = elem->index_s_x;
		double n1_x = node_ss[id_x].x;
		double n2_x = node_ss[id_x + 2].x;
		na_coord[0] = (2.0 * coord[0] - (n2_x + n1_x)) / (n2_x - n1_x);
	}

	inline void naturalCoord_f(Element *elem, double *coord, double *na_coord)
	{
		size_t id_x = elem->index_f_x;
		double n1_x = node_fs[id_x].x;
		double n2_x = node_fs[id_x + 1].x;
		na_coord[0] = (2.0 * coord[0] - (n2_x + n1_x)) / (n2_x - n1_x);
	}

	// For displacement
#define _Ns1(xi) (((xi)*(xi) - (xi)) / 2.0)
#define _Ns2(xi) (1.0 - (xi)*(xi))
#define _Ns3(xi) (((xi)*(xi) + (xi)) / 2.0)
#define _dNs1_dxi(xi) ((xi) - 0.5)
#define _dNs2_dxi(xi) (-2.0 * (xi))
#define _dNs3_dxi(xi) ((xi) + 0.5)
// For pore pressure
#define _Nf1(xi) ((1.0 - (xi))/2.0)
#define _Nf2(xi) ((1.0 + (xi))/2.0)
#define _dNf1_dxi(xi) (-0.5)
#define _dNf2_dxi(xi) (0.5)

	inline double dx_dxi_s(Element *elem, double xi)
	{
		size_t id_x = elem->index_s_x;
		//std::cout << _dNs1_dxi(xi) << " " << _dNs2_dxi(xi) << " " << _dNs3_dxi(xi) << std::endl;
		return node_ss[id_x].x     * _dNs1_dxi(xi)
			+ node_ss[id_x + 1].x * _dNs2_dxi(xi)
			+ node_ss[id_x + 2].x * _dNs3_dxi(xi);
	}

	inline double dx_dxi_f(Element *elem, double xi)
	{
		size_t id_x = elem->index_f_x;
		return node_fs[id_x].x     * _dNf1_dxi(xi)
			+ node_fs[id_x + 1].x * _dNf2_dxi(xi);
	}


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
			node_ss[acceleration_bc[i].node_id].a = acceleration_bc[i].a;
		}
		// Velocity 
		for (size_t i = 0; i < velocity_num; i++)
		{
			node_ss[velocity_bc[i].node_id].v = velocity_bc[i].v;
		}
		// Displacement 
		for (size_t i = 0; i < displacement_num; i++)
		{
			node_ss[displacement_bc[i].node_id].u = displacement_bc[i].u;
			// need to update stress here
		}
		// Acceleration of pore pressure
		for (size_t i = 0; i < p_acceleration_num; i++)
		{
			node_fs[p_acceleration_bc[i].node_id].pa = p_acceleration_bc[i].pa;
		}
		// Pore pressure
		for (size_t i = 0; i < pore_pressure_num; i++)
		{
			node_fs[pore_pressure_bc[i].node_id].p = pore_pressure_bc[i].p;
		}

		return 0;
	};

	// Calculate each substeps
	int substep(void);
	int output(std::ofstream &out);

	// calculate the whole step
	// simplify ver.
	int cal(std::ofstream &out)
	{
		while (cur_t < total_t)
		{
			substep();
			cur_t += dt;
			output(out);
		}

		return 0;
	}

	int hyme1D_FEM2(void)
	{
		/* init mesh data */
		elem_num = 10;
		node_s_num = 2 * elem_num + 1;
		node_f_num = elem_num + 1;
		elems = (Element *)malloc(sizeof(Element) * elem_num);
		node_ss = (Node_s *)malloc(sizeof(Node_s) * node_s_num);
		node_fs = (Node_f *)malloc(sizeof(Node_f) * node_f_num);
		for (size_t i = 0; i < node_s_num; i++)
		{
			node_ss[i].x = (double)i * 0.5;

			// default initialization
			node_ss[i].a = 0.0;
			node_ss[i].v = 0.0;
			node_ss[i].u = 0.0;
		}
		for (size_t i = 0; i < node_f_num; i++)
		{
			node_fs[i].x = (double)i;

			// default initialization
			node_fs[i].pa = 0.0;
			node_fs[i].pv = 0.0;
			node_fs[i].p = 0.0;
		}
		for (size_t i = 0; i < elem_num; i++)
		{
			elems[i].index_s_x = 2 * i;
			elems[i].index_f_x = i;

			elems[i].density_s = 3000.0;
			elems[i].density_f = 1000.0;

			elems[i].E = 40.0e6;
			elems[i].k = 1.0e-5;
			elems[i].miu = 1.0;

			// init gauss point
			elems[i].gps1.s11 = 0.0;
			elems[i].gps1.e11 = 0.0;
			elems[i].gps1.n = 0.5;
			elems[i].gps2.s11 = 0.0;
			elems[i].gps2.e11 = 0.0;
			elems[i].gps2.n = 0.5;

			// default initialization
			elems[i].k_div_miu = elems[i].k / elems[i].miu;
			elems[i].gps1.density_avg = (1.0 - elems[i].gps1.n) * elems[i].density_s
				+ elems[i].gps1.n  * elems[i].density_f;
			elems[i].gps2.density_avg = (1.0 - elems[i].gps2.n) * elems[i].density_s
				+ elems[i].gps2.n  * elems[i].density_f;
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
		surface_force[0].x = node_ss[node_s_num - 1].x;

		// Acceleration BC
		acceleration_num = 1;
		acceleration_bc = (AccelerationBC *)malloc(sizeof(AccelerationBC) * acceleration_num);
		acceleration_bc[0].a = 0.0;
		acceleration_bc[0].node_id = 0;

		// Pore pressure BC
		// free flow boundary
		//pore_pressure_num = 1;
		//pore_pressure_bc = (PorePressureBC *)malloc(sizeof(PorePressureBC) * pore_pressure_num);
		//pore_pressure_bc[0].p = 0.0;
		//pore_pressure_bc[0].node_id = node_f_num - 1;
		// impermeable boundary
		pore_pressure_num = 0;
		pore_pressure_bc = NULL;

		// Velocity BC 
		velocity_num = 0;
		velocity_bc = NULL;
		// Displacement BC
		displacement_num = 0;
		displacement_bc = NULL;

		std::ofstream out("res.txt");

		total_dof = node_s_num + node_f_num;
		// v + (1 - gamma) * a * dt
		tmp_coef3 = (double *)malloc(sizeof(double) * node_s_num);
		// v * dt + (0.5 - beta) * a * dt * dt
		tmp_coef4 = (double *)malloc(sizeof(double) * total_dof);

		// Initialize calculation
		init();

		total_t = 1.0;
		dt = 0.002;
		cal(out);

		// free flow boundary conditions
		pore_pressure_num = 1;
		pore_pressure_bc = (PorePressureBC *)malloc(sizeof(PorePressureBC) * pore_pressure_num);
		pore_pressure_bc[0].p = 0.0;
		pore_pressure_bc[0].node_id = node_f_num - 1;
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

		// Modify initial conditions
		// Acceleration of pore pressure
		for (size_t i = 0; i < p_acceleration_num; i++)
		{
			node_fs[p_acceleration_bc[i].node_id].pa = p_acceleration_bc[i].pa;
		}
		// Pore pressure
		for (size_t i = 0; i < pore_pressure_num; i++)
		{
			node_fs[pore_pressure_bc[i].node_id].p = pore_pressure_bc[i].p;
		}

		total_t = 3.0;
		dt = 0.002;
		cal(out);

		// Delete data
		if (node_ss) free(node_ss);
		if (node_fs) free(node_fs);
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
		double elem_st_mat[5][5];
		// Elemental force vector
		double elem_f_vec[5];
		Element *cur_elem;
		Node_s *cur_nodes1, *cur_nodes2, *cur_nodes3;
		Node_f *cur_nodef1, *cur_nodef2;
		GaussPoint *cur_gps1, *cur_gps2, *cur_gpf;
		/* Index mapping array
		used to map coefficient of elemental stiffness
		matrix into glabal stiffness matrix.
		*/
		size_t cur_g_id[5];

		// Temporary calculation param
		const static double gamma_dt = gamma * dt;
		const static double one_minus_gamma_dt = (1.0 - gamma) * dt;
		const static double beta_dt_square = beta * dt * dt;
		const static double point_five_minus_beta_dt_square = (0.5 - beta) * dt * dt;
		for (size_t i = 0; i < node_s_num; i++)
		{
			// v + (1 - gamma) * dt * a
			tmp_coef3[i] = node_ss[i].v + one_minus_gamma_dt * node_ss[i].a;
			// v * dt + (0.5 - beta) * dt * dt * a
			tmp_coef4[i] = dt * node_ss[i].v + point_five_minus_beta_dt_square * node_ss[i].a;
		}
		for (size_t i = 0; i < node_f_num; i++)
		{
			tmp_coef4[node_s_num + i] = dt * node_fs[i].pv + point_five_minus_beta_dt_square * node_fs[i].pa;
		}

		g_f_vec.setZero();

		// Calculate elemental stiffness matrix and force vector
		for (size_t elem_id = 0; elem_id < elem_num; elem_id++)
		{
			cur_elem = elems + elem_id;
			cur_nodes1 = node_ss + cur_elem->index_s_x;
			cur_nodes2 = node_ss + cur_elem->index_s_x + 1;
			cur_nodes3 = node_ss + cur_elem->index_s_x + 2;
			cur_gps1 = &(cur_elem->gps1);
			cur_gps2 = &(cur_elem->gps2);
			cur_nodef1 = node_fs + cur_elem->index_f_x;
			cur_nodef2 = node_fs + cur_elem->index_f_x + 1;
			cur_gpf = &(cur_elem->gpf);

			// init index mapping array
			cur_g_id[0] = cur_elem->index_s_x;
			cur_g_id[1] = cur_elem->index_s_x + 1;
			cur_g_id[2] = cur_elem->index_s_x + 2;
			cur_g_id[3] = node_s_num + cur_elem->index_f_x;
			cur_g_id[4] = node_s_num + cur_elem->index_f_x + 1;

			// init integration point
#define INIT_INTEGRATION_POINT(elem, gp, gp_coord, gp_w)        \
		{                                                       \
			(gp).xi = gp_coord;                                 \
			(gp).w = gp_w;                                      \
			/* solid phase */                                   \
			(gp).Ns1 = _Ns1((gp).xi);                           \
			(gp).Ns2 = _Ns2((gp).xi);                           \
			(gp).Ns3 = _Ns3((gp).xi);                           \
			(gp).dx_dxi_s = dx_dxi_s(&(elem), (gp).xi);         \
			(gp).dxi_dx_s = 1.0 / (gp).dx_dxi_s;                \
			(gp).dNs1_dx = _dNs1_dxi(gp_coord) * (gp).dxi_dx_s; \
			(gp).dNs2_dx = _dNs2_dxi(gp_coord) * (gp).dxi_dx_s; \
			(gp).dNs3_dx = _dNs3_dxi(gp_coord) * (gp).dxi_dx_s; \
			/* fluid phase */                                   \
			(gp).Nf1 = _Nf1((gp).xi);                           \
			(gp).Nf2 = _Nf2((gp).xi);                           \
			(gp).dx_dxi_f = dx_dxi_f(&(elem), (gp).xi);         \
			(gp).dxi_dx_f = 1.0 / (gp).dx_dxi_f;                \
			(gp).dNf1_dx = _dNf1_dxi(gp_coord) * (gp).dxi_dx_f; \
			(gp).dNf2_dx = _dNf2_dxi(gp_coord) * (gp).dxi_dx_f; \
		} while (false)
			INIT_INTEGRATION_POINT(*cur_elem, *cur_gps1, -0.57735, 1.0);
			INIT_INTEGRATION_POINT(*cur_elem, *cur_gps2, 0.57735, 1.0);
			INIT_INTEGRATION_POINT(*cur_elem, *cur_gpf, 0.0, 2.0);

			// Init elemental stiffness matrix
			memset(elem_st_mat, 0, sizeof(double) * 5 * 5);
			// M
			elem_st_mat[0][0] = cur_gps1->density_avg * cur_gps1->Ns1 * cur_gps1->Ns1 * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->density_avg * cur_gps2->Ns1 * cur_gps2->Ns1 * cur_gps2->dx_dxi_s * cur_gps2->w;
			elem_st_mat[0][1] = cur_gps1->density_avg * cur_gps1->Ns1 * cur_gps1->Ns2 * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->density_avg * cur_gps2->Ns1 * cur_gps2->Ns2 * cur_gps2->dx_dxi_s * cur_gps2->w;
			elem_st_mat[0][2] = cur_gps1->density_avg * cur_gps1->Ns1 * cur_gps1->Ns3 * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->density_avg * cur_gps2->Ns1 * cur_gps2->Ns3 * cur_gps2->dx_dxi_s * cur_gps2->w;
			elem_st_mat[1][0] = cur_gps1->density_avg * cur_gps1->Ns2 * cur_gps1->Ns1 * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->density_avg * cur_gps2->Ns2 * cur_gps2->Ns1 * cur_gps2->dx_dxi_s * cur_gps2->w;
			elem_st_mat[1][1] = cur_gps1->density_avg * cur_gps1->Ns2 * cur_gps1->Ns2 * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->density_avg * cur_gps2->Ns2 * cur_gps2->Ns2 * cur_gps2->dx_dxi_s * cur_gps2->w;
			elem_st_mat[1][2] = cur_gps1->density_avg * cur_gps1->Ns2 * cur_gps1->Ns3 * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->density_avg * cur_gps2->Ns2 * cur_gps2->Ns3 * cur_gps2->dx_dxi_s * cur_gps2->w;
			elem_st_mat[2][0] = cur_gps1->density_avg * cur_gps1->Ns3 * cur_gps1->Ns1 * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->density_avg * cur_gps2->Ns3 * cur_gps2->Ns1 * cur_gps2->dx_dxi_s * cur_gps2->w;
			elem_st_mat[2][1] = cur_gps1->density_avg * cur_gps1->Ns3 * cur_gps1->Ns2 * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->density_avg * cur_gps2->Ns3 * cur_gps2->Ns2 * cur_gps2->dx_dxi_s * cur_gps2->w;
			elem_st_mat[2][2] = cur_gps1->density_avg * cur_gps1->Ns3 * cur_gps1->Ns3 * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->density_avg * cur_gps2->Ns3 * cur_gps2->Ns3 * cur_gps2->dx_dxi_s * cur_gps2->w;
			// Terms can cause ossilliation in result 
			/*
			elem_st_mat[3][0] = cur_elem->density_f * cur_elem->k_div_miu * cur_gpf->dNf1_dx * cur_gpf->Ns1 * cur_gpf->dx_dxi_f * cur_gpf->w;
			elem_st_mat[3][1] = cur_elem->density_f * cur_elem->k_div_miu * cur_gpf->dNf1_dx * cur_gpf->Ns2 * cur_gpf->dx_dxi_f * cur_gpf->w;
			elem_st_mat[3][2] = cur_elem->density_f * cur_elem->k_div_miu * cur_gpf->dNf1_dx * cur_gpf->Ns3 * cur_gpf->dx_dxi_f * cur_gpf->w;
			elem_st_mat[4][0] = cur_elem->density_f * cur_elem->k_div_miu * cur_gpf->dNf2_dx * cur_gpf->Ns1 * cur_gpf->dx_dxi_f * cur_gpf->w;
			elem_st_mat[4][1] = cur_elem->density_f * cur_elem->k_div_miu * cur_gpf->dNf2_dx * cur_gpf->Ns2 * cur_gpf->dx_dxi_f * cur_gpf->w;
			elem_st_mat[4][2] = cur_elem->density_f * cur_elem->k_div_miu * cur_gpf->dNf2_dx * cur_gpf->Ns3 * cur_gpf->dx_dxi_f * cur_gpf->w;
			*/
			//dis_mat((double *)elem_st_mat, 5, 5);
			// + gamma * dt * C
			elem_st_mat[3][0] += gamma_dt * cur_gpf->Nf1 * cur_gpf->dNs1_dx * cur_gpf->dx_dxi_f * cur_gpf->w;
			elem_st_mat[3][1] += gamma_dt * cur_gpf->Nf1 * cur_gpf->dNs2_dx * cur_gpf->dx_dxi_f * cur_gpf->w;
			elem_st_mat[3][2] += gamma_dt * cur_gpf->Nf1 * cur_gpf->dNs3_dx * cur_gpf->dx_dxi_f * cur_gpf->w;
			elem_st_mat[4][0] += gamma_dt * cur_gpf->Nf2 * cur_gpf->dNs1_dx * cur_gpf->dx_dxi_f * cur_gpf->w;
			elem_st_mat[4][1] += gamma_dt * cur_gpf->Nf2 * cur_gpf->dNs2_dx * cur_gpf->dx_dxi_f * cur_gpf->w;
			elem_st_mat[4][2] += gamma_dt * cur_gpf->Nf2 * cur_gpf->dNs3_dx * cur_gpf->dx_dxi_f * cur_gpf->w;
			//dis_mat((double *)elem_st_mat, 5, 5);
			// + beta * dt * dt * K
			double tmp1;
			tmp1 = (cur_gps1->dNs1_dx * cur_elem->E * cur_gps1->dNs1_dx * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->dNs1_dx * cur_elem->E * cur_gps2->dNs1_dx * cur_gps2->dx_dxi_s * cur_gps2->w);
			elem_st_mat[0][0] += beta_dt_square * tmp1;
			tmp1 = (cur_gps1->dNs1_dx * cur_elem->E * cur_gps1->dNs2_dx * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->dNs1_dx * cur_elem->E * cur_gps2->dNs2_dx * cur_gps2->dx_dxi_s * cur_gps2->w);
			elem_st_mat[0][1] += beta_dt_square * tmp1;
			tmp1 = (cur_gps1->dNs1_dx * cur_elem->E * cur_gps1->dNs3_dx * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->dNs1_dx * cur_elem->E * cur_gps2->dNs3_dx * cur_gps2->dx_dxi_s * cur_gps2->w);
			elem_st_mat[0][2] += beta_dt_square * tmp1;
			tmp1 = -(cur_gps1->dNs1_dx * cur_gps1->Nf1 * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->dNs1_dx * cur_gps2->Nf1 * cur_gps2->dx_dxi_s * cur_gps2->w);
			elem_st_mat[0][3] += beta_dt_square * tmp1;
			tmp1 = -(cur_gps1->dNs1_dx * cur_gps1->Nf2 * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->dNs1_dx * cur_gps2->Nf2 * cur_gps2->dx_dxi_s * cur_gps2->w);
			elem_st_mat[0][4] += beta_dt_square * tmp1;

			tmp1 = (cur_gps1->dNs2_dx * cur_elem->E * cur_gps1->dNs1_dx * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->dNs2_dx * cur_elem->E * cur_gps2->dNs1_dx * cur_gps2->dx_dxi_s * cur_gps2->w);
			elem_st_mat[1][0] += beta_dt_square * tmp1;
			tmp1 = (cur_gps1->dNs2_dx * cur_elem->E * cur_gps1->dNs2_dx * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->dNs2_dx * cur_elem->E * cur_gps2->dNs2_dx * cur_gps2->dx_dxi_s * cur_gps2->w);
			elem_st_mat[1][1] += beta_dt_square * tmp1;
			tmp1 = (cur_gps1->dNs2_dx * cur_elem->E * cur_gps1->dNs3_dx * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->dNs2_dx * cur_elem->E * cur_gps2->dNs3_dx * cur_gps2->dx_dxi_s * cur_gps2->w);
			elem_st_mat[1][2] += beta_dt_square * tmp1;
			tmp1 = -(cur_gps1->dNs2_dx * cur_gps1->Nf1 * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->dNs2_dx * cur_gps2->Nf1 * cur_gps2->dx_dxi_s * cur_gps2->w);
			elem_st_mat[1][3] += beta_dt_square * tmp1;
			tmp1 = -(cur_gps1->dNs2_dx * cur_gps1->Nf2 * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->dNs2_dx * cur_gps2->Nf2 * cur_gps2->dx_dxi_s * cur_gps2->w);
			elem_st_mat[1][4] += beta_dt_square * tmp1;

			tmp1 = (cur_gps1->dNs3_dx * cur_elem->E * cur_gps1->dNs1_dx * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->dNs3_dx * cur_elem->E * cur_gps2->dNs1_dx * cur_gps2->dx_dxi_s * cur_gps2->w);
			elem_st_mat[2][0] += beta_dt_square * tmp1;
			tmp1 = (cur_gps1->dNs3_dx * cur_elem->E * cur_gps1->dNs2_dx * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->dNs3_dx * cur_elem->E * cur_gps2->dNs2_dx * cur_gps2->dx_dxi_s * cur_gps2->w);
			elem_st_mat[2][1] += beta_dt_square * tmp1;
			tmp1 = (cur_gps1->dNs3_dx * cur_elem->E * cur_gps1->dNs3_dx * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->dNs3_dx * cur_elem->E * cur_gps2->dNs3_dx * cur_gps2->dx_dxi_s * cur_gps2->w);
			elem_st_mat[2][2] += beta_dt_square * tmp1;
			tmp1 = -(cur_gps1->dNs3_dx * cur_gps1->Nf1 * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->dNs3_dx * cur_gps2->Nf1 * cur_gps2->dx_dxi_s * cur_gps2->w);
			elem_st_mat[2][3] += beta_dt_square * tmp1;
			tmp1 = -(cur_gps1->dNs3_dx * cur_gps1->Nf2 * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->dNs3_dx * cur_gps2->Nf2 * cur_gps2->dx_dxi_s * cur_gps2->w);
			elem_st_mat[2][4] += beta_dt_square * tmp1;

			tmp1 = cur_elem->k_div_miu * cur_gpf->dNf1_dx * cur_gpf->dNf1_dx * cur_gpf->dx_dxi_f * cur_gpf->w;
			elem_st_mat[3][3] += beta_dt_square * tmp1;
			tmp1 = cur_elem->k_div_miu * cur_gpf->dNf1_dx * cur_gpf->dNf2_dx * cur_gpf->dx_dxi_f * cur_gpf->w;
			elem_st_mat[3][4] += beta_dt_square * tmp1;
			tmp1 = cur_elem->k_div_miu * cur_gpf->dNf2_dx * cur_gpf->dNf1_dx * cur_gpf->dx_dxi_f * cur_gpf->w;
			elem_st_mat[4][3] += beta_dt_square * tmp1;
			tmp1 = cur_elem->k_div_miu * cur_gpf->dNf2_dx * cur_gpf->dNf2_dx * cur_gpf->dx_dxi_f * cur_gpf->w;
			elem_st_mat[4][4] += beta_dt_square * tmp1;

			//dis_mat((double *)elem_st_mat, 5, 5);

			// Init elemental force vector
			memset(elem_f_vec, 0, sizeof(double) * 5);
			// - C(v+ 0.5 * dt * a)
			elem_f_vec[3] -= cur_gpf->Nf1 * cur_gpf->dx_dxi_f * cur_gpf->w *
				(cur_gpf->dNs1_dx * tmp_coef3[cur_g_id[0]]    // node1
					+ cur_gpf->dNs2_dx * tmp_coef3[cur_g_id[1]]    // node2
					+ cur_gpf->dNs3_dx * tmp_coef3[cur_g_id[2]]);  // node3
			elem_f_vec[4] -= cur_gpf->Nf2 * cur_gpf->dx_dxi_f * cur_gpf->w *
				(cur_gpf->dNs1_dx * tmp_coef3[cur_g_id[0]]    // node1
					+ cur_gpf->dNs2_dx * tmp_coef3[cur_g_id[1]]    // node2
					+ cur_gpf->dNs3_dx * tmp_coef3[cur_g_id[2]]);  // node3
	 //std::cout << "force vector matrix: " << std::endl << elem_f_vec << std::endl;
	 // -N (dn)
			elem_f_vec[0] -= cur_gps1->dNs1_dx * (cur_gps1->s11 - cur_gps1->Nf1 * cur_nodef1->p - cur_gps1->Nf2 * cur_nodef2->p) * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->dNs1_dx * (cur_gps2->s11 - cur_gps2->Nf1 * cur_nodef1->p - cur_gps2->Nf2 * cur_nodef2->p) * cur_gps2->dx_dxi_s * cur_gps2->w;
			elem_f_vec[1] -= cur_gps1->dNs2_dx * (cur_gps1->s11 - cur_gps1->Nf1 * cur_nodef1->p - cur_gps1->Nf2 * cur_nodef2->p) * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->dNs2_dx * (cur_gps2->s11 - cur_gps2->Nf1 * cur_nodef1->p - cur_gps2->Nf2 * cur_nodef2->p) * cur_gps2->dx_dxi_s * cur_gps2->w;
			elem_f_vec[2] -= cur_gps1->dNs3_dx * (cur_gps1->s11 - cur_gps1->Nf1 * cur_nodef1->p - cur_gps1->Nf2 * cur_nodef2->p) * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->dNs3_dx * (cur_gps2->s11 - cur_gps2->Nf1 * cur_nodef1->p - cur_gps2->Nf2 * cur_nodef2->p) * cur_gps2->dx_dxi_s * cur_gps2->w;
			elem_f_vec[3] -= cur_elem->k_div_miu * cur_gpf->dNf1_dx * (cur_gpf->dNf1_dx * cur_nodef1->p + cur_gpf->dNf2_dx * cur_nodef2->p) * cur_gpf->dx_dxi_f * cur_gpf->w;
			elem_f_vec[4] -= cur_elem->k_div_miu * cur_gpf->dNf2_dx * (cur_gpf->dNf1_dx * cur_nodef1->p + cur_gpf->dNf2_dx * cur_nodef2->p) * cur_gpf->dx_dxi_f * cur_gpf->w;
			//std::cout << "force vector matrix: " << std::endl << elem_f_vec << std::endl;
			// - dN(v * dt + 0.25 * dt * dt * a)
			double gps1_ds11 = cur_elem->E * (cur_gps1->dNs1_dx * tmp_coef4[cur_g_id[0]]
				+ cur_gps1->dNs2_dx * tmp_coef4[cur_g_id[1]]
				+ cur_gps1->dNs3_dx * tmp_coef4[cur_g_id[2]]);
			double gps2_ds11 = cur_elem->E * (cur_gps2->dNs1_dx * tmp_coef4[cur_g_id[0]]
				+ cur_gps2->dNs2_dx * tmp_coef4[cur_g_id[1]]
				+ cur_gps2->dNs3_dx * tmp_coef4[cur_g_id[2]]);
			double cur_node1_dp = tmp_coef4[cur_g_id[3]];
			double cur_node2_dp = tmp_coef4[cur_g_id[4]];
			elem_f_vec[0] -= cur_gps1->dNs1_dx * (gps1_ds11 - cur_gps1->Nf1 * cur_node1_dp - cur_gps1->Nf2 * cur_node2_dp) * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->dNs1_dx * (gps2_ds11 - cur_gps2->Nf1 * cur_node1_dp - cur_gps2->Nf2 * cur_node2_dp) * cur_gps2->dx_dxi_s * cur_gps2->w;
			elem_f_vec[1] -= cur_gps1->dNs2_dx * (gps1_ds11 - cur_gps1->Nf1 * cur_node1_dp - cur_gps1->Nf2 * cur_node2_dp) * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->dNs2_dx * (gps2_ds11 - cur_gps2->Nf1 * cur_node1_dp - cur_gps2->Nf2 * cur_node2_dp) * cur_gps2->dx_dxi_s * cur_gps2->w;
			elem_f_vec[2] -= cur_gps1->dNs3_dx * (gps1_ds11 - cur_gps1->Nf1 * cur_node1_dp - cur_gps1->Nf2 * cur_node2_dp) * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->dNs3_dx * (gps2_ds11 - cur_gps2->Nf1 * cur_node1_dp - cur_gps2->Nf2 * cur_node2_dp) * cur_gps2->dx_dxi_s * cur_gps2->w;
			elem_f_vec[3] -= cur_elem->k_div_miu * cur_gpf->dNf1_dx * (cur_gpf->dNf1_dx * cur_node1_dp + cur_gpf->dNf2_dx * cur_node2_dp) * cur_gpf->dx_dxi_f * cur_gpf->w;
			elem_f_vec[4] -= cur_elem->k_div_miu * cur_gpf->dNf2_dx * (cur_gpf->dNf1_dx * cur_node1_dp + cur_gpf->dNf2_dx * cur_node2_dp) * cur_gpf->dx_dxi_f * cur_gpf->w;

			//dis_vec(elem_f_vec, 5);

			// Map elemental stiffness matrix to global stiffness matrix.
#define ElementalStiffnessMatrixToGlobalStiffnessMatrix(id_x, id_y) \
		g_st_mat_coef.push_back(Eigen::Triplet<double>(cur_g_id[id_x], cur_g_id[id_y], elem_st_mat[id_x][id_y]))
			ElementalStiffnessMatrixToGlobalStiffnessMatrix(0, 0);
			ElementalStiffnessMatrixToGlobalStiffnessMatrix(0, 1);
			ElementalStiffnessMatrixToGlobalStiffnessMatrix(0, 2);
			ElementalStiffnessMatrixToGlobalStiffnessMatrix(0, 3);
			ElementalStiffnessMatrixToGlobalStiffnessMatrix(0, 4);
			ElementalStiffnessMatrixToGlobalStiffnessMatrix(1, 0);
			ElementalStiffnessMatrixToGlobalStiffnessMatrix(1, 1);
			ElementalStiffnessMatrixToGlobalStiffnessMatrix(1, 2);
			ElementalStiffnessMatrixToGlobalStiffnessMatrix(1, 3);
			ElementalStiffnessMatrixToGlobalStiffnessMatrix(1, 4);
			ElementalStiffnessMatrixToGlobalStiffnessMatrix(2, 0);
			ElementalStiffnessMatrixToGlobalStiffnessMatrix(2, 1);
			ElementalStiffnessMatrixToGlobalStiffnessMatrix(2, 2);
			ElementalStiffnessMatrixToGlobalStiffnessMatrix(2, 3);
			ElementalStiffnessMatrixToGlobalStiffnessMatrix(2, 4);
			ElementalStiffnessMatrixToGlobalStiffnessMatrix(3, 0);
			ElementalStiffnessMatrixToGlobalStiffnessMatrix(3, 1);
			ElementalStiffnessMatrixToGlobalStiffnessMatrix(3, 2);
			ElementalStiffnessMatrixToGlobalStiffnessMatrix(3, 3);
			ElementalStiffnessMatrixToGlobalStiffnessMatrix(3, 4);
			ElementalStiffnessMatrixToGlobalStiffnessMatrix(4, 0);
			ElementalStiffnessMatrixToGlobalStiffnessMatrix(4, 1);
			ElementalStiffnessMatrixToGlobalStiffnessMatrix(4, 2);
			ElementalStiffnessMatrixToGlobalStiffnessMatrix(4, 3);
			ElementalStiffnessMatrixToGlobalStiffnessMatrix(4, 4);

			// Map elemental force vector to global force vector.
#define ElementalForceVectorToGlobalForceVector(id) g_f_vec[cur_g_id[id]] += elem_f_vec[id]
			ElementalForceVectorToGlobalForceVector(0);
			ElementalForceVectorToGlobalForceVector(1);
			ElementalForceVectorToGlobalForceVector(2);
			ElementalForceVectorToGlobalForceVector(3);
			ElementalForceVectorToGlobalForceVector(4);
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
			cur_gps1 = &(cur_elem->gps1);
			cur_gps2 = &(cur_elem->gps2);
			cur_gpf = &(cur_elem->gpf);

			// Init index mapping array
			cur_g_id[0] = cur_elem->index_s_x;
			cur_g_id[1] = cur_elem->index_s_x + 1;
			cur_g_id[2] = cur_elem->index_s_x + 2;
			cur_g_id[3] = node_s_num + cur_elem->index_f_x;
			cur_g_id[4] = node_s_num + cur_elem->index_f_x + 1;

			// Form external force vector for this element
			memset(elem_f_vec, 0, sizeof(double) * 5);
			elem_f_vec[0] = cur_gps1->density_avg * cur_gps1->Ns1 * body_force_tmp * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->density_avg * cur_gps2->Ns1 * body_force_tmp * cur_gps2->dx_dxi_s * cur_gps2->w;
			elem_f_vec[1] = cur_gps1->density_avg * cur_gps1->Ns2 * body_force_tmp * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->density_avg * cur_gps2->Ns2 * body_force_tmp * cur_gps2->dx_dxi_s * cur_gps2->w;
			elem_f_vec[2] = cur_gps1->density_avg * cur_gps1->Ns3 * body_force_tmp * cur_gps1->dx_dxi_s * cur_gps1->w
				+ cur_gps2->density_avg * cur_gps2->Ns3 * body_force_tmp * cur_gps2->dx_dxi_s * cur_gps2->w;
			elem_f_vec[3] = cur_elem->density_f * cur_elem->k_div_miu * cur_gpf->dNf1_dx * body_force_tmp * cur_gpf->dx_dxi_f * cur_gpf->w;
			elem_f_vec[4] = cur_elem->density_f * cur_elem->k_div_miu * cur_gpf->dNf2_dx * body_force_tmp * cur_gpf->dx_dxi_f * cur_gpf->w;

			// Map to global force vector
			ElementalForceVectorToGlobalForceVector(0);
			ElementalForceVectorToGlobalForceVector(1);
			ElementalForceVectorToGlobalForceVector(2);
			ElementalForceVectorToGlobalForceVector(3);
			ElementalForceVectorToGlobalForceVector(4);
			//std::cout << "elemental force vector matrix: " << std::endl << elem_f_vec << std::endl;
		}
		//std::cout << "global force vector matrix: " << std::endl << g_f_vec << std::endl;

		double xi_tmp;
		for (size_t i = 0; i < surface_force_num; i++)
		{
			cur_elem = elems + surface_force[i].elem_id;;

			// Init index mapping array
			cur_g_id[0] = cur_elem->index_s_x;
			cur_g_id[1] = cur_elem->index_s_x + 1;
			cur_g_id[2] = cur_elem->index_s_x + 2;
			cur_g_id[3] = node_s_num + cur_elem->index_f_x;
			cur_g_id[4] = node_s_num + cur_elem->index_f_x + 1;

			// Form external force vector for this element
			memset(elem_f_vec, 0, sizeof(double) * 5);
			naturalCoord_s(cur_elem, &surface_force[i].x, &xi_tmp);
			elem_f_vec[0] = _Ns1(xi_tmp) * surface_force[i].sf_x;
			elem_f_vec[1] = _Ns2(xi_tmp) * surface_force[i].sf_x;
			elem_f_vec[2] = _Ns3(xi_tmp) * surface_force[i].sf_x;

			// Map to global force vector
			ElementalForceVectorToGlobalForceVector(0);
			ElementalForceVectorToGlobalForceVector(1);
			ElementalForceVectorToGlobalForceVector(2);
			ElementalForceVectorToGlobalForceVector(3);
			ElementalForceVectorToGlobalForceVector(4);
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
			g_id_tmp = node_s_num + p_acceleration_bc[i].node_id;
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
		for (size_t i = 0; i < node_s_num; i++)
		{
			g_id_tmp = i;
			node_ss[i].du = dt * node_ss[i].v + point_five_minus_beta_dt_square * node_ss[i].a + beta_dt_square * g_a_vec[g_id_tmp];
			node_ss[i].u += node_ss[i].du;
			node_ss[i].v += one_minus_gamma_dt * node_ss[i].a + gamma_dt * g_a_vec[g_id_tmp];
			node_ss[i].a = g_a_vec[g_id_tmp];
		}
		for (size_t i = 0; i < node_f_num; i++)
		{
			g_id_tmp = node_s_num + i;
			node_fs[i].p += dt * node_fs[i].pv + point_five_minus_beta_dt_square * node_fs[i].pa + beta_dt_square * g_a_vec[g_id_tmp];
			node_fs[i].pv += one_minus_gamma_dt * node_fs[i].pa + gamma_dt * g_a_vec[g_id_tmp];
			node_fs[i].pa = g_a_vec[g_id_tmp];
		}

		// Update variables on Gauss points
		double de11;
		double de_vol_tmp;
		for (size_t i = 0; i < elem_num; i++)
		{
			cur_elem = elems + i;
			cur_nodes1 = node_ss + cur_elem->index_s_x;
			cur_nodes2 = node_ss + cur_elem->index_s_x + 1;
			cur_nodes3 = node_ss + cur_elem->index_s_x + 2;
			cur_gps1 = &(elems[i].gps1);
			cur_gps2 = &(elems[i].gps2);

			de11 = cur_nodes1->du * cur_gps1->dNs1_dx
				+ cur_nodes2->du * cur_gps1->dNs2_dx
				+ cur_nodes3->du * cur_gps1->dNs3_dx;
			cur_gps1->e11 += de11;
			cur_gps1->s11 += de11 * elems[i].E;
			de_vol_tmp = de11;
			cur_gps1->n = (de_vol_tmp + cur_gps1->n) / (1.0 + de_vol_tmp);
			cur_gps1->density_avg = (1.0 - cur_gps1->n) * cur_elem->density_s + cur_gps1->n * cur_elem->density_f;

			de11 = cur_nodes1->du * cur_gps2->dNs1_dx
				+ cur_nodes2->du * cur_gps2->dNs2_dx
				+ cur_nodes3->du * cur_gps2->dNs3_dx;
			cur_gps2->e11 += de11;
			cur_gps2->s11 += de11 * elems[i].E;
			de_vol_tmp = de11;
			cur_gps2->n = (de_vol_tmp + cur_gps2->n) / (1.0 + de_vol_tmp);
			cur_gps2->density_avg = (1.0 - cur_gps1->n) * cur_elem->density_s + cur_gps1->n * cur_elem->density_f;

			//std::cout << "n: " << cur_gp1->n << " de11: " << de11
			//	      << " s11: " << cur_gp1->s11 << std::endl;

			// need to calculate w in the future
		}

		return 0;
	}


	int output(std::ofstream &out)
	{
		// output current time
		out << cur_t << ", ";
		// output displacement
		for (size_t i = 0; i < node_s_num; i++)
		{
			out << node_ss[i].u << ", ";
		}
		for (size_t i = 0; i < node_f_num; i++)
		{
			out << node_fs[i].p << ", ";
		}
		out << std::endl;

		return 0;
	}

}