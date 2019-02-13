#include <cstddef>
#include <cstdlib>

#include <fstream>
#include <iostream>

#include <Eigen/Dense>

#include "hyme1D_MPM.h"

/*
Note:
1. Quadratic shape function for displacement;
2. Linear shape function for pore pressure;
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

namespace hyme1D_MPM
{
	size_t node_s_num;
	Node_s *node_ss;
	size_t node_f_num;
	Node_f *node_fs;

	size_t elem_num;
	Element *elems;

	size_t pcl_num;
	Particle *pcls;
	
	// We assume that the same body force is applied to both solid and fluid phase.
	size_t body_force_num;
	BodyForce *body_force;

	// Neumann boundary conditions
	size_t surface_force_num;
	SurfaceForce *surface_force;

	// Dirichlet boundary conditions (specifed on nodes)
	size_t acceleration_num;
	AccelerationBC *acceleration_bc;

	size_t p_acceleration_num;
	PAccelerationBC *p_acceleration_bc;

	//Initial boundary conditions (specifed on material points)
	size_t displacement_num;
	DisplacementBC *displacement_bc;

	size_t velocity_num;
	VelocityBC *velocity_bc;

	size_t p_velocity_num;
	PVelocityBC *p_velocity_bc;

	size_t pore_pressure_num;
	PorePressureBC *pore_pressure_bc;

	// Step time length
	double total_t, cur_t;
	// Substep time increment
	double dt;

	size_t *node_s_id_map, *node_f_id_map; // from node id to dof id
	size_t *dof_id_map; // from dof id to node id 
	size_t dof_s_num, dof_f_num, dof_num;
	size_t total_possible_dof_num; // used to identified whether a node is valid

	// Temporary calculation variables
	// v + 0.5 * a * dt
	double *tmp_coef3;
	// v * dt + 0.25 * a * dt * dt
	double *tmp_coef4;
	// constants
	double gamma_dt;
	double one_minus_gamma_dt;
	double beta_dt_square;
	double point_five_minus_beta_dt_square;

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

	// Needed to be optimized ...
	Element *findInWhichElement(double x)
	{
		size_t node_id_tmp;
		for (size_t i = 0; i < elem_num; i++)
		{
			node_id_tmp = elems[i].index_f_x;
			if (node_fs[node_id_tmp].x < x && x <= node_fs[node_id_tmp + 1].x)
			{
				return elems + i;
			}
		}
		return NULL;
	}

	int init(void);
	// Calculate each substeps
	int substep_init(void);
	int substep_map_to_node(void);
	int substep_solve(void);
	int substep_map_to_pcl(void);
	int output(std::ofstream &out);

	// calculate the whole step
	// simplify ver.
	int cal(std::ofstream &out);

	int hyme1D_MPM(void)
	{
		/* init mesh data */
		elem_num = 10;
		node_s_num = 2 * elem_num + 1;
		node_f_num = elem_num + 1;
		node_ss = (Node_s *) malloc(sizeof(Node_s)  * node_s_num);
		node_fs = (Node_f *) malloc(sizeof(Node_f)  * node_f_num);
		elems   = (Element *)malloc(sizeof(Element) * elem_num);
		for (size_t i = 0; i < node_s_num; i++)
		{
			node_ss[i].x = (double)i * 0.5;
		}
		for (size_t i = 0; i < node_f_num; i++)
		{
			node_fs[i].x = (double)i;
		}
		for (size_t i = 0; i < elem_num; i++)
		{
			elems[i].index_s_x = 2 * i;
			elems[i].index_f_x = i;
		}

		/* init particles data */
		pcl_num = elem_num * 5;
		pcls = (Particle *)malloc(sizeof(Particle) * pcl_num);
		for (size_t i = 0; i < pcl_num; i++)
		{
			pcls[i].x = 0.1 + i * 0.2;
			
			pcls[i].n = 0.5;
			pcls[i].density_s = 3000.0;
			pcls[i].density_f = 1000.0;
			pcls[i].vol = 0.2;

			pcls[i].E = 40.0e6;
			pcls[i].k = 1.0e-5;
			pcls[i].miu = 1.0;

			pcls[i].s11 = 0.0;
			pcls[i].p = 0.0;
			pcls[i].e11 = 0.0;

			// default initialization
			pcls[i].mass_s = pcls[i].vol * pcls[i].density_s * (1.0 - pcls[i].n);
			pcls[i].mass_f = pcls[i].vol * pcls[i].density_f * pcls[i].n;
			pcls[i].k_div_miu = pcls[i].k / pcls[i].miu;
			pcls[i].density_avg = (1.0 - pcls[i].n) * pcls[i].density_s + pcls[i].n  * pcls[i].density_f;
			pcls[i].a  = 0.0;
			pcls[i].v  = 0.0;
			pcls[i].pa = 0.0;
			pcls[i].pv = 0.0;
			pcls[i].p  = 0.0;
		}

		/* init BCs */
		// Body force
		body_force_num = pcl_num;
		body_force = (BodyForce *)malloc(sizeof(BodyForce) * body_force_num);
		for (size_t i = 0; i < body_force_num; i++)
		{
			body_force[i].bf_x = 0.0;
			body_force[i].pcl_id = i;
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
		// Acceleration of pore pressure BC
		//p_acceleration_num = 1;
		//p_acceleration_bc = (PAccelerationBC *)malloc(sizeof(PAccelerationBC) * p_acceleration_num);
		//p_acceleration_bc[0].pa = 0.0;
		//p_acceleration_bc[0].node_id = node_f_num - 1;
		// impermeable boundary
		p_acceleration_num = 0;
		p_acceleration_bc = NULL;

		// Initial bounday conditions at nodes
		velocity_num = 0;
		velocity_bc = NULL;
		displacement_num = 0;
		displacement_bc = NULL;
		p_velocity_num = 0;
		p_velocity_bc = NULL;
		pore_pressure_num = 0;
		pore_pressure_bc = NULL;

		std::ofstream out("res.txt");


		init();

		total_t = 1.0;
		dt = 0.01;
		cal(out);

		// pa
		p_acceleration_num = 1;
		p_acceleration_bc = (PAccelerationBC *)malloc(sizeof(PAccelerationBC) * p_acceleration_num);
		p_acceleration_bc[0].pa = 0.0;
		p_acceleration_bc[0].node_id = node_f_num - 1;
		// pv
		//p_velocity_num = 1;
		//p_velocity_bc = (PVelocityBC *)malloc(sizeof(PVelocityBC) * p_velocity_num);
		//p_velocity_bc[0].pv = 0.0;
		//p_velocity_bc[0].node_id = node_f_num - 1;
		// p
		pore_pressure_num = 1;
		pore_pressure_bc = (PorePressureBC *)malloc(sizeof(PorePressureBC) * pore_pressure_num);
		pore_pressure_bc[0].p = 0.0;
		pore_pressure_bc[0].node_id = node_f_num - 1;

		total_t = 3.0;
		dt = 0.01;
		cal(out);

		// Delete data
		if (node_ss) free(node_ss);
		if (node_fs) free(node_fs);
		if (elems) free(elems);
		if (pcls) free(pcls);
		if (body_force) free(body_force);
		// Dirichlet boundary conditions
		if (acceleration_bc) free(acceleration_bc);
		if (p_acceleration_bc) free(p_acceleration_bc);
		// Neumann boundary conditions
		if (surface_force) free(surface_force);
		// Initial boundary conditions
		if (velocity_bc) free(velocity_bc);
		if (displacement_bc) free(displacement_bc);
		if (p_velocity_bc) free(p_velocity_bc);
		if (pore_pressure_bc) free(pore_pressure_bc);

		// Node id mapping
		if (node_s_id_map) free(node_s_id_map);
		if (node_f_id_map) free(node_f_id_map);
		if (dof_id_map) free(dof_id_map);

		// Output file
		out.close();

		if (tmp_coef3) free(tmp_coef3);
		if (tmp_coef4) free(tmp_coef4);

		return 0;
	}


	// initialize calculation
	int init(void)
	{
		cur_t = 0.0;

		// id mapping between node and dof
		total_possible_dof_num = node_s_num + node_f_num;
		node_s_id_map = (size_t *)malloc(sizeof(size_t) * node_s_num);
		node_f_id_map = (size_t *)malloc(sizeof(size_t) * node_f_num);
		dof_id_map    = (size_t *)malloc(sizeof(size_t) * total_possible_dof_num);

		// v + (1 - gamma) * a * dt
		tmp_coef3 = (double *)malloc(sizeof(double) * node_s_num);
		// v * dt + (0.5 - beta) * a * dt * dt
		tmp_coef4 = (double *)malloc(sizeof(double) * (node_s_num + node_f_num));

		//Apply initial boundary conditions
		/*
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
		*/

		return 0;
	};

	int substep_init(void)
	{
		gamma_dt = gamma * dt;
		one_minus_gamma_dt = (1.0 - gamma) * dt;
		beta_dt_square = beta * dt * dt;
		point_five_minus_beta_dt_square = (0.5 - beta) * dt * dt;

		for (size_t i = 0; i < elem_num; i++) elems[i].clear();

		Particle *cur_pcl;
		for (size_t i = 0; i < pcl_num; i++)
		{
			cur_pcl = pcls + i;

			cur_pcl->inElem = findInWhichElement(cur_pcl->x);
			//std::cout << "elem index: " << cur_pcl->inElem->index_f_x << std::endl;
			if (cur_pcl->inElem)
			{
				// Add particles to element
				cur_pcl->inElem->add(cur_pcl);

				// Cal particle shape functions
				naturalCoord_s(cur_pcl->inElem, &(cur_pcl->x), &(cur_pcl->xi));
				/* solid phase */
				cur_pcl->Ns1 = _Ns1(cur_pcl->xi);
				cur_pcl->Ns2 = _Ns2(cur_pcl->xi);
				cur_pcl->Ns3 = _Ns3(cur_pcl->xi);
				cur_pcl->dx_dxi_s = dx_dxi_s(cur_pcl->inElem, cur_pcl->xi);
				cur_pcl->dxi_dx_s = 1.0 / cur_pcl->dx_dxi_s;
				cur_pcl->dNs1_dx = _dNs1_dxi(cur_pcl->xi) * cur_pcl->dxi_dx_s;
				cur_pcl->dNs2_dx = _dNs2_dxi(cur_pcl->xi) * cur_pcl->dxi_dx_s;
				cur_pcl->dNs3_dx = _dNs3_dxi(cur_pcl->xi) * cur_pcl->dxi_dx_s;
				/* fluid phase */
				cur_pcl->Nf1 = _Nf1(cur_pcl->xi);
				cur_pcl->Nf2 = _Nf2(cur_pcl->xi);
				cur_pcl->dx_dxi_f = dx_dxi_f(cur_pcl->inElem, cur_pcl->xi);
				cur_pcl->dxi_dx_f = 1.0 / cur_pcl->dx_dxi_f;
				cur_pcl->dNf1_dx = _dNf1_dxi(cur_pcl->xi) * cur_pcl->dxi_dx_f;
				cur_pcl->dNf2_dx = _dNf2_dxi(cur_pcl->xi) * cur_pcl->dxi_dx_f;
			}
		}

		//elems[9].print_pcl();

		total_possible_dof_num = node_s_num + node_f_num; // To identify invalid node
		for (size_t i = 0; i < node_s_num; i++)
		{
			node_s_id_map[i] = total_possible_dof_num;
		}
		for (size_t i = 0; i < node_f_num; i++)
		{
			node_f_id_map[i] = total_possible_dof_num;
		}

		// Form node id mapping
		size_t dof_s_id, dof_f_id;
		size_t node_id_tmp;
		dof_s_id = 0;
		dof_f_id = 0;
		// Count valid nodes and assign id
		// Valid nodes refer to nodes of elements that contains particles. 
		for (size_t i = 0; i < elem_num; i++)
		{
			if (elems[i].top) // if this element contains particles
			{
				// nodes for solid phase
				node_id_tmp = elems[i].index_s_x;
				if (node_s_id_map[node_id_tmp] == total_possible_dof_num)
				{
					node_s_id_map[node_id_tmp] = dof_s_id;
					++dof_s_id;
					// init node variables
					node_ss[node_id_tmp].m = 0.0;
					node_ss[node_id_tmp].ma = 0.0;
					node_ss[node_id_tmp].mv = 0.0;
					node_ss[node_id_tmp].a = 0.0;
					node_ss[node_id_tmp].v = 0.0;
					node_ss[node_id_tmp].dv = 0.0;
					node_ss[node_id_tmp].du = 0.0;
				}
				++node_id_tmp;
				if (node_s_id_map[node_id_tmp] == total_possible_dof_num)
				{
					node_s_id_map[node_id_tmp] = dof_s_id;
					++dof_s_id;
					// init node variables
					node_ss[node_id_tmp].m = 0.0;
					node_ss[node_id_tmp].ma = 0.0;
					node_ss[node_id_tmp].mv = 0.0;
					node_ss[node_id_tmp].a = 0.0;
					node_ss[node_id_tmp].v = 0.0;
					node_ss[node_id_tmp].dv = 0.0;
					node_ss[node_id_tmp].du = 0.0;
				}
				++node_id_tmp;
				if (node_s_id_map[node_id_tmp] == total_possible_dof_num)
				{
					node_s_id_map[node_id_tmp] = dof_s_id;
					++dof_s_id;
					// init node variables
					node_ss[node_id_tmp].m = 0.0;
					node_ss[node_id_tmp].ma = 0.0;
					node_ss[node_id_tmp].mv = 0.0;
					node_ss[node_id_tmp].a = 0.0;
					node_ss[node_id_tmp].v = 0.0;
					node_ss[node_id_tmp].dv = 0.0;
					node_ss[node_id_tmp].du = 0.0;
				}
				// nodes for fluid phase
				node_id_tmp = elems[i].index_f_x;
				if (node_f_id_map[node_id_tmp] == total_possible_dof_num)
				{
					node_f_id_map[node_id_tmp] = dof_f_id;
					++dof_f_id;
					// init node variables
					node_fs[node_id_tmp].m = 0.0;
					node_fs[node_id_tmp].mpa = 0.0;
					node_fs[node_id_tmp].mpv = 0.0;
					node_fs[node_id_tmp].mp = 0.0;
					node_fs[node_id_tmp].pa = 0.0;
					node_fs[node_id_tmp].pv = 0.0;
					node_fs[node_id_tmp].p = 0.0;
					node_fs[node_id_tmp].dpv = 0.0;
					node_fs[node_id_tmp].dp = 0.0;
				}
				++node_id_tmp;
				if (node_f_id_map[node_id_tmp] == total_possible_dof_num)
				{
					node_f_id_map[node_id_tmp] = dof_f_id;
					++dof_f_id;
					// init node variables
					node_fs[node_id_tmp].m = 0.0;
					node_fs[node_id_tmp].mpa = 0.0;
					node_fs[node_id_tmp].mpv = 0.0;
					node_fs[node_id_tmp].mp = 0.0;
					node_fs[node_id_tmp].pa = 0.0;
					node_fs[node_id_tmp].pv = 0.0;
					node_fs[node_id_tmp].p = 0.0;
					node_fs[node_id_tmp].dpv = 0.0;
					node_fs[node_id_tmp].dp = 0.0;
				}
			}
		}
		dof_s_num = dof_s_id;
		dof_f_num = dof_f_id;
		dof_num = dof_s_num + dof_f_num;
		// form dof to node id map
		for (size_t i = 0; i < node_s_num; i++)
		{
			if (node_s_id_map[i] != total_possible_dof_num)
			{
				dof_id_map[node_s_id_map[i]] = i;
			}
		}
		for (size_t i = 0; i < node_f_num; i++)
		{
			if (node_f_id_map[i] != total_possible_dof_num)
			{
				node_f_id_map[i] += dof_s_num;
				dof_id_map[node_f_id_map[i]] = i;
			}
		}

		/*
		for (size_t i = 0; i < node_s_num; i++)
		{
			std::cout << node_s_id_map[i] << " ";
		}
		std::cout << std::endl;
		for (size_t i = 0; i < node_f_num; i++)
		{
			std::cout << node_f_id_map[i] << " ";
		}
		std::cout << std::endl;
		for (size_t i = 0; i < dof_num; i++)
		{
			std::cout << dof_id_map[i] << " ";
		}
		std::cout << std::endl;
		*/

		return 0;
	}

	int substep_map_to_node(void)
	{
		Particle *cur_pcl;
		Element *cur_elem;
		Node_s *node_s1, *node_s2, *node_s3;
		Node_f *node_f1, *node_f2;

		for (size_t pcl_id = 0; pcl_id < pcl_num; pcl_id++)
		{
			cur_pcl = pcls + pcl_id;
			cur_elem = cur_pcl->inElem;
			node_s1 = node_ss + cur_elem->index_s_x;
			node_s2 = node_ss + cur_elem->index_s_x + 1;
			node_s3 = node_ss + cur_elem->index_s_x + 2;
			node_f1 = node_fs + cur_elem->index_f_x;
			node_f2 = node_fs + cur_elem->index_f_x + 1;

			// map acceleration and velocity
			node_s1->m  += cur_pcl->mass_s * cur_pcl->Ns1;
			node_s1->ma += cur_pcl->a * cur_pcl->mass_s * cur_pcl->Ns1;
			node_s1->mv += cur_pcl->v * cur_pcl->mass_s * cur_pcl->Ns1;
			node_s2->m  += cur_pcl->mass_s * cur_pcl->Ns2;
			node_s2->ma += cur_pcl->a * cur_pcl->mass_s * cur_pcl->Ns2;
			node_s2->mv += cur_pcl->v * cur_pcl->mass_s * cur_pcl->Ns2;
			node_s3->m  += cur_pcl->mass_s * cur_pcl->Ns3;
			node_s3->ma += cur_pcl->a * cur_pcl->mass_s * cur_pcl->Ns3;
			node_s3->mv += cur_pcl->v * cur_pcl->mass_s * cur_pcl->Ns3;

			// map pore pressure "acceleration" and "velocity"
			node_f1->m   += cur_pcl->mass_f * cur_pcl->Nf1;
			node_f1->mpa += cur_pcl->pa * cur_pcl->mass_f * cur_pcl->Nf1;
			node_f1->mpv += cur_pcl->pv * cur_pcl->mass_f * cur_pcl->Nf1;
			node_f1->mp  += cur_pcl->p  * cur_pcl->mass_f * cur_pcl->Nf1;
			node_f2->m   += cur_pcl->mass_f * cur_pcl->Nf2;
			node_f2->mpa += cur_pcl->pa * cur_pcl->mass_f * cur_pcl->Nf2;
			node_f2->mpv += cur_pcl->pv * cur_pcl->mass_f * cur_pcl->Nf2;
			node_f2->mp  += cur_pcl->p  * cur_pcl->mass_f * cur_pcl->Nf2;
		}


		for (size_t node_id = 0; node_id < node_s_num; node_id++)
		{
			if (node_s_id_map[node_id] != total_possible_dof_num) // if this node is valid
			{
				node_ss[node_id].a = node_ss[node_id].ma / node_ss[node_id].m;
				node_ss[node_id].v = node_ss[node_id].mv / node_ss[node_id].m;
			}
		}

		/*
		for (size_t i = 0; i < node_s_num; i++)
			std::cout << node_ss[i].a << " ";
		std::cout << std::endl;
		for (size_t i = 0; i < node_s_num; i++)
			std::cout << node_ss[i].v << " ";
		std::cout << std::endl;
		*/

		for (size_t node_id = 0; node_id < node_f_num; node_id++)
		{
			if (node_f_id_map[node_id] != total_possible_dof_num) // if this node is valid
			{
				node_fs[node_id].pa = node_fs[node_id].mpa / node_fs[node_id].m;
				node_fs[node_id].pv = node_fs[node_id].mpv / node_fs[node_id].m;
				node_fs[node_id].p  = node_fs[node_id].mp  / node_fs[node_id].m;
			}
		}

		/*
		for (size_t i = 0; i < node_f_num; i++)
			std::cout << node_fs[i].pa << " ";
		std::cout << std::endl;
		for (size_t i = 0; i < node_f_num; i++)
			std::cout << node_fs[i].pv << " ";
		std::cout << std::endl;
		for (size_t i = 0; i < node_f_num; i++)
			std::cout << node_fs[i].p << " ";
		std::cout << std::endl;
		*/

		// Apply acceleration boundary condition
		for (size_t i = 0; i < acceleration_num; i++)
		{
			node_ss[acceleration_bc[i].node_id].a = acceleration_bc[i].a;
		}
		for (size_t i = 0; i < p_acceleration_num; i++)
		{
			node_fs[p_acceleration_bc[i].node_id].pa = p_acceleration_bc[i].pa;
		}

		// Modify initial conditions at nodes
		// Velocity of pore pressure
		for (size_t i = 0; i < p_velocity_num; i++)
		{
			node_fs[p_velocity_bc[i].node_id].p = p_velocity_bc[i].pv;
			// simple method
			node_fs[p_velocity_bc[i].node_id].dpv += p_velocity_bc[i].pv - node_fs[p_velocity_bc[i].node_id].pv;
			// more complex method (may be better)
		}
		// Pore pressure
		for (size_t i = 0; i < pore_pressure_num; i++)
		{
			node_fs[pore_pressure_bc[i].node_id].p = pore_pressure_bc[i].p;
			// simple method
			node_fs[pore_pressure_bc[i].node_id].dp += pore_pressure_bc[i].p - node_fs[pore_pressure_bc[i].node_id].p;
			// more complex method (may be better)
		}

		return 0;
	}

	int substep_solve(void)
	{
		// Global stiffness matrix
		Eigen::SparseMatrix<double> g_st_mat(dof_num, dof_num);
		std::list<Eigen::Triplet<double> > g_st_mat_coef;
		// Global force vector
		Eigen::VectorXd g_f_vec(dof_num);
		// Global acceleration vector
		Eigen::VectorXd g_a_vec(dof_num);

		// Elemental stiffness matrix
		double elem_st_mat[5][5];
		double elem_m_mat[5][5], elem_c_mat[5][5], elem_k_mat[5][5];

		// Elemental force vector
		double elem_f_vec[5];
		Element *cur_elem;
		Node_s *cur_nodes1, *cur_nodes2, *cur_nodes3;
		Node_f *cur_nodef1, *cur_nodef2;
		Particle *cur_pcl;
		/* Index mapping array
		used to map coefficient of elemental stiffness
		matrix into glabal stiffness matrix.
		*/
		size_t cur_g_id[5];

		// Temporary calculation param
		size_t dof_to_node_tmp;
		for (size_t i = 0; i < dof_s_num; i++)
		{
			dof_to_node_tmp = dof_id_map[i];
			// v + (1 - gamma) * dt * a
			tmp_coef3[i] = node_ss[dof_to_node_tmp].v + one_minus_gamma_dt * node_ss[dof_to_node_tmp].a;
			// v * dt + (0.5 - beta) * dt * dt * a
			tmp_coef4[i] = dt * node_ss[dof_to_node_tmp].v + point_five_minus_beta_dt_square * node_ss[dof_to_node_tmp].a;
		}
		for (size_t i = 0; i < dof_f_num; i++)
		{
			dof_to_node_tmp = dof_id_map[dof_s_num + i];
			tmp_coef4[dof_s_num + i] = dt * node_fs[dof_to_node_tmp].pv + point_five_minus_beta_dt_square * node_fs[dof_to_node_tmp].pa;
		}

		g_f_vec.setZero();
		
		/* form element stiffness matrix and force vector */
		for (size_t elem_id = 0; elem_id < elem_num; elem_id++)
		{
			cur_elem = elems + elem_id;
			if (!(cur_elem->top)) continue; // if this element has no particles

			cur_nodes1 = node_ss + cur_elem->index_s_x;
			cur_nodes2 = node_ss + cur_elem->index_s_x + 1;
			cur_nodes3 = node_ss + cur_elem->index_s_x + 2;
			cur_nodef1 = node_fs + cur_elem->index_f_x;
			cur_nodef2 = node_fs + cur_elem->index_f_x + 1;

			// init index mapping array
			cur_g_id[0] = node_s_id_map[cur_elem->index_s_x];
			cur_g_id[1] = node_s_id_map[cur_elem->index_s_x + 1];
			cur_g_id[2] = node_s_id_map[cur_elem->index_s_x + 2];
			cur_g_id[3] = node_f_id_map[cur_elem->index_f_x];
			cur_g_id[4] = node_f_id_map[cur_elem->index_f_x + 1];

			// Init elemental stiffness matrix
			memset(elem_st_mat, 0, sizeof(double) * 5 * 5);
			memset(elem_m_mat,  0, sizeof(double) * 5 * 5);
			memset(elem_c_mat,  0, sizeof(double) * 5 * 5);
			memset(elem_k_mat,  0, sizeof(double) * 5 * 5);

#define ElementFirstParticle(elem) ((elem).top)
#define ElementNextParticle(pcl) ((pcl).next)
			for (cur_pcl = ElementFirstParticle(*cur_elem); cur_pcl; cur_pcl = ElementNextParticle(*cur_pcl))
			{
				// M
				elem_m_mat[0][0] += cur_pcl->density_avg * cur_pcl->Ns1 * cur_pcl->Ns1 * cur_pcl->vol;
				elem_m_mat[0][1] += cur_pcl->density_avg * cur_pcl->Ns1 * cur_pcl->Ns2 * cur_pcl->vol;
				elem_m_mat[0][2] += cur_pcl->density_avg * cur_pcl->Ns1 * cur_pcl->Ns3 * cur_pcl->vol;
				elem_m_mat[1][0] += cur_pcl->density_avg * cur_pcl->Ns2 * cur_pcl->Ns1 * cur_pcl->vol;
				elem_m_mat[1][1] += cur_pcl->density_avg * cur_pcl->Ns2 * cur_pcl->Ns2 * cur_pcl->vol;
				elem_m_mat[1][2] += cur_pcl->density_avg * cur_pcl->Ns2 * cur_pcl->Ns3 * cur_pcl->vol;
				elem_m_mat[2][0] += cur_pcl->density_avg * cur_pcl->Ns3 * cur_pcl->Ns1 * cur_pcl->vol;
				elem_m_mat[2][1] += cur_pcl->density_avg * cur_pcl->Ns3 * cur_pcl->Ns2 * cur_pcl->vol;
				elem_m_mat[2][2] += cur_pcl->density_avg * cur_pcl->Ns3 * cur_pcl->Ns3 * cur_pcl->vol;
				/*
				// may be ignored?
				elem_m_mat[3][0] = cur_pcl->density_f * cur_pcl->k_div_miu * cur_pcl->dNf1_dx * cur_pcl->Ns1 * cur_pcl->w;
				elem_m_mat[3][1] = cur_pcl->density_f * cur_pcl->k_div_miu * cur_pcl->dNf1_dx * cur_pcl->Ns2 * cur_pcl->w;
				elem_m_mat[3][2] = cur_pcl->density_f * cur_pcl->k_div_miu * cur_pcl->dNf1_dx * cur_pcl->Ns3 * cur_pcl->w;
				elem_m_mat[4][0] = cur_pcl->density_f * cur_pcl->k_div_miu * cur_pcl->dNf2_dx * cur_pcl->Ns1 * cur_pcl->w;
				elem_m_mat[4][1] = cur_pcl->density_f * cur_pcl->k_div_miu * cur_pcl->dNf2_dx * cur_pcl->Ns2 * cur_pcl->w;
				elem_m_mat[4][2] = cur_pcl->density_f * cur_pcl->k_div_miu * cur_pcl->dNf2_dx * cur_pcl->Ns3 * cur_pcl->w;
				*/
				// C
				elem_c_mat[3][0] += cur_pcl->Nf1 * cur_pcl->dNs1_dx * cur_pcl->vol;
				elem_c_mat[3][1] += cur_pcl->Nf1 * cur_pcl->dNs2_dx * cur_pcl->vol;
				elem_c_mat[3][2] += cur_pcl->Nf1 * cur_pcl->dNs3_dx * cur_pcl->vol;
				elem_c_mat[4][0] += cur_pcl->Nf2 * cur_pcl->dNs1_dx * cur_pcl->vol;
				elem_c_mat[4][1] += cur_pcl->Nf2 * cur_pcl->dNs2_dx * cur_pcl->vol;
				elem_c_mat[4][2] += cur_pcl->Nf2 * cur_pcl->dNs3_dx * cur_pcl->vol;
				// K
				elem_k_mat[0][0] += cur_pcl->dNs1_dx * cur_pcl->E * cur_pcl->dNs1_dx * cur_pcl->vol;
				elem_k_mat[0][1] += cur_pcl->dNs1_dx * cur_pcl->E * cur_pcl->dNs2_dx * cur_pcl->vol;
				elem_k_mat[0][2] += cur_pcl->dNs1_dx * cur_pcl->E * cur_pcl->dNs3_dx * cur_pcl->vol;
				elem_k_mat[0][3] += -cur_pcl->dNs1_dx * cur_pcl->Nf1 * cur_pcl->vol;
				elem_k_mat[0][4] += -cur_pcl->dNs1_dx * cur_pcl->Nf2 * cur_pcl->vol;
				elem_k_mat[1][0] += cur_pcl->dNs2_dx * cur_pcl->E * cur_pcl->dNs1_dx * cur_pcl->vol;
				elem_k_mat[1][1] += cur_pcl->dNs2_dx * cur_pcl->E * cur_pcl->dNs2_dx * cur_pcl->vol;
				elem_k_mat[1][2] += cur_pcl->dNs2_dx * cur_pcl->E * cur_pcl->dNs3_dx * cur_pcl->vol;
				elem_k_mat[1][3] += -cur_pcl->dNs2_dx * cur_pcl->Nf1 * cur_pcl->vol;
				elem_k_mat[1][4] += -cur_pcl->dNs2_dx * cur_pcl->Nf2 * cur_pcl->vol;
				elem_k_mat[2][0] += cur_pcl->dNs3_dx * cur_pcl->E * cur_pcl->dNs1_dx * cur_pcl->vol;
				elem_k_mat[2][1] += cur_pcl->dNs3_dx * cur_pcl->E * cur_pcl->dNs2_dx * cur_pcl->vol;
				elem_k_mat[2][2] += cur_pcl->dNs3_dx * cur_pcl->E * cur_pcl->dNs3_dx * cur_pcl->vol;
				elem_k_mat[2][3] += -cur_pcl->dNs3_dx * cur_pcl->Nf1 * cur_pcl->vol;
				elem_k_mat[2][4] += -cur_pcl->dNs3_dx * cur_pcl->Nf2 * cur_pcl->vol;
				elem_k_mat[3][3] += cur_pcl->k_div_miu * cur_pcl->dNf1_dx * cur_pcl->dNf1_dx * cur_pcl->vol;
				elem_k_mat[3][4] += cur_pcl->k_div_miu * cur_pcl->dNf1_dx * cur_pcl->dNf2_dx * cur_pcl->vol;
				elem_k_mat[4][3] += cur_pcl->k_div_miu * cur_pcl->dNf2_dx * cur_pcl->dNf1_dx * cur_pcl->vol;
				elem_k_mat[4][4] += cur_pcl->k_div_miu * cur_pcl->dNf2_dx * cur_pcl->dNf2_dx * cur_pcl->vol;
			}
			// M + gamma * dt * C + beta * dt * dt * K
#define FormStiffnessMatrix(i, j) elem_st_mat[i][j] = elem_m_mat[i][j] + gamma_dt * elem_c_mat[i][j] + beta_dt_square * elem_k_mat[i][j]
			FormStiffnessMatrix(0, 0);
			FormStiffnessMatrix(0, 1);
			FormStiffnessMatrix(0, 2);
			FormStiffnessMatrix(0, 3);
			FormStiffnessMatrix(0, 4);
			FormStiffnessMatrix(1, 0);
			FormStiffnessMatrix(1, 1);
			FormStiffnessMatrix(1, 2);
			FormStiffnessMatrix(1, 3);
			FormStiffnessMatrix(1, 4);
			FormStiffnessMatrix(2, 0);
			FormStiffnessMatrix(2, 1);
			FormStiffnessMatrix(2, 2);
			FormStiffnessMatrix(2, 3);
			FormStiffnessMatrix(2, 4);
			FormStiffnessMatrix(3, 0);
			FormStiffnessMatrix(3, 1);
			FormStiffnessMatrix(3, 2);
			FormStiffnessMatrix(3, 3);
			FormStiffnessMatrix(3, 4);
			FormStiffnessMatrix(4, 0);
			FormStiffnessMatrix(4, 1);
			FormStiffnessMatrix(4, 2);
			FormStiffnessMatrix(4, 3);
			FormStiffnessMatrix(4, 4);

			// Init elemental force vector
			memset(elem_f_vec, 0, sizeof(double) * 5);
			for (cur_pcl = ElementFirstParticle(*cur_elem); cur_pcl; cur_pcl = ElementNextParticle(*cur_pcl))
			{
				// - C(v+ 0.5 * dt * a)
				double cur_node1_v = tmp_coef3[cur_g_id[0]];
				double cur_node2_v = tmp_coef3[cur_g_id[1]];
				double cur_node3_v = tmp_coef3[cur_g_id[2]];
				double cur_pcl_v = (cur_pcl->dNs1_dx * cur_node1_v   // node1
								  + cur_pcl->dNs2_dx * cur_node2_v   // node2
								  + cur_pcl->dNs3_dx * cur_node3_v); // node3
				elem_f_vec[3] -= cur_pcl->Nf1 * cur_pcl_v * cur_pcl->vol;
				elem_f_vec[4] -= cur_pcl->Nf2 * cur_pcl_v * cur_pcl->vol;
				// -N (dn)
				elem_f_vec[0] -= cur_pcl->dNs1_dx * (cur_pcl->s11 - cur_pcl->Nf1 * cur_nodef1->p - cur_pcl->Nf2 * cur_nodef2->p) * cur_pcl->vol;
				elem_f_vec[1] -= cur_pcl->dNs2_dx * (cur_pcl->s11 - cur_pcl->Nf1 * cur_nodef1->p - cur_pcl->Nf2 * cur_nodef2->p) * cur_pcl->vol;
				elem_f_vec[2] -= cur_pcl->dNs3_dx * (cur_pcl->s11 - cur_pcl->Nf1 * cur_nodef1->p - cur_pcl->Nf2 * cur_nodef2->p) * cur_pcl->vol;
				elem_f_vec[3] -= cur_pcl->k_div_miu * cur_pcl->dNf1_dx * (cur_pcl->dNf1_dx * cur_nodef1->p + cur_pcl->dNf2_dx * cur_nodef2->p) * cur_pcl->vol;
				elem_f_vec[4] -= cur_pcl->k_div_miu * cur_pcl->dNf2_dx * (cur_pcl->dNf1_dx * cur_nodef1->p + cur_pcl->dNf2_dx * cur_nodef2->p) * cur_pcl->vol;
				//std::cout << "force vector matrix: " << std::endl << elem_f_vec << std::endl;
				// - dN(v * dt + 0.25 * dt * dt * a)
				double cur_pcl_ds11 = cur_pcl->E *
					     (cur_pcl->dNs1_dx * tmp_coef4[cur_g_id[0]]
						+ cur_pcl->dNs2_dx * tmp_coef4[cur_g_id[1]]
						+ cur_pcl->dNs3_dx * tmp_coef4[cur_g_id[2]]);
				double cur_node1_dp = tmp_coef4[cur_g_id[3]];
				double cur_node2_dp = tmp_coef4[cur_g_id[4]];
				elem_f_vec[0] -= cur_pcl->dNs1_dx * (cur_pcl_ds11 - cur_pcl->Nf1 * cur_node1_dp - cur_pcl->Nf2 * cur_node2_dp) * cur_pcl->vol;
				elem_f_vec[1] -= cur_pcl->dNs2_dx * (cur_pcl_ds11 - cur_pcl->Nf1 * cur_node1_dp - cur_pcl->Nf2 * cur_node2_dp) * cur_pcl->vol;
				elem_f_vec[2] -= cur_pcl->dNs3_dx * (cur_pcl_ds11 - cur_pcl->Nf1 * cur_node1_dp - cur_pcl->Nf2 * cur_node2_dp) * cur_pcl->vol;
				elem_f_vec[3] -= cur_pcl->k_div_miu * cur_pcl->dNf1_dx * (cur_pcl->dNf1_dx * cur_node1_dp + cur_pcl->dNf2_dx * cur_node2_dp) * cur_pcl->vol;
				elem_f_vec[4] -= cur_pcl->k_div_miu * cur_pcl->dNf2_dx * (cur_pcl->dNf1_dx * cur_node1_dp + cur_pcl->dNf2_dx * cur_node2_dp) * cur_pcl->vol;
			}

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

		// Apply Neumann BC (directly modify the total force vector)
		double body_force_tmp;
		for (size_t i = 0; i < body_force_num; i++)
		{
			body_force_tmp = body_force[i].bf_x;
			cur_pcl = pcls + body_force[i].pcl_id;
			cur_elem = cur_pcl->inElem;

			// Init index mapping array
			cur_g_id[0] = node_s_id_map[cur_elem->index_s_x];
			cur_g_id[1] = node_s_id_map[cur_elem->index_s_x + 1];
			cur_g_id[2] = node_s_id_map[cur_elem->index_s_x + 2];
			cur_g_id[3] = node_f_id_map[cur_elem->index_f_x];
			cur_g_id[4] = node_f_id_map[cur_elem->index_f_x + 1];

			// Form external force vector for this element
			memset(elem_f_vec, 0, sizeof(double) * 5);
			elem_f_vec[0] = cur_pcl->density_avg * cur_pcl->Ns1 * body_force_tmp * cur_pcl->vol;
			elem_f_vec[1] = cur_pcl->density_avg * cur_pcl->Ns2 * body_force_tmp * cur_pcl->vol;
			elem_f_vec[2] = cur_pcl->density_avg * cur_pcl->Ns3 * body_force_tmp * cur_pcl->vol;
			elem_f_vec[3] = cur_pcl->density_f * cur_pcl->k_div_miu * cur_pcl->dNf1_dx * body_force_tmp * cur_pcl->vol;
			elem_f_vec[4] = cur_pcl->density_f * cur_pcl->k_div_miu * cur_pcl->dNf2_dx * body_force_tmp * cur_pcl->vol;

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
			if (!(cur_elem->top)) continue;

			// Init index mapping array
			cur_g_id[0] = node_s_id_map[cur_elem->index_s_x];
			cur_g_id[1] = node_s_id_map[cur_elem->index_s_x + 1];
			cur_g_id[2] = node_s_id_map[cur_elem->index_s_x + 2];
			cur_g_id[3] = node_f_id_map[cur_elem->index_f_x];
			cur_g_id[4] = node_f_id_map[cur_elem->index_f_x + 1];

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
		}
		// currently, the s * w type BC is ignored.
		//std::cout << "global force vector matrix: " << std::endl << g_f_vec << std::endl;

		/** Apply Dirichlet BC into force vector and stiffness matrix. **/
		// (The performance of this part can be enhanced.)
		// Acceleration
		size_t dof_id_tmp;
		double Knn_tmp;
		for (size_t i = 0; i < acceleration_num; i++)
		{
			dof_id_tmp = node_s_id_map[acceleration_bc[i].node_id];
			// modify stiffness matrix
			for (std::list<Eigen::Triplet<double> >::iterator iter = g_st_mat_coef.begin();
				iter != g_st_mat_coef.end();)
			{
				if (dof_id_tmp == iter->row() || dof_id_tmp == iter->col())
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
			g_f_vec[dof_id_tmp] = Knn_tmp * acceleration_bc[i].a;
		}
		// Pore pressure
		for (size_t i = 0; i < p_acceleration_num; i++)
		{
			dof_id_tmp = node_f_id_map[p_acceleration_bc[i].node_id];
			// modify stiffness matrix
			for (std::list<Eigen::Triplet<double> >::iterator iter = g_st_mat_coef.begin();
				iter != g_st_mat_coef.end();)
			{
				if (dof_id_tmp == iter->row() || dof_id_tmp == iter->col())
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
			g_f_vec[dof_id_tmp] = Knn_tmp * p_acceleration_bc[i].pa;
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

		// update node variables
		// scheme 1 ("more conventional")
		size_t node_id_tmp;
		for (size_t i = 0; i < dof_s_num; i++)
		{
			node_id_tmp = dof_id_map[i];
			node_ss[node_id_tmp].du += node_ss[node_id_tmp].v * dt + point_five_minus_beta_dt_square * node_ss[node_id_tmp].a + beta_dt_square * g_a_vec[i];
			node_ss[node_id_tmp].dv += one_minus_gamma_dt * node_ss[node_id_tmp].a + gamma_dt * g_a_vec[i];
			node_ss[node_id_tmp].a = g_a_vec[i];
		}
		for (size_t i = 0; i < dof_f_num; i++)
		{
			dof_id_tmp = dof_s_num + i;
			node_id_tmp = dof_id_map[dof_id_tmp];
			node_fs[node_id_tmp].dp  += node_fs[node_id_tmp].pv * dt
				+ point_five_minus_beta_dt_square * node_fs[node_id_tmp].pa + beta_dt_square * g_a_vec[dof_id_tmp];
			node_fs[node_id_tmp].dpv += one_minus_gamma_dt * node_fs[node_id_tmp].pa + gamma_dt * g_a_vec[dof_id_tmp];
			node_fs[node_id_tmp].pa = g_a_vec[dof_id_tmp];
		}

		// scheme 2
		// note: this scheme only need acceleration in nodal variables
		/*
		size_t dof_id_tmp;
		for (size_t i = 0; i < dof_s_num; i++)
		{
			dof_id_tmp = dof_id_map[i];
			node_ss[dof_id_tmp].a = g_a_vec[i];
		}
		for (size_t i = 0; i < dof_f_num; i++)
		{
			dof_id_tmp = dof_id_map[dof_s_num + i];
			node_fs[dof_id_tmp].pa = g_a_vec[dof_s_num + i];
		}
		*/

		return 0;
	}

	// Map variables back to particles
	int substep_map_to_pcl(void)
	{
		Particle *cur_pcl;
		Element *cur_elem;
		Node_s *node_s1, *node_s2, *node_s3;
		Node_f *node_f1, *node_f2;

		for (size_t pcl_id = 0; pcl_id < pcl_num; pcl_id++)
		{
			cur_pcl = pcls + pcl_id;
			cur_elem = cur_pcl->inElem;
			if (!cur_elem) continue;
			
			node_s1 = node_ss + cur_elem->index_s_x;
			node_s2 = node_ss + cur_elem->index_s_x + 1;
			node_s3 = node_ss + cur_elem->index_s_x + 2;
			node_f1 = node_fs + cur_elem->index_f_x;
			node_f2 = node_fs + cur_elem->index_f_x + 1;

			// scheme 1 (more "conventional")
			cur_pcl->a   = cur_pcl->Ns1 * node_s1->a  + cur_pcl->Ns2 * node_s2->a  + cur_pcl->Ns3 * node_s3->a;
			cur_pcl->v  += cur_pcl->Ns1 * node_s1->dv + cur_pcl->Ns2 * node_s2->dv + cur_pcl->Ns3 * node_s3->dv;
			cur_pcl->x  += cur_pcl->Ns1 * node_s1->du + cur_pcl->Ns2 * node_s2->du + cur_pcl->Ns3 * node_s3->du;
			cur_pcl->pa  = cur_pcl->Nf1 * node_f1->pa  + cur_pcl->Nf2 * node_f2->pa;
			cur_pcl->pv += cur_pcl->Nf1 * node_f1->dpv + cur_pcl->Nf2 * node_f2->dpv;
			cur_pcl->p  += cur_pcl->Nf1 * node_f1->dp  + cur_pcl->Nf2 * node_f2->dp;

			// scheme 2
			//a_tmp = cur_pcl->Ns1 * node_s1->a + cur_pcl->Ns2 * node_s2->a + cur_pcl->Ns3 * node_s3->a;
			//cur_pcl->du = cur_pcl->v * dt + point_five_minus_beta_dt_square * cur_pcl->a + beta_dt_square * a_tmp;
			//cur_pcl->v += one_minus_gamma_dt * cur_pcl->a + gamma_dt * a_tmp;
			//cur_pcl->a = a_tmp;
			//cur_pcl->x += cur_pcl->du;

			
			// update variables
			double de_vol_tmp; // volumetric strain
			double de11 = cur_pcl->dNs1_dx * node_s1->du + cur_pcl->dNs2_dx * node_s2->du + cur_pcl->dNs3_dx * node_s3->du;
			cur_pcl->e11 += de11;
			cur_pcl->s11 += cur_pcl->E * de11;
			de_vol_tmp = de11;
			cur_pcl->n = (de_vol_tmp + cur_pcl->n) / (1.0 + de_vol_tmp);
			cur_pcl->density_avg = (1.0 - cur_pcl->n) * cur_pcl->density_s + cur_pcl->n * cur_pcl->density_f;
		}
	
		return 0;
	}

	int cal(std::ofstream &out)
	{
		while (cur_t < total_t)
		{
			substep_init();
			substep_map_to_node();
			substep_solve();
			substep_map_to_pcl();
			cur_t += dt;
			output(out);
		}

		return 0;
	}

	int output(std::ofstream &out)
	{
		// output current time
		out << cur_t << ", ";
		// output displacement
		for (size_t i = 0; i < pcl_num; i++)
		{
			out << pcls[i].x << ", ";
		}
		for (size_t i = 0; i < pcl_num; i++)
		{
			out << pcls[i].p << ", ";
		}
		out << std::endl;

		return 0;
	}

}