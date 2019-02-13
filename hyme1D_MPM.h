#ifndef _HYME1D_MPM_H_
#define _HYME1D_MPM_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

namespace hyme1D_MPM
{
	struct Element;

	struct Particle
	{
		double x;

		double density_s, density_f;
		double vol; // the volume this particle represent 

		double a, v;
		double pa, pv, p;
		double s11;
		double e11;
		double n;

		double E;
		// permeability
		double k;
		// dynamic viscosity
		double miu;
		// bulk modulus
		//double K_f;

		Element *inElem;
		double xi; // natural coordinates

		// for solid phase
		double Ns1, Ns2, Ns3;
		double dNs1_dx, dNs2_dx, dNs3_dx;
		double dx_dxi_s;
		double dxi_dx_s;

		// for fluid phase
		double Nf1, Nf2;
		double dNf1_dx, dNf2_dx;
		double dx_dxi_f;
		double dxi_dx_f;

		// Temporary variables
		double mass_s, mass_f;
		double density_avg;
		double k_div_miu;

		// for particle stack in each elements 
		Particle *next;
	};

	struct Element
	{
		// topology
		size_t index_s_x; // node index of solid phase
		size_t index_f_x; // node index of fluid phase

						  // particles in this element
		Particle *top;

		Element() : top(nullptr) {}
		inline void add(Particle *pcl)
		{
			pcl->next = top;
			top = pcl;
		}
		inline void clear(void) { top = nullptr; }
		inline Particle *first(void) { return top; }
		inline Particle *next(Particle *prev) { return prev->next; }
		void print_pcl(void)
		{
			for (Particle *pcl_tmp = top; pcl_tmp; pcl_tmp = pcl_tmp->next)
			{
				std::cout << "pcl x: " << pcl_tmp->x << std::endl;
			}
		}
	};

	struct Node_s
	{
		double x;
		double a, v, dv, du;
		double m, ma, mv;
	};

	// Node in the middle of the element
	struct Node_f
	{
		double x;
		double pa, pv, p, dpv, dp;
		double m, mpa, mpv, mp;
	};

	struct BodyForce
	{
		double bf_x;
		size_t pcl_id;
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

	// Pore pressure "acceleration" boundary conditions (Dirchlet boundary conditions)
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
	double u;
	size_t node_id;
	};

	// Pore pressure "velocity" initial boundary conditions
	struct PVelocityBC
	{
		double pv;
		size_t node_id;
	};

	// Pore pressure initial boundary conditions
	struct PorePressureBC
	{
	double p;
	size_t node_id;
	};

	int hyme1D_MPM(void);
}

#endif // _HYME1D_MPM_H_