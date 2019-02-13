#ifndef _HYME1D_FEM2_H_
#define _HYME1D_FEM2_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

namespace hyme1D_FEM2
{
	// Gauss point
	struct GaussPoint
	{
		double n;

		double s11;
		double e11;

		double xi; // gauss coordinates
		double w; // weight

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
		double density_avg;
	};


	struct Node_s
	{
		double x;
		double du; // used to update stress
		double a, v, u;
	};

	// Node in the middle of the element
	struct Node_f
	{
		double x;
		double pa, pv, p;
	};

	struct Element
	{
		// topology
		size_t index_s_x; // node index of solid phase
		size_t index_f_x; // node index of fluid phase 

		double density_s, density_f;
		double E;
		double k; // permeability
		double miu; // dynamic viscosity

					/*
					All physical parameters needed to be
					updated are contained in gauss points.
					*/
		GaussPoint gps1, gps2, gpf;

		// Temporary variables
		double k_div_miu; // k divided by miu
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
		double u;
		size_t node_id;
	};

	// pore pressure initial boundary conditions
	struct PorePressureBC
	{
		double p;
		size_t node_id;
	};

	int hyme1D_FEM2(void);
}

#endif