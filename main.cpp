#include <cstdlib>
#include <iostream>
#include <vector>
#include <list>

#include "hyme1D_FEM.h"
#include "hyme1D_FEM2.h"
#include "hyme1D_FEM3.h"
#include "hyme1D_FEM4.h"
#include "hyme1D_MPM.h"


int main(int argc, void** argv)
{
	//hyme1D_FEM::hyme1D_FEM();
	//hyme1D_FEM2::hyme1D_FEM2();
	//hyme1D_FEM3::hyme1D_FEM3(); // bad
	hyme1D_FEM4::hyme1D_FEM4();

	//hyme1D_MPM::hyme1D_MPM();

	//test_list();
	//eigen_test2();

	system("pause");
	return 0;
}


void test_list(void)
{
	std::list<int> int_list;
	std::list<int>::iterator iter;

	for (size_t i = 0; i < 10; i++)
	{
		int_list.push_back(i);
	}
	for (iter = int_list.begin(); iter != int_list.end(); iter++)
	{
		std::cout << *iter << " ";
	}
	std::cout << std::endl;

	iter = int_list.begin();
	iter++;
	std::cout << "del: " << *iter << std::endl;
	iter = int_list.erase(iter);
	std::cout << "now: " << *iter << std::endl;
	iter = int_list.erase(iter);
	iter++;
	iter = int_list.erase(iter);
	iter = int_list.erase(iter);


	for (iter = int_list.begin(); iter != int_list.end(); iter++)
	{
		std::cout << *iter << " ";
	}

}

void eigen_test(void)
{
	Eigen::Matrix<double, 4, 4> A;
	Eigen::Matrix<double, 4, 1> b;
	A << 1, 1, 2, 3,
		3, -1, -1, -2,
		2, 3, -1, -1,
		1, 2, 3, -1;
	b << 1, -4, -6, -4;
	std::cout << "Here is the matrix A:\n" << A << std::endl;
	std::cout << "Here is the right hand side b:\n" << b << std::endl;
	Eigen::Matrix<double, 4, 1> x = A.colPivHouseholderQr().solve(b);
	std::cout << "The solution is:\n" << x << std::endl;

	Eigen::SparseMatrix<double> a2(4, 4);
	std::vector<Eigen::Triplet<double> > coefficients;

	Eigen::VectorXd b2(4);
	Eigen::VectorXd x2(4);

	coefficients.push_back(Eigen::Triplet<double>(0, 0, 1.0));
	coefficients.push_back(Eigen::Triplet<double>(0, 1, 1.0));
	coefficients.push_back(Eigen::Triplet<double>(0, 2, 2.0));
	coefficients.push_back(Eigen::Triplet<double>(0, 3, 3.0));
	coefficients.push_back(Eigen::Triplet<double>(1, 0, 3.0));
	coefficients.push_back(Eigen::Triplet<double>(1, 1, -1.0));
	coefficients.push_back(Eigen::Triplet<double>(1, 2, -1.0));
	coefficients.push_back(Eigen::Triplet<double>(1, 3, -2.0));
	coefficients.push_back(Eigen::Triplet<double>(2, 0, 2.0));
	coefficients.push_back(Eigen::Triplet<double>(2, 1, 3.0));
	coefficients.push_back(Eigen::Triplet<double>(2, 2, -1.0));
	coefficients.push_back(Eigen::Triplet<double>(2, 3, -1.0));
	coefficients.push_back(Eigen::Triplet<double>(3, 0, 1.0));
	coefficients.push_back(Eigen::Triplet<double>(3, 1, 2.0));
	coefficients.push_back(Eigen::Triplet<double>(3, 2, 3.0));
	coefficients.push_back(Eigen::Triplet<double>(3, 3, -1.0));
	//coefficients.push_back(Eigen::Triplet<double>(3, 3, -2.0));
	//coefficients.push_back(Eigen::Triplet<double>(3, 3, 3.0));


	a2.setFromTriplets(coefficients.begin(), coefficients.end());
	a2.makeCompressed();

	//b2 << 1.0, -4.0, -6.0, -4.0;
	for (size_t i = 0; i < 4; i++)
	{
		b2(i) = (double)i;
	}

	Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<Eigen::Index> > solver;
	solver.analyzePattern(a2);
	solver.factorize(a2);
	x2 = solver.solve(b2);
	std::cout << "Here is the matrix A:\n" << a2 << std::endl;
	std::cout << "Here is the right hand side b:\n" << b2 << std::endl;
	std::cout << "The solution is:\n" << x2 << std::endl;

	return;
}

void eigen_test2(void)
{
	Eigen::Matrix<double, 4, 4> A;
	Eigen::Matrix<double, 4, 1> b;

	Eigen::VectorXd vec(5);

	A.setZero();
	b.setConstant(3.0);
	vec.setOnes();

	std::cout << A << std::endl;
	std::cout << b << std::endl;
	std::cout << vec << std::endl;
}
