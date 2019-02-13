#include <iostream>

#include "utility_func.h"

void dis_mat(double *mat, size_t row, size_t col)
{
	for (size_t i = 0; i < row; i++)
	{
		for (size_t j = 0; j < col; j++)
		{
			std::cout << mat[i * col + j] << " ";
		}
		std::cout << std::endl;
	}
}

void dis_vec(double *vec, size_t num)
{
	for (size_t i = 0; i < num; i++)
	{
		std::cout << vec[i] << " ";
	}
	std::cout << std::endl;
}