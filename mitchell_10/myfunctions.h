#ifndef MYFUNCTIONS_H
#define MYFUNCTIONS_H

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/derivative_approximation.h>

#include <fstream>
#include <cmath>


using namespace dealii;

/*******************************************************************************
 * Define forcing function, boundary values and exact solution
 ******************************************************************************/
template<int dim>
class ForcingFunction : public Function<dim>
{
    public:
		double a, b;
        ForcingFunction() : Function<dim>(1) {}
        virtual double value(const Point<dim> & point,
						const unsigned int component = 0 ) 	const;
};

template<int dim>
double ForcingFunction<dim>::value(const Point<dim> & point,
						const unsigned int ) 	const
{  
	double x = point[0];
	double y = point[1];
	
	double f;
	
	if (x <= b*(y + 1))
	{
		f = pow(M_PI,2)*cos(M_PI*y/2)/4;
	}
	else
	{
		f = pow(M_PI,2)*cos(M_PI*y/2)/4 - a*(a-1)*pow(x - b*(y+1), a-2)*(1 + pow(b,2));
	}
	
	return f;
}

template<int dim>
class BoundaryFunction : public Function<dim>
{
    public:
		double a, b;
        BoundaryFunction() : Function<dim>(1) {}
        virtual double value(const Point<dim> & point,
						const unsigned int component = 0 ) 	const;
};

template<int dim>
double BoundaryFunction<dim>::value(const Point<dim> & point,
						const unsigned int ) 	const
{  
	double x = point[0];
	double y = point[1];
	
	double value;
	
	if (x <= b*(y + 1))
	{
		value = cos(M_PI*y/2);
	}
	else
	{
		value = cos(M_PI*y/2) + pow(x - b*(y+1), a);
	}
	
	return value;
}


template<int dim>
class ExactSolution : public Function<dim>
{
    public:
		double a, b;
        ExactSolution() : Function<dim>(1) {};       
        virtual double value(const Point<dim> & point,
						const unsigned int component = 0 ) 	const;
		virtual Tensor<1,dim> gradient (const Point<dim> &point,
						const unsigned int component = 0) const;
};

template<int dim>
double ExactSolution<dim>::value(const Point<dim> & point,
						const unsigned int ) 	const
{
	double x = point[0];
	double y = point[1];
	
	double u;
	
	if (x <= b*(y + 1))
	{
		u = cos(M_PI*y/2);
	}
	else
	{
		u = cos(M_PI*y/2) + pow(x - b*(y+1), a);
	}
	
	return u;
}


template <int dim>
Tensor<1,dim> ExactSolution<dim>::gradient (const Point<dim> &point,
			const unsigned int) const
{
	double x = point[0];
	double y = point[1];

	Tensor<1,dim> grad_u;
	
	if (x <= b*(y+1))
	{
		grad_u[0] = 0;
		grad_u[1] = -M_PI*sin(M_PI*y/2)/2;
	}
	else
	{
		grad_u[0] = a*pow(x-b*(y+1),a-1);
		grad_u[1] = -M_PI*sin(M_PI*y/2)/2 - a*b*pow(x - b*(y+1), a-1);
	}
	
	return grad_u;
}


#endif
