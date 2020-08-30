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
#include <deal.II/base/table_handler.h>
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
#include <deal.II/lac/affine_constraints.h>
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
 * Define forcing function and exact solution
 ******************************************************************************/
template<int dim>
class ForcingFunction : public Function<dim>
{
    public:
		double p;
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
	
	if (abs(p) < 1e-14)
	{
		f = 0;
	}     
	else if (abs(p - 1) < 1e-14)
	{
		f = 32*(y*(1-y) + x*(1-x));
	}
	else
	{
		f = -pow(2,4*p)*p*(pow(y,p)*pow(1-y,p)*((p-1)*pow(x,p-2)*pow(1-x,p) 
					- 2*p*pow(x,p-1)*pow(1-x,p-1) + (p-1)*pow(x,p)*pow(1-x,p-2))
					+ pow(x,p)*pow(1-x,p)*((p-1)*pow(y,p-2)*pow(1-y,p)
					- 2*p*pow(y,p-1)*pow(1-y,p-1) + (p-1)*pow(y,p)*pow(1-y,p-2)));
	
	}
	
	return f;
}

template<int dim>
class ExactSolution : public Function<dim>
{
    public:
		double p;
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
	
	u = pow(2,4*p)*pow(x,p)*pow(1-x,p)*pow(y,p)*pow(1-y,p);
	
	return u;
}


template <int dim>
Tensor<1,dim> ExactSolution<dim>::gradient (const Point<dim> &point,
			const unsigned int) const
{
	double x = point[0];
	double y = point[1];
	Tensor<1,dim> grad_u;
	
	if (abs(p) < 1e-14)
	{
		grad_u[0] = 1;
		grad_u[1] = 1;
	}
	else
	{
		grad_u[0] = pow(2,4*p)*pow(y,p)*pow(1-y,p)*(p*pow(x,p-1)*pow(1-x,p) 
				- p*pow(x,p)*pow(1-x,p-1));
		grad_u[1] = pow(2,4*p)*pow(x,p)*pow(1-x,p)*(p*pow(y,p-1)*pow(1-y,p) 
				- p*pow(y,p)*pow(1-y,p-1));
	}
	
	return grad_u;
}


#endif
