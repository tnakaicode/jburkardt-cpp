#ifndef MYFUNCTIONS_H
#define MYFUNCTIONS_H

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
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

inline double calculate_theta(double x, double y)
{
	const double EPSILON = 1e-10;
	double theta;
	
	double phi = atan(fabs(y/x));
	
	if (fabs(x) < EPSILON && fabs(y) < EPSILON)
	{
		theta = 0;
	}
	else if (fabs(y) < EPSILON && x > 0)
	{
		theta = 0;
	}
	else if (fabs(x) < EPSILON && y > 0)
	{
		theta = M_PI/2;
	}
	else if (fabs(y) < EPSILON && x < 0)
	{
		theta = M_PI;
	}
	else if (fabs(x) < EPSILON && y < 0)
	{
		theta = 3*M_PI/2;
	}
	else if (x >= 0 && y >= 0)
	{
		theta = phi;
	}
	else if (y >= 0 && x <= 0)
	{
		theta = M_PI - phi;
	}
	else if (y < 0 && x < 0)
	{
		theta = M_PI + phi;
	}
	else
	{
		theta = 2*M_PI - phi;
	}
	
	return theta;
}

/*******************************************************************************
 * Define boundary values and exact solution
 ******************************************************************************/

template<int dim>
class BoundaryFunction : public Function<dim>
{
    public:
		double alpha;
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
	double r = sqrt(pow(x,2) + pow(y,2));
	double theta = calculate_theta(x,y);
	
	double value = pow(r,alpha)*sin(alpha*theta);
	
	return value;
}


template<int dim>
class ExactSolution : public Function<dim>
{
    public:
		double alpha;
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
	double r = sqrt(pow(x,2) + pow(y,2));
	double theta = calculate_theta(x,y);
	
	double value = pow(r,alpha)*sin(alpha*theta);
	
	return value;
}


template <int dim>
Tensor<1,dim> ExactSolution<dim>::gradient (const Point<dim> &point,
			const unsigned int) const
{
	double x = point[0];
	double y = point[1];
	double r = sqrt(pow(x,2) + pow(y,2));
	double theta = calculate_theta(x,y);
	Tensor<1,dim> grad_u;	
	
	grad_u[0] = alpha*pow(r,alpha-2)*(x*sin(alpha*theta) - y*cos(alpha*theta));
	grad_u[1] = alpha*pow(r,alpha-2)*(y*sin(alpha*theta) + x*cos(alpha*theta));
	
	return grad_u;
}


#endif
