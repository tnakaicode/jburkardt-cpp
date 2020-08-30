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
 * Define forcing function, boundary values and exact solution
 ******************************************************************************/
template<int dim>
class ForcingFunction : public Function<dim>
{
    public:
        ForcingFunction() : Function<dim>(1) {}
        virtual double value(const Point<dim> & point,
						const unsigned int component = 0 ) 	const;
};

template<int dim>
double ForcingFunction<dim>::value(const Point<dim> & point,
						const unsigned int ) 	const
{  
	return 0;
}

template<int dim>
class BoundaryFunction : public Function<dim>
{
    public:
		double a, alpha, beta, p1, p2, p3, p4;
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
	double theta, value, mu;
	
	theta = calculate_theta(x,y);
	
	if (theta <= M_PI/2)
	{
		mu = cos(a*(M_PI/2 - beta))*cos(a*(theta-M_PI/2+alpha));
	}
	else if (theta <= M_PI)
	{		
		mu = cos(alpha*a)*cos(a*(theta - M_PI + beta));
	}
	else if (theta <= 3*M_PI/2)
	{
		mu = cos(beta*a)*cos(a*(theta - M_PI - alpha));
	}
	else
	{
		mu = cos(a*(M_PI/2-alpha))*cos(a*(theta-3*M_PI/2-beta));
	}
	
	value = pow(r,a)*mu;
	
	return value;
}


template<int dim>
class ExactSolution : public Function<dim>
{
    public:
		double a, alpha, beta, p1, p2, p3, p4;
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
	double theta, u, mu;
	
	theta = calculate_theta(x,y);
	
	if (theta <= M_PI/2)
	{
		mu = cos(a*(M_PI/2 - beta))*cos(a*(theta-M_PI/2+alpha));
	}
	else if (theta <= M_PI)
	{
		mu = cos(alpha*a)*cos(a*(theta - M_PI + beta));
	}
	else if (theta <= 3*M_PI/2)
	{
		mu = cos(beta*a)*cos(a*(theta-M_PI-alpha));
	}
	else
	{
		mu = cos(a*(M_PI/2-alpha))*cos(a*(theta-3*M_PI/2-beta));
	}
	
	u = pow(r,a)*mu;
	
	return u;
}

template <int dim>
Tensor<1,dim> ExactSolution<dim>::gradient (const Point<dim> &point,
			const unsigned int) const
{
//double x = point[0];
//double y = point[1];

	Tensor<1,dim> grad_u;
	
	grad_u[0] = 0;
	grad_u[1] = 0;
	
	return grad_u;
}

#endif
