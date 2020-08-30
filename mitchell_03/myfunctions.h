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

const double nu = 0.3;
const double E = 1;

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
		int mode;
        BoundaryFunction() : Function<dim>(2) {}
        virtual void vector_value(const Point<dim> &p, 
                Vector<double> &values) const;
};

template<int dim>
void BoundaryFunction<dim>::vector_value(const Point<dim> &point, 
                Vector<double> &values) const
{  
	const double G = E/(2*(1+nu));
	const double k = 3 - 4*nu;
		
	double x = point[0];
	double y = point[1];
	double r = sqrt(pow(x,2) + pow(y,2));
	double theta = calculate_theta(x,y);
	double lambda, Q;
	
	if (mode == 1)
	{
		lambda = 0.5444837367825;
		Q = 0.5430755788367;
		
		values(0) = 1/(2*G)*pow(r,lambda)*((k - Q*(lambda+1))*cos(lambda*theta)-lambda*cos((lambda-2)*theta));
		values(1) = 1/(2*G)*pow(r,lambda)*((k + Q*(lambda+1))*sin(lambda*theta)+lambda*sin((lambda-2)*theta));
	}
	else
	{
		lambda = 0.9085291898461;
		Q = -0.2189232362488;

		values(0) = 1/(2*G)*pow(r,lambda)*((k - Q*(lambda+1))*cos(lambda*theta)-lambda*cos((lambda-2)*theta));
		values(1) = -1/(2*G)*pow(r,lambda)*((k + Q*(lambda+1))*sin(lambda*theta)+lambda*sin((lambda-2)*theta));
	}
}


template<int dim>
class ExactSolution : public Function<dim>
{
    public:
		int mode;
        ExactSolution() : Function<dim>(2) {};       
        virtual void vector_value(const Point<dim> &p, 
                Vector<double> &values) const;
		virtual void vector_gradient(const Point<dim> &p, 
                std::vector<Tensor<1,dim> > &gradients) const;
};

template<int dim>
void ExactSolution<dim>::vector_value(const Point<dim> &point, 
                Vector<double> &values) const
{  
	const double G = E/(2*(1+nu));
	const double k = 3 - 4*nu;
	
	double x = point[0];
	double y = point[1];
	double r = sqrt(pow(x,2) + pow(y,2));
	double theta = calculate_theta(x,y);
	double lambda, Q;
	
	if (mode == 1)
	{
		lambda = 0.5444837367825;
		Q = 0.5430755788367;
		
		values(0) = 1/(2*G)*pow(r,lambda)*((k - Q*(lambda+1))*cos(lambda*theta)-lambda*cos((lambda-2)*theta));
		values(1) = 1/(2*G)*pow(r,lambda)*((k + Q*(lambda+1))*sin(lambda*theta)+lambda*sin((lambda-2)*theta));
	}
	else
	{
		lambda = 0.9085291898461;
		Q = -0.2189232362488;

		values(0) = 1/(2*G)*pow(r,lambda)*((k - Q*(lambda+1))*cos(lambda*theta)-lambda*cos((lambda-2)*theta));
		values(1) = -1/(2*G)*pow(r,lambda)*((k + Q*(lambda+1))*sin(lambda*theta)+lambda*sin((lambda-2)*theta));
	}
};


template<int dim>
void ExactSolution<dim>::vector_gradient(const Point<dim> &point, 
                std::vector<Tensor<1,dim> > & gradients) const
{	
    gradients[0][0] = 0;
    gradients[0][1] = 0;
    
    gradients[1][0] = 0;
    gradients[1][1] = 0;
}



#endif
