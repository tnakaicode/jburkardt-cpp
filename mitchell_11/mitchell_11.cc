#include "myfunctions.h"

template<int dim>
class Problem11
{
	public:
		Problem11();
		
		void run();
	
		ForcingFunction<dim> *forcing_function;
        ExactSolution<dim>   *exact_solution;
        BoundaryFunction<dim>  *boundary_function;		
		
	private:
        void read_inputs();
		void setup_geometry(int cycle);
		void assemble();
		void solve();
		void refine_grid(int cycle, const char* s);
		void output_results(int cycle, const char* s);
		void calculate_error(int cycle, ConvergenceTable &error_table);        
		void print_errors(bool convergence_study, ConvergenceTable error_table);    
        
		struct InputFile
        {
            int        nx, ny, nc_con, nc_ad;
            double     a, alpha, beta, p1, p2, p3, p4;
            bool       convergence_test, adaptive_grad, adaptive_hess;
        };
        
        InputFile              input;        
		Triangulation<dim>     mesh;
		FESystem<dim>          fe;
		DoFHandler<dim>        dof_handler;
		
		ConstraintMatrix       constraints; //for boundary conditions
		SparsityPattern        sparsity_pattern;
		SparseMatrix<double>   system_matrix;
		Vector<double>         system_rhs;
		Vector<double>         solution;

		ConvergenceTable	   convergence_table;
		ConvergenceTable	   adaptive_error_g_table;
		ConvergenceTable	   adaptive_error_h_table;	
};

template<int dim>
Problem11<dim>::Problem11():
    fe(FE_Q<dim>(2), 1),
    dof_handler(mesh)
{}

template<int dim>
void Problem11<dim>::read_inputs()
{
    std::ifstream inputFile;
    std::string   line;
    std::string   tmp;
    
    inputFile.open("input.in");
    
    while (getline(inputFile, line))
    {
        //commented lines start with #
        if (!line.find("#")==0)
        {            
            if (line.find("nx")==0)
            {
                tmp      = line.substr(3);
                input.nx = atoi(tmp.c_str());
            }
            else if (line.find("ny")==0)
            {
                tmp = line.substr(3);
                input.ny = atoi(tmp.c_str());
            }
            else if (line.find("alpha")==0)
            {
                tmp      = line.substr(6);
                input.alpha = atof(tmp.c_str());
            }
            else if (line.find("beta")==0)
            {
                tmp = line.substr(5);
                input.beta = atof(tmp.c_str());
            }
            else if (line.find("p1")==0)
            {
                tmp = line.substr(3);
                input.p1 = atof(tmp.c_str());
            }
            else if (line.find("p2")==0)
            {
                tmp = line.substr(3);
                input.p2 = atof(tmp.c_str());
            }
            else if (line.find("p3")==0)
            {
                tmp = line.substr(3);
                input.p3 = atof(tmp.c_str());
            }
            else if (line.find("p4")==0)
            {
                tmp = line.substr(3);
                input.p4 = atof(tmp.c_str());
            }
            else if (line.find("convergence_test")==0)
            {
                tmp = line.substr(17);
                input.convergence_test = (atoi(tmp.c_str()) == 1);
            }
            else if (line.find("adaptive_gradient")==0)
            {
                tmp = line.substr(18);
                input.adaptive_grad = (atoi(tmp.c_str()) == 1);
            }
            else if (line.find("adaptive_hessian")==0)
            {
                tmp = line.substr(17);
                input.adaptive_hess = (atoi(tmp.c_str()) == 1);
            }
            else if (line.find("convergence_cycles")==0)
            {
                tmp = line.substr(19);
                input.nc_con = atoi(tmp.c_str());
            }
            else if (line.find("adaptive_cycles")==0)
            {
                tmp = line.substr(16);
                input.nc_ad = atoi(tmp.c_str());
            }
            else if (line.find("a")==0)
            {
                tmp      = line.substr(2);
                input.a = atof(tmp.c_str());
            }
        }
    }    
    inputFile.close();    
}

template<int dim>
void Problem11<dim>::setup_geometry (int cycle)
{  		
	if (cycle == 0)
	{
		mesh.clear();
		std::vector<unsigned int> number_elements(2);
		number_elements[0] = input.nx-1;
		number_elements[1] = input.ny-1;
  
		GridGenerator::subdivided_hyper_rectangle(mesh, number_elements,
			Point<dim>(-1, -1), Point<dim>(1, 1), false);
	}
  
    std::printf("Number of active cells:%d\n", mesh.n_active_cells());
  
    dof_handler.distribute_dofs(fe);  
    std::printf("Number of degrees of freedom:%d\n", dof_handler.n_dofs()); 
		
    constraints.clear ();    
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);			   
    VectorTools::interpolate_boundary_values(dof_handler, 0, *boundary_function, constraints);	   
    constraints.close();
            
	//calculate sparsity pattern 
	DynamicSparsityPattern c_sparsity(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, c_sparsity, constraints, false);
    
    constraints.condense(c_sparsity);
	sparsity_pattern.copy_from(c_sparsity);
    
    system_matrix.reinit(sparsity_pattern);
    system_rhs.reinit(dof_handler.n_dofs()); 
    solution.reinit(dof_handler.n_dofs());         
}

template<int dim>
void Problem11<dim>::assemble()
{
	//Assemble stiffness matrix and rhs, apply boundary conditions
	printf("Assembling matrix and right hand side...");
	Timer timer;
	timer.start();
	
	const QGauss<dim> quadrature_formula(3);
	FEValues<dim> fe_values (fe, quadrature_formula, update_values 
			| update_gradients | update_quadrature_points | update_JxW_values);
			
	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_q_points = quadrature_formula.size();
	
	FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
	Vector<double>     cell_rhs(dofs_per_cell);
	
	std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

	typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
												   endc = dof_handler.end();
	for (; cell!=endc; ++cell)
	{
		cell_matrix = 0;
		cell_rhs = 0;
		fe_values.reinit (cell);

		for (unsigned int q_index=0; q_index<n_q_points; q_index++)
		{
			for (unsigned int i=0; i<dofs_per_cell; i++)
			{
				Point<dim> q_point = fe_values.quadrature_point(q_index);
				double x = q_point[0];
				double y = q_point[1];
				double p;
				
				if (y >=0 && x >= 0)
				{
					p = input.p1;
				}
				else if (y >= 0 && x < 0)
				{
					p = input.p2;
				}
				else if (y < 0 && x < 0)
				{
					p = input.p3;
				}
				else
				{
					p = input.p4;
				}
				
				for (unsigned int j=0; j<dofs_per_cell; j++)
				{
					cell_matrix(i,j) += (p*fe_values.shape_grad(i,q_index)*fe_values.shape_grad(j,q_index))
								*fe_values.JxW(q_index);
	
				}
				
				cell_rhs(i) += forcing_function->value(fe_values.quadrature_point(q_index))
											*fe_values.shape_value(i,q_index)*fe_values.JxW(q_index);
			}
		}
		
		cell->get_dof_indices (local_dof_indices);
		constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices,
								system_matrix, system_rhs);
	}
	
    timer.stop();
    printf("done (%gs)\n", timer());                           
}

template<int dim>
void Problem11<dim>::solve()
{
	Timer timer;
	printf("Solving linear system... ");
	timer.start ();
	
	SparseDirectUMFPACK A_direct;
	A_direct.initialize(system_matrix);

	A_direct.vmult(solution, system_rhs);
    
    constraints.distribute(solution);
    timer.stop ();
	
    printf("done (%gs)\n", timer());
}


template<int dim>
void Problem11<dim>::refine_grid(int cycle, const char* s)
{
	if (cycle == 0)
	{
		setup_geometry(0);
	}
	else
	{	
		printf("Refining mesh...");
		
		Timer timer;
		timer.start ();	
		Vector<float> estimated_error_per_cell(mesh.n_active_cells());
		
		if (!strcmp(s, "gradient"))
		{
			DerivativeApproximation::approximate_gradient(dof_handler, solution, estimated_error_per_cell);
			typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
														   endc = dof_handler.end();
														   
			for(unsigned int cell_no=0; cell!=endc; ++cell, ++cell_no)
			{
				estimated_error_per_cell(cell_no)*=std::pow(cell->diameter(), 1+1.0*dim/2);
			}			
		}
		else
		{		
			
			KellyErrorEstimator<dim>::estimate(dof_handler, QGauss<1>(3), 
					typename FunctionMap<dim>::type(), solution, estimated_error_per_cell);			
		}
		
		GridRefinement::refine_and_coarsen_fixed_number(mesh, estimated_error_per_cell,
					0.3, 0.03);
		
		mesh.execute_coarsening_and_refinement();
		
		timer.stop();
		printf("done (%gs)\n",timer());
		
		setup_geometry(cycle);
	}	
}

template<int dim>
void Problem11<dim>::output_results(int cycle, const char *s)
{	
			
	DataOut<dim> data_out;
	data_out.attach_dof_handler (dof_handler);
	data_out.add_data_vector (solution, "solution");
	data_out.build_patches ();
	
	std::ostringstream filename;
	filename << s << cycle << ".vtk";	
	
	std::ofstream output (filename.str().c_str());
	data_out.write_vtk (output);
}

template<int dim>
void Problem11<dim>::calculate_error(int cycle, ConvergenceTable &error_table)
{	
    Vector<float> difference_per_cell (mesh.n_active_cells());
    VectorTools::integrate_difference (dof_handler,
                                       solution,
                                       *exact_solution,
                                       difference_per_cell,
                                       QGauss<dim>(3),
                                       VectorTools::L2_norm);
    const double L2_error = difference_per_cell.l2_norm();

	const QTrapez<1> q_trapez;
	const QIterated<dim> q_iterated (q_trapez, 5);
	VectorTools::integrate_difference (dof_handler,
										solution,
										*exact_solution,
										difference_per_cell,
										q_iterated,
										VectorTools::Linfty_norm);
	const double Linfty_error = difference_per_cell.linfty_norm();

	error_table.add_value("cycle", cycle);
	error_table.add_value("cells", mesh.n_active_cells());
	error_table.add_value("dofs", dof_handler.n_dofs());
	error_table.add_value("L2", L2_error);
	error_table.add_value("Linfty", Linfty_error);
}

template<int dim>
void Problem11<dim>::run()
{	
	read_inputs(); 
	
	boundary_function->a = input.a;	
	boundary_function->alpha = input.alpha;	
	boundary_function->beta = input.beta;	
	boundary_function->p1 = input.p1;	
	boundary_function->p2 = input.p2;	
	boundary_function->p3 = input.p3;
	boundary_function->p4 = input.p4;	
	exact_solution->a = input.a;	
	exact_solution->alpha = input.alpha;	
	exact_solution->beta = input.beta;		
	exact_solution->p1 = input.p1;	
	exact_solution->p2 = input.p2;	
	exact_solution->p3 = input.p3;
	exact_solution->p4 = input.p4;
	
	printf("a = %f\n", input.a);
	printf("alpha = %f\n", input.alpha);
	printf("beta = %f\n", input.beta);
	printf("p1 = %f\n", input.p1);
	printf("p2 = %f\n", input.p2);
	printf("p3 = %f\n", input.p3);
	printf("p4 = %f\n", input.p4);
	
	if (input.convergence_test)
	{		
		printf("\nBeginning convergence study...\n\n");
		
		for (int cycle = 0; cycle <=input.nc_con; cycle++)
		{			
			Timer timer;
			timer.start();
				
			if (cycle!=0)
			{
				mesh.refine_global(1);
			}
			
			printf("Cycle %i\n", cycle);	
			
			setup_geometry(cycle);			
			assemble();
			solve();			
			output_results(cycle, "convergence");	
			
			calculate_error(cycle, convergence_table);
			
			timer.stop();			
			convergence_table.add_value("Time", timer());
		}
		
		print_errors(true, convergence_table);	
	}
	if (input.adaptive_grad)
	{
		printf("\nBeginning adaptive meshing using gradient error estimate...\n\n");
		
		for (int cycle = 0; cycle <= input.nc_ad; cycle++)
		{
			Timer timer;
			timer.start();
			
			printf("Cycle %i\n", cycle);
			
			refine_grid(cycle, "gradient");									
			assemble();
			solve();
			output_results(cycle, "adaptive_grad");
			
			calculate_error(cycle, adaptive_error_g_table);
			timer.stop();			
			adaptive_error_g_table.add_value("Time", timer());
		}
		
		print_errors(false, adaptive_error_g_table);			
	}	
	
	if (input.adaptive_hess)
	{
		printf("\nBeginning adaptive meshing using hessian error estimate...\n\n");
		
		for (int cycle = 0; cycle <= input.nc_ad; cycle++)
		{
			Timer timer;
			timer.start();
			
			printf("Cycle %i\n", cycle);
			
			refine_grid(cycle, "hessian");									
			assemble();
			solve();
			output_results(cycle, "adaptive_hess");
			
			calculate_error(cycle, adaptive_error_h_table);
			timer.stop();			
			adaptive_error_h_table.add_value("Time", timer());
		}
		
		print_errors(false, adaptive_error_h_table);		
	}
}

template<int dim>
void Problem11<dim>::print_errors(bool convergence_study, ConvergenceTable error_table)
{
	error_table.set_precision("L2", 3);
	error_table.set_precision("Linfty", 3);
	error_table.set_scientific("L2", true);
	error_table.set_scientific("Linfty", true);
	
	if (convergence_study)
	{
		error_table.evaluate_convergence_rates("L2", ConvergenceTable::reduction_rate_log2);
	}
	
	printf("\n");
	printf("Error analysis:\n");
	error_table.write_text(std::cout);
}

int main ()
{
	const int dim = 2;
	
	Problem11<dim> problem;
	ForcingFunction<dim>    ff;
	ExactSolution<dim>      ex;
	BoundaryFunction<dim>   bf;
	
	problem.forcing_function   = &ff;
	problem.exact_solution     = &ex;
	problem.boundary_function  = &bf;
    problem.run();
}
