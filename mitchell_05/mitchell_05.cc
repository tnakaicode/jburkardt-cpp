#include "myfunctions.h"

template<int dim>
class Problem5
{
	public:
		Problem5();
		
		void run();
		
	private:
        void read_inputs();
		void setup_geometry(int cycle);
		void assemble();
		void solve();
		void refine_grid(int cycle, const char* s);
		void output_results(int cycle, const char* s);
              
        struct InputFile
        {
            int        nc_con, nc_ad;
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
};

template<int dim>
Problem5<dim>::Problem5():
    fe(FE_Q<dim>(2), 1),
    dof_handler(mesh)
{}

template<int dim>
void Problem5<dim>::read_inputs()
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
			if (line.find("convergence_test")==0)
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
        }
    }    
    inputFile.close();    
}

template<int dim>
void Problem5<dim>::setup_geometry (int cycle)
{  		
	if (cycle == 0)
	{
		mesh.clear();
		GridIn<dim> grid_in;
		grid_in.attach_triangulation(mesh);
		std::ifstream input_file("mesh.msh");
		grid_in.read_msh(input_file);		
	}
  
    std::printf("Number of active cells:%d\n", mesh.n_active_cells());
  
    dof_handler.distribute_dofs(fe);  
    std::printf("Number of degrees of freedom:%d\n", dof_handler.n_dofs()); 
		
    constraints.clear ();    
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);			     
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
void Problem5<dim>::assemble()
{
	//Assemble stiffness matrix and rhs, apply boundary conditions
	printf("Assembling matrix and right hand side...");
	Timer timer;
	timer.start();
	
	const QGauss<dim> quadrature_formula(3);
	const QGauss<dim-1> face_quadrature_formula(3);
	
	FEValues<dim> fe_values (fe, quadrature_formula, update_values 
			| update_gradients | update_quadrature_points | update_JxW_values);
	
	
	FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula, update_values 
			| update_quadrature_points | update_gradients | update_JxW_values);
		
	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_q_points_cell = quadrature_formula.size();
	const unsigned int n_q_points_face = face_quadrature_formula.size();
	
	FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
	Vector<double> cell_rhs (dofs_per_cell);

	std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

	typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
												   endc = dof_handler.end();
	for (; cell!=endc; ++cell)
	{
		cell_matrix = 0;
		cell_rhs = 0;
		fe_values.reinit (cell);

		int k = cell->material_id();		
		double p, q, f;
		
		switch (k)
		{
			case 1: //region 1
				p = 25;
				q = 25;
				f = 0;
				break;
				
			case 2: //region 2
			case 3:
				p = 7;
				q = 0.8;
				f = 1;
				break;
				
			case 4: //region 3
				p = 5;
				q = 0.0001;
				f = 1;
				break;
			
			case 5://region 4
				p = 0.2;
				q = 0.2;
				f = 0;
				break;
				
			case 6://region 5
			case 7:
				p = 0.05;
				q = 0.05;
				f = 0;
				break;						
		}
		
		for (unsigned int q_index=0; q_index<n_q_points_cell; q_index++)
		{			
			for (unsigned int i=0; i<dofs_per_cell; i++)
			{
				Tensor<1,dim> grad_phi_i = fe_values.shape_grad(i,q_index);
				double phi_i_x = grad_phi_i[0];
				double phi_i_y = grad_phi_i[1];
				
				for (unsigned int j=0; j<dofs_per_cell; j++)
				{
					Tensor<1,dim> grad_phi_j = fe_values.shape_grad(j,q_index);
					double phi_j_x = grad_phi_j[0];
					double phi_j_y = grad_phi_j[1];
					
					cell_matrix(i,j) += (p*phi_i_x*phi_j_x + q*phi_i_y*phi_j_y)
								*fe_values.JxW(q_index);
						
				}
				
				cell_rhs(i) += f*fe_values.shape_value(i,q_index)*fe_values.JxW(q_index);
			}
		}
		
		//add line integral term to matrix and rhs	
		for (unsigned int face=0; face < GeometryInfo<dim>::faces_per_cell; face++)
		{						
			if (cell->face(face)->at_boundary())
			{
				fe_face_values.reinit (cell, face);
							
				for (unsigned int q_boundary = 0; q_boundary < n_q_points_face; q_boundary++)
				{
					for (unsigned int i=0; i<dofs_per_cell; i++)
					{
						double phi_i = fe_face_values.shape_value(i, q_boundary);
						
						for (unsigned int j=0; j<dofs_per_cell; j++)
						{
							Tensor<1,dim> grad_phi_j = fe_face_values.shape_grad(j,q_boundary);
							double phi_j_x = grad_phi_j[0];
							double phi_j_y = grad_phi_j[1];
							double phi_j = fe_face_values.shape_value(j, q_boundary);
							
							int boundary_number = cell->face(face)->boundary_id();							
							int c, g;
							
							switch (boundary_number)
							{
								case 1://bottom
									c = 3;
									g = 1;
									cell_matrix(i,j) += (-p*phi_j_x*phi_i - c*phi_i*phi_j)
															*fe_face_values.JxW(q_boundary);
									break;
								
								case 2://right
									c = 2;
									g = 2;
									cell_matrix(i,j) += (c*phi_i*phi_j - q*phi_j_y*phi_i)
															*fe_face_values.JxW(q_boundary);
									break;
									
								case 3://top
									c = 1;
									g = 3;
									cell_matrix(i,j) += (-p*phi_j_x*phi_i + c*phi_i*phi_j)
															*fe_face_values.JxW(q_boundary);
									break;
								
								default:
									c = 0;
									g = 0;
									cell_matrix(i,j) += (-q*phi_j_y*phi_i - c*phi_i*phi_j)
															*fe_face_values.JxW(q_boundary);
									
							}
						}
						
						int boundary_number = cell->face(face)->boundary_id();
						int c, g;
						
						switch (boundary_number)
						{
							case 1:
								c = 3;
								g = 1;
								cell_rhs(i) += -g*phi_i*fe_face_values.JxW(q_boundary);	
								break;
								
							case 2: 
								c = 2;
								g = 2;
								cell_rhs(i) += g*phi_i*fe_face_values.JxW(q_boundary);	
								break;
								
							case 3:
								c = 1;
								g = 3;
								cell_rhs(i) += g*phi_i*fe_face_values.JxW(q_boundary);	
								break;
							
							default:
									c = 0;
									g = 0;
									cell_rhs(i) += -g*phi_i*fe_face_values.JxW(q_boundary);
						}								
					}
				}
			}
		}	
					
		cell->get_dof_indices(local_dof_indices);
		constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices,
								system_matrix, system_rhs);
	}
	
    timer.stop();
    printf("done (%gs)\n", timer());                           
}

template<int dim>
void Problem5<dim>::solve()
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
void Problem5<dim>::refine_grid(int cycle, const char *s)
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
					0.2, 0.05);
		
		mesh.execute_coarsening_and_refinement();
		
		timer.stop();
		printf("done (%gs)\n",timer());
		
		setup_geometry(cycle);
	}	
}

template<int dim>
void Problem5<dim>::output_results(int cycle, const char *s)
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
void Problem5<dim>::run()
{		
	read_inputs();
	
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
			//calculate_error(cycle, convergence_table);
			
			timer.stop();			
			//convergence_table.add_value("Time", timer());
		}
		
		//print_errors(true, convergence_table);	
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
			
			//calculate_error(cycle, adaptive_error_g_table);
			timer.stop();			
			//adaptive_error_g_table.add_value("Time", timer());
		}
		
		//print_errors(false, adaptive_error_g_table);			
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
			
			//calculate_error(cycle, adaptive_error_h_table);
			timer.stop();			
			//adaptive_error_h_table.add_value("Time", timer());
		}
		
		//print_errors(false, adaptive_error_h_table);		
	}
}

int main ()
{
	const int dim = 2;
	
	Problem5<dim> problem;
	
    problem.run();
}
