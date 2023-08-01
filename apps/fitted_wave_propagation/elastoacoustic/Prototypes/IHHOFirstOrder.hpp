void IHHOFirstOrder(int argc, char **argv);

void IHHOFirstOrder(int argc, char **argv){
  
  // ##################################################
  // ################################################## Simulation paramaters 
  // ##################################################
  
  using RealType = double;
  simulation_data sim_data = preprocessor::process_args(argc, argv);
  sim_data.print_simulation_data();
  
  // ##################################################
  // ################################################## Mesh generation 
  // ##################################################

  timecounter tc;
  tc.tic();
  
  // Lenght of the domain
  RealType lx = 2.0;
  RealType ly = 1.0;
  
  // Number of element in each direction
  size_t nx = 3;
  size_t ny = 3;
  
  typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
  typedef disk::BoundaryConditions<mesh_type, false> e_boundary_type;
  typedef disk::BoundaryConditions<mesh_type, true> a_boundary_type;
  mesh_type msh;
  
  cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
  mesh_builder.refine_mesh(sim_data.m_n_divs);
  mesh_builder.set_translation_data(-1.0, 0.0);
  mesh_builder.build_mesh();
  mesh_builder.move_to_mesh_storage(msh);
  
  tc.toc();
  std::cout << bold << red << "   MESH GENERATION : ";
  std::cout << tc.to_double() << " seconds" << reset << std::endl << std::endl << std::endl;
  
  // ##################################################
  // ################################################## Time controls 
  // ##################################################
  
  size_t nt = 10;
  for (unsigned int i = 0; i < sim_data.m_nt_divs; i++) {
    nt *= 2;
  }
  
  RealType ti = 0.0;
  RealType tf = 1.0;
  RealType dt = (tf-ti)/nt;

  RealType t = ti;

  // ##################################################
  // ################################################## Manufactured solution 
  // ##################################################

  scal_vec_analytic_functions functions;
  //functions.set_function_type(scal_vec_analytic_functions::EFunctionType::EFunctionNonPolynomial);
  functions.set_function_type(scal_vec_analytic_functions::EFunctionType::EFunctionQuadraticInTime);
  //functions.set_function_type(scal_vec_analytic_functions::EFunctionType::EFunctionQuadraticInSpace);

  // Elastic analytical functions
  auto u_fun    = functions.Evaluate_u(t);
  auto v_fun    = functions.Evaluate_v(t);
  auto a_fun    = functions.Evaluate_a(t);
  auto f_fun    = functions.Evaluate_f(t);
  auto flux_fun = functions.Evaluate_sigma(t);

  // Acoustic analytical functions
  auto s_u_fun    = functions.Evaluate_s_u(t);
  auto s_v_fun    = functions.Evaluate_s_v(t);
  auto s_a_fun    = functions.Evaluate_s_a(t);
  auto s_f_fun    = functions.Evaluate_s_f(t);
  auto s_flux_fun = functions.Evaluate_s_q(t);

  // ##################################################
  // ################################################## HHO setting 
  // ##################################################
  
  // Creating HHO approximation spaces and corresponding linear operator
  size_t cell_k_degree = sim_data.m_k_degree;
  if(sim_data.m_hdg_stabilization_Q){
    cell_k_degree++;
  }

  // Stabilization type 
  disk::hho_degree_info hho_di(cell_k_degree,sim_data.m_k_degree);

  // ##################################################
  // ################################################## Material data 
  // ##################################################
  
  // Classify cells per material data and bc faces
  auto elastic_mat_fun = [](const typename mesh_type::point_type& pt)
                                                                -> elastic_material_data<RealType> {
    double x,y;
    x = pt.x();
    y = pt.y();
    RealType rho, vp, vs;
    rho = 1.0;            // Solid mass density
    vp  = std::sqrt(3.0); // Seismic compressional velocity vp
    vs  = 1.0;            // Seismic shear velocity vs
    elastic_material_data<RealType> material(rho,vp,vs);
    return material;
  };
  
  auto acoustic_mat_fun = [](const typename mesh_type::point_type& pt)
                                                               -> acoustic_material_data<RealType> {
    double x,y;
    x = pt.x();
    y = pt.y();
    RealType rho, vp;
    rho = 1.0; // Fluid mass density
    vp  = 1.0; // Seismic compressional velocity vp
    acoustic_material_data<RealType> material(rho,vp);
    return material;
  };
  
  std::map<size_t,elastic_material_data<RealType>> e_material;
  std::map<size_t,acoustic_material_data<RealType>> a_material;
  std::set<size_t> elastic_bc_face_indexes, acoustic_bc_face_indexes, interface_face_indexes;
  std::map<size_t,std::pair<size_t,size_t>> interface_cell_pair_indexes;
  
  RealType eps = 1.0e-10;
  for (auto face_it = msh.faces_begin(); face_it != msh.faces_end(); face_it++){
    const auto face = *face_it;
    mesh_type::point_type bar = barycenter(msh, face);
    auto fc_id = msh.lookup(face);
    if (std::fabs(bar.x()) < eps) {
      interface_face_indexes.insert(fc_id);
      continue;
    } 
  }
  
  for (auto & cell : msh ) {
    
    auto cell_ind = msh.lookup(cell);
    mesh_type::point_type bar = barycenter(msh, cell);
    
    // Assigning the material properties
    if (bar.x() > 0) {
      acoustic_material_data<RealType> material = acoustic_mat_fun(bar);
      a_material.insert(std::make_pair(cell_ind,material));
    }
    
    else {
      elastic_material_data<RealType> material = elastic_mat_fun(bar);
      e_material.insert(std::make_pair(cell_ind,material));
    }
    
    // Detection of faces on the interfaces
    auto cell_faces = faces(msh,cell);
    for (auto face :cell_faces) {
      auto fc_id = msh.lookup(face);
      bool is_member_Q = interface_face_indexes.find(fc_id) != interface_face_indexes.end();
      if (is_member_Q) {
	if (bar.x() > 0) {
	  interface_cell_pair_indexes[fc_id].second = cell_ind;
	}
	else {
	  interface_cell_pair_indexes[fc_id].first = cell_ind;
	}
      }
    }
  }
  
  size_t bc_elastic_id  = 0;
  size_t bc_acoustic_id = 1;
  for (auto face_it = msh.boundary_faces_begin(); face_it != msh.boundary_faces_end(); face_it++){
    auto face = *face_it;
    mesh_type::point_type bar = barycenter(msh, face);
    auto fc_id = msh.lookup(face);
    if (bar.x() > 0) {
      disk::bnd_info bi{bc_acoustic_id, true};
      msh.backend_storage()->boundary_info.at(fc_id) = bi;
      acoustic_bc_face_indexes.insert(fc_id);
    }
    else {
      disk::bnd_info bi{bc_elastic_id, true};
      msh.backend_storage()->boundary_info.at(fc_id) = bi;
      elastic_bc_face_indexes.insert(fc_id);
    }   
  }
  
  // detect interface elastic - acoustic
  e_boundary_type e_bnd(msh);
  a_boundary_type a_bnd(msh);
  e_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_elastic_id, u_fun);
  a_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_acoustic_id, s_u_fun);
  
  // ##################################################
  // ################################################## Solving a primal HHO mixed problem 
  // ##################################################
  
  tc.tic();
  auto assembler = elastoacoustic_four_fields_assembler<mesh_type>(msh, hho_di, e_bnd, a_bnd,
								   e_material, a_material);

  assembler.set_interface_cell_indexes(interface_cell_pair_indexes);

  if(sim_data.m_hdg_stabilization_Q){
    assembler.set_hdg_stabilization();
  }
  
  if(sim_data.m_scaled_stabilization_Q){
    assembler.set_scaled_stabilization();
  }
  
  tc.toc();
  std::cout << bold << red << "   ASSEMBLY 1 : " << std::endl;
  std::cout << bold << cyan << "      Assembler generation : ";
  std::cout << tc.to_double() << " seconds" << reset << std::endl;
  
  tc.tic();
  assembler.assemble_mass(msh);
  tc.toc();
  std::cout << bold << cyan << "      Mass Assembly : ";
  std::cout << tc.to_double() << " seconds" << reset << std::endl;
  
  tc.tic();
  assembler.assemble_coupling_terms(msh);
  tc.toc();
  std::cout << bold << cyan << "      Coupling assembly : ";
  std::cout << tc.to_double() << " seconds" << reset << std::endl << std::endl << std::endl;
  
  // ##################################################
  // ################################################## Projecting initial data 
  // ##################################################
  
  std::cout << bold << red << "   PROJECTION OF INITIAL DATA : ";

  Matrix<RealType, Dynamic, 1> x_dof;
  assembler.project_over_cells(msh, x_dof, v_fun, flux_fun, s_v_fun, s_flux_fun);
  assembler.project_over_faces(msh, x_dof, v_fun, s_v_fun);
  
  if (sim_data.m_render_silo_files_Q) {
    size_t it = 0;
    std::string silo_file_name = "elasto_acoustic_four_fields_";
    postprocessor<mesh_type>::write_silo_four_fields_elastoacoustic(silo_file_name, it, msh, hho_di,
								    x_dof, e_material, a_material,
								    false);
    std::cout << std::endl;
  }
  
  std::ofstream simulation_log("elasto_acoustic_four_fields.txt");
  
  if (sim_data.m_report_energy_Q) {
    //postprocessor<mesh_type>::compute_elastic_energy_three_fields(msh, hho_di, assembler,
    //                                                              t, x_dof,simulation_log);
  }
  
  // Solving a first order equation HDG/HHO propagation problem
  Matrix<RealType, Dynamic, Dynamic> a;
  Matrix<RealType, Dynamic, 1> b;
  Matrix<RealType, Dynamic, 1> c;
  
  // DIRK(s) schemes
  int s = 3;
  bool is_sdirk_Q = true;
  
  if (is_sdirk_Q) {
    dirk_butcher_tableau::sdirk_tables(s, a, b, c);
  }
 
  else {
    dirk_butcher_tableau::dirk_tables(s, a, b, c);
  }
  
  std::cout << std::endl << std::endl;
  std::cout << bold << red << "   ASSEMBLY 2 : " << std::endl;
  std::cout << bold << cyan << "      First stiffness assembly completed: ";
  tc.tic();
  assembler.assemble(msh, f_fun, s_f_fun);
  tc.toc();
  std::cout << tc << " seconds" << reset << std::endl;
  assembler.LHS += assembler.COUPLING; 
  dirk_hho_scheme<RealType> dirk_an(assembler.LHS,assembler.RHS,assembler.MASS);
  
  if (sim_data.m_sc_Q) {
    std::vector<std::pair<size_t,size_t>> vec_cell_basis_data(2);
    vec_cell_basis_data[0] = std::make_pair(assembler.get_e_material_data().size(),
					    assembler.get_e_cell_basis_data());
    vec_cell_basis_data[1] = std::make_pair(assembler.get_a_material_data().size(),
					    assembler.get_a_cell_basis_data());
    dirk_an.set_static_condensation_data(vec_cell_basis_data, assembler.get_n_face_dof());
  }
  
  if (is_sdirk_Q) {
    double scale = a(0,0) * dt;
    dirk_an.SetScale(scale);
    std::cout << bold << cyan << "      Matrix decomposed: "; 
    tc.tic();
    dirk_an.ComposeMatrix();
    // dirk_an.setIterativeSolver();
    dirk_an.DecomposeMatrix();
    tc.toc();
    std::cout << bold << cyan << tc << " seconds" << reset << std::endl;
  }

  Matrix<RealType, Dynamic, 1> x_dof_n;
  
  // ##################################################
  // ################################################## Time marching
  // ##################################################
  std::cout << std::endl << std::endl;
  
  for(size_t it = 1; it <= nt; it++){
    
    std::cout << bold << red << "   Time step number " << it << reset << std::endl;
    
    RealType tn = dt*(it-1)+ti;
    
    // DIRK step
    tc.tic();
    {
      size_t n_dof = x_dof.rows();
      Matrix<double, Dynamic, Dynamic> k = Matrix<double, Dynamic, Dynamic>::Zero(n_dof, s);
      Matrix<double, Dynamic, 1> Fg, Fg_c,xd;
      xd = Matrix<double, Dynamic, 1>::Zero(n_dof, 1);
      
      RealType t;
      Matrix<double, Dynamic, 1> yn, ki;
      
      x_dof_n = x_dof;
      for (int i = 0; i < s; i++) {
	
	yn = x_dof;
	for (int j = 0; j < s - 1; j++) {
	  yn += a(i,j) * dt * k.block(0, j, n_dof, 1);
	}
	
	t = tn + c(i,0) * dt;
	// Manufactured solution
	auto v_fun      = functions.Evaluate_v(t);
	auto f_fun      = functions.Evaluate_f(t);
	auto s_v_fun    = functions.Evaluate_s_v(t);
	auto s_f_fun    = functions.Evaluate_s_f(t);
        
	assembler.get_e_bc_conditions().updateDirichletFunction(v_fun, 0);
	assembler.get_a_bc_conditions().updateDirichletFunction(s_v_fun, 0);
	assembler.assemble_rhs(msh, f_fun, s_f_fun);
	dirk_an.SetFg(assembler.RHS);
	dirk_an.irk_weight(yn, ki, dt, a(i,i),is_sdirk_Q);
        
	// Accumulated solution
	x_dof_n += dt*b(i,0)*ki;
	k.block(0, i, n_dof, 1) = ki;
      }
    }
    
    tc.toc();
    std::cout << bold << cyan << "      DIRK step completed: " << tc << " seconds" << reset;
    x_dof = x_dof_n;
    t = tn + dt;
    auto v_fun      = functions.Evaluate_v(t);
    auto flux_fun   = functions.Evaluate_sigma(t);
    auto s_v_fun    = functions.Evaluate_s_v(t);
    auto s_flux_fun = functions.Evaluate_s_q(t);
    
    if (sim_data.m_render_silo_files_Q) {
      std::string silo_file_name = "elasto_acoustic_four_fields_";
      postprocessor<mesh_type>::write_silo_four_fields_elastoacoustic(silo_file_name, it, msh,
								      hho_di, x_dof, e_material,
								      a_material, false);
      std::cout << std::endl;
    }
    
    if (sim_data.m_report_energy_Q) {
      //postprocessor<mesh_type>::compute_elastic_energy_three_fields(msh, hho_di, assembler, t,
      //                                                              x_dof, simulation_log);
    }

    std::cout << std::endl << std::endl;
    
    if(it == nt){
      // Computing errors
      postprocessor<mesh_type>::compute_errors_four_fields_elastoacoustic(msh, hho_di, assembler,
									  x_dof, v_fun, flux_fun,
									  s_v_fun, s_flux_fun,
									  simulation_log);
      std::cout << std::endl << std::endl;
    }
    
  }
  
  simulation_log << "Number of equations : " << dirk_an.DirkAnalysis().n_equations() << std::endl;
  simulation_log << "Number of DIRK steps =  " << s << std::endl;
  simulation_log << "Number of time steps =  " << nt << std::endl;
  simulation_log << "Step size =  " << dt << std::endl;
  simulation_log.flush();
}

