
void HeterogeneousIHHOFirstOrder(int argc, char **argv);

void HeterogeneousIHHOFirstOrder(int argc, char **argv){
    
  // ##################################################
  // ################################################## Simulation paramaters 
  // ##################################################
  
  using RealType = double;
  simulation_data sim_data = preprocessor::process_args(argc, argv);
  sim_data.print_simulation_data();
  
  // ##################################################
  // ################################################## Mesh generation 
  // ##################################################

  // ################################ Building a cartesian mesh
  timecounter tc;
  tc.tic();
  // Lenght of the domain
  RealType lx = 3.0;
  RealType ly = 3.0;
  // Number of element in each direction
  size_t nx = 2;
  size_t ny = 2;
  typedef disk::mesh<RealType, 2, disk::generic_mesh_storage<RealType, 2>>  mesh_type;
  typedef disk::BoundaryConditions<mesh_type, false> e_boundary_type;
  typedef disk::BoundaryConditions<mesh_type, true> a_boundary_type;
  mesh_type msh;
  cartesian_2d_mesh_builder<RealType> mesh_builder(lx,ly,nx,ny);
  mesh_builder.refine_mesh(sim_data.m_n_divs);
  mesh_builder.set_translation_data(-1.5, -1.5);
  mesh_builder.build_mesh();
  mesh_builder.move_to_mesh_storage(msh);
  tc.toc();
  std::cout << bold << red << "   MESH GENERATION : ";
  std::cout << tc.to_double() << " seconds" << reset << std::endl << std::endl << std::endl;

  // ################################ Polyhedral meshes
  
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
  // ################################################## HHO setting 
  // ##################################################
  
  // Creating HHO approximation spaces and corresponding linear operator
  size_t cell_k_degree = sim_data.m_k_degree;
  if(sim_data.m_hdg_stabilization_Q){
    cell_k_degree++;
  }
  
  disk::hho_degree_info hho_di(cell_k_degree,sim_data.m_k_degree);
  
  // ##################################################
  // ################################################## Analytical (vectorial and scalar) functions 
  // ##################################################               (source term + IC) 
  // ##################################################
  
  auto null_fun = [](const mesh_type::point_type& pt) -> static_vector<RealType, 2> {
    static_vector<RealType, 2> f{0,0};
    return f;
  };
  
  auto null_flux_fun = [](const typename mesh_type::point_type& pt) -> static_matrix<RealType,2,2> {
    double x,y;
    x = pt.x();
    y = pt.y();
    static_matrix<RealType, 2, 2> sigma = static_matrix<RealType,2,2>::Zero(2,2);
    return sigma;
  };
  
  auto null_s_fun = [](const mesh_type::point_type& pt) -> RealType {
    return 0.0;
  };
  
  // Initial condition Gar6more
  auto v_s_fun = [](const mesh_type::point_type& pt) -> static_vector<RealType, 2> {
    RealType x,y,xc,yc,r,wave,vx,vy,c,lp;
    x    = pt.x();
    y    = pt.y();
    xc   = 0.0;
    yc   = 2.0/3.0;
    c    = 10.0;
    lp   = std::sqrt(9.0)/10.0;
    r    = std::sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));
    wave = (c)/(std::exp((1.0/(lp*lp))*r*r*M_PI*M_PI));
    vx   = -wave*(x-xc);
    vy   = -wave*(y-yc);
    static_vector<RealType, 2> v{vx,vy};
    return v;
   };

  auto test = [](const mesh_type::point_type& pt) -> RealType {
    RealType x,y,xc,yc,r,wave;
    x = pt.x();
    y = pt.y();
    xc = 0.0;
    yc = 2.0/3.0;
    r = std::sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));
    wave = 0.1*(-4*std::sqrt(10.0/3.0)*(-1 + 1600.0*r*r))/(std::exp(800*r*r)*std::pow(M_PI,0.25));
    return wave;
  };

  // ##################################################
  // ################################################## Material data 
  // ##################################################
  
  auto elastic_mat_fun = [](const typename mesh_type::point_type& pt)
                                                                 -> elastic_material_data<RealType> {
    double x,y;
    x = pt.x();
    y = pt.y();
    RealType rho, vp, vs;
    rho = 1.0;            // Solid mass density                (ex ground : 2800 kg/m^3)
    vp  = std::sqrt(3.0); // Seismic compressional velocity vp (ex ground : 5000 m/s   ) 
    vs  = 1.0;            // Seismic shear velocity vs         (ex ground : 3000 m/s   )
    elastic_material_data<RealType> material(rho,vp,vs);
    return material;
  };
  
  auto acoustic_mat_fun = [](const typename mesh_type::point_type& pt)
                                                                -> acoustic_material_data<RealType> {
    double x,y;
    x = pt.x();
    y = pt.y();
    RealType rho, vp;
    rho = 1.0;            // Fluid mass density (ex air : 1.292 kg.m^-3)
    vp  = std::sqrt(9.0); // seismic compressional velocity vp (ex air : 340 m/s)
    acoustic_material_data<RealType> material(rho,vp);
    return material;
  };
  
  std::map<size_t,elastic_material_data<RealType>>  e_material;
  std::map<size_t,acoustic_material_data<RealType>> a_material;
  std::set<size_t> elastic_bc_face_indexes, acoustic_bc_face_indexes, interface_face_indexes;
  std::map<size_t,std::pair<size_t,size_t>> interface_cell_pair_indexes;
  
  RealType eps = 1.0e-10;
  for (auto face_it = msh.faces_begin(); face_it != msh.faces_end(); face_it++) {
    const auto face = *face_it;
    mesh_type::point_type bar = barycenter(msh, face);
    auto fc_id = msh.lookup(face);
    if (std::fabs(bar.y()) < eps) {
      interface_face_indexes.insert(fc_id);
      continue;
    }      
  }
  
  for (auto & cell : msh ) {
    auto cell_ind = msh.lookup(cell);
    mesh_type::point_type bar = barycenter(msh, cell);
    
    // Assigning the material properties
    if (bar.y() > 0) {
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
	if (bar.y() > 0) {
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
  for (auto face_it = msh.boundary_faces_begin(); face_it != msh.boundary_faces_end(); face_it++) {
    auto face = *face_it;
    mesh_type::point_type bar = barycenter(msh, face);
    auto fc_id = msh.lookup(face);
    if (bar.y() > 0) {
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
  
  // Detect interface elastic - acoustic
  e_boundary_type e_bnd(msh);
  a_boundary_type a_bnd(msh);
  e_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_elastic_id, null_fun);
  a_bnd.addDirichletBC(disk::DirichletType::DIRICHLET, bc_acoustic_id, null_s_fun);
  
  // ##################################################
  // ################################################## Solving a primal HHO mixed problem 
  // ##################################################
  
  tc.tic();
  auto assembler = elastoacoustic_four_fields_assembler<mesh_type>(msh, hho_di, e_bnd, a_bnd,
								   e_material, a_material);
  assembler.set_interface_cell_indexes(interface_cell_pair_indexes);

  // Stabilization type 
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
  assembler.project_over_cells(msh, x_dof, null_fun, null_flux_fun, null_s_fun, v_s_fun);
  assembler.project_over_faces(msh, x_dof, null_fun, null_s_fun);
  // Tests
  //assembler.project_over_cells(msh, x_dof, null_fun, null_flux_fun, test, v_s_fun);
  //assembler.project_over_faces(msh, x_dof, null_fun, test);
  
  ////////// Post process of the initial data 
  if (sim_data.m_render_silo_files_Q) {
    size_t it = 0;
    std::string silo_file_name = "elasto_acoustic_inhomogeneous_four_fields_";
    postprocessor<mesh_type>::write_silo_four_fields_elastoacoustic(silo_file_name, it, msh,
								    hho_di, x_dof,
								    e_material, a_material, false);
  }
  
  std::ofstream simulation_log("elasto_acoustic_inhomogeneous_four_fields.txt");
  std::ofstream sensor_1_log("s1_elasto_acoustic_four_fields_h.csv");
  std::ofstream sensor_2_log("s2_elasto_acoustic_four_fields_h.csv");
  std::ofstream sensor_3_log("s3_elasto_acoustic_four_fields_h.csv");

  bool e_side_Q = false;
  typename mesh_type::point_type s1_pt(-1.0/3.0, +1.0/3.0);
  typename mesh_type::point_type s2_pt( 0.0, +1.0/3.0);
  typename mesh_type::point_type s3_pt(+1.0/3.0, +1.0/3.0);
  std::pair<typename mesh_type::point_type,size_t> s1_pt_cell = std::make_pair(s1_pt, -1);
  std::pair<typename mesh_type::point_type,size_t> s2_pt_cell = std::make_pair(s2_pt, -1);
  std::pair<typename mesh_type::point_type,size_t> s3_pt_cell = std::make_pair(s3_pt, -1);
  
  std::cout << bold << cyan << "      " << "Sensor at (-1/3,1/3)" << reset << std::endl; 
  postprocessor<mesh_type>::record_velocity_data_elasto_acoustic_four_fields
    (0, s1_pt_cell, msh, hho_di, assembler, x_dof, e_side_Q, sensor_1_log);
  std::cout << bold << cyan << "      " << "Sensor at (0,1/3)" << reset << std::endl;
  postprocessor<mesh_type>::record_velocity_data_elasto_acoustic_four_fields
    (0, s2_pt_cell, msh, hho_di, assembler, x_dof, e_side_Q, sensor_2_log);
  std::cout << bold << cyan << "      " << "Sensor at (1/3,1/3)" << reset << std::endl;
  postprocessor<mesh_type>::record_velocity_data_elasto_acoustic_four_fields
    (0, s3_pt_cell, msh, hho_di, assembler, x_dof, e_side_Q, sensor_3_log);
  
  // if (sim_data.m_report_energy_Q) {
  //     postprocessor<mesh_type>::compute_elastic_energy_three_fields
  //       (msh, hho_di, assembler, t, x_dof, simulation_log);
  // }    
  
  // Solving a first order equation HDG/HHO propagation problem
  Matrix<RealType, Dynamic, Dynamic> a;
  Matrix<RealType, Dynamic, 1> b;
  Matrix<RealType, Dynamic, 1> c;
  
  // DIRK(s) schemes
  int s = 3;
  bool is_sdirk_Q = true;
  
  if (is_sdirk_Q) {
    dirk_butcher_tableau::sdirk_tables(s, a, b, c);
  }else{
    dirk_butcher_tableau::dirk_tables(s, a, b, c);
  }
  
  std::cout << std::endl << std::endl;
  std::cout << bold << red << "   ASSEMBLY 2 : " << std::endl;
  std::cout << bold << cyan << "      First stiffness assembly completed: ";
  tc.tic();
  assembler.assemble(msh, null_fun, null_s_fun);
  tc.toc();
  std::cout << bold << cyan << tc << " seconds" << reset << std::endl;
  assembler.LHS += assembler.COUPLING; 
  dirk_hho_scheme<RealType> dirk_an(assembler.LHS,assembler.RHS,assembler.MASS);
  
  if (sim_data.m_sc_Q) {
    std::vector<std::pair<size_t,size_t>> vec_cell_basis_data(2);
    vec_cell_basis_data[0] = std::make_pair(assembler.get_e_material_data().size(), assembler.get_e_cell_basis_data());
    vec_cell_basis_data[1] = std::make_pair(assembler.get_a_material_data().size(), assembler.get_a_cell_basis_data());
    dirk_an.set_static_condensation_data(vec_cell_basis_data, assembler.get_n_face_dof());
  }
  
  if (is_sdirk_Q) {
    double scale = a(0,0) * dt;
    dirk_an.SetScale(scale);
    std::cout << bold << cyan << "      Matrix decomposed: "; 
    tc.tic();
    dirk_an.ComposeMatrix();
    //        dirk_an.setIterativeSolver();
    dirk_an.DecomposeMatrix();
    tc.toc();
    std::cout << tc << " seconds" << reset << std::endl;
  }
  
  // ##################################################
  // ################################################## Time marching
  // ##################################################
  
  std::cout << std::endl << std::endl;
  
  Matrix<RealType, Dynamic, 1> x_dof_n;
  
  for(size_t it = 1; it <= nt; it++) {
    
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
	assembler.RHS.setZero();
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

    // ##################################################
    // ################################################## Last postprocess
    // ##################################################

    if (sim_data.m_render_silo_files_Q) {
      std::string silo_file_name = "elasto_acoustic_inhomogeneous_four_fields_";
      postprocessor<mesh_type>::write_silo_four_fields_elastoacoustic(silo_file_name, it, msh, hho_di, x_dof, e_material, a_material, false);
    }

    postprocessor<mesh_type>::record_velocity_data_elasto_acoustic_four_fields(it, s1_pt_cell, msh, hho_di, assembler, x_dof, e_side_Q, sensor_1_log);
    postprocessor<mesh_type>::record_velocity_data_elasto_acoustic_four_fields(it, s2_pt_cell, msh, hho_di, assembler, x_dof, e_side_Q, sensor_2_log);
    postprocessor<mesh_type>::record_velocity_data_elasto_acoustic_four_fields(it, s3_pt_cell, msh, hho_di, assembler, x_dof, e_side_Q, sensor_3_log);
    
    //  if (sim_data.m_report_energy_Q) {
    //      postprocessor<mesh_type>::compute_elastic_energy_three_fields(msh, hho_di, assembler, t, x_dof, simulation_log);
    //  }
    
    std::cout << std::endl << std::endl;
    
  }
  
  simulation_log << "Number of equations : " << dirk_an.DirkAnalysis().n_equations() << std::endl;
  simulation_log << "Number of DIRK steps =  " << s << std::endl;
  simulation_log << "Number of time steps =  " << nt << std::endl;
  simulation_log << "Step size =  " << dt << std::endl;
  simulation_log.flush();

}





























