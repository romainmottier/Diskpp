
void scatter_e_data(size_t e_cell_ind,
		    const Mesh& msh,
		    const typename Mesh::cell_type& cl,
		    const Matrix<T,Dynamic,Dynamic>& lhs,
		    const Matrix<T,Dynamic,1>& rhs) {

  auto fcs = faces(msh, cl);
  size_t n_cbs = get_e_cell_basis_data();
  size_t n_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
  std::vector<assembly_index> asm_map;
  asm_map.reserve(n_cbs + n_fbs*fcs.size());
  
  auto cell_LHS_offset = e_cell_ind * n_cbs;
  
  for (size_t i = 0; i < n_cbs; i++)
    asm_map.push_back(assembly_index(cell_LHS_offset+i,true));
  
  for (size_t face_i = 0; face_i < fcs.size(); face_i++) {
    auto fc = fcs[face_i];
    auto fc_id = msh.lookup(fc);
    auto face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof
                                                + m_e_compress_indexes.at(fc_id)*n_fbs;
    bool dirichlet = m_e_bnd.is_dirichlet_face(fc_id);	
    for (size_t i = 0; i < n_fbs; i++)
      asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );
  }
  
  assert( asm_map.size() == lhs.rows() && asm_map.size() == lhs.cols() );
  
  for (size_t i = 0; i < lhs.rows(); i++) {
    if (!asm_map[i].assemble())
      continue;
    for (size_t j = 0; j < lhs.cols(); j++) {
      if (asm_map[j].assemble())
	m_triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], lhs(i,j)) );
    }
  }
  
  for (size_t i = 0; i < rhs.rows(); i++) {
    if (!asm_map[i].assemble())
      continue;
    RHS(int(asm_map[i])) += rhs(i);
  }
  
}


void scatter_a_data(size_t a_cell_ind,
		    const Mesh& msh,
		    const typename Mesh::cell_type& cl,
		    const Matrix<T, Dynamic, Dynamic>& lhs,
		    const Matrix<T, Dynamic, 1>& rhs) {

  auto fcs = faces(msh, cl);
  size_t n_cbs = get_a_cell_basis_data();
  size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
  std::vector<assembly_index> asm_map;
  asm_map.reserve(n_cbs + n_fbs*fcs.size());
  
  auto cell_LHS_offset    = a_cell_ind * n_cbs + m_n_elastic_cell_dof;
  
  for (size_t i = 0; i < n_cbs; i++)
    asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );
  
  for (size_t face_i = 0; face_i < fcs.size(); face_i++) {
    auto fc = fcs[face_i];
    auto fc_id = msh.lookup(fc);
    auto face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_n_elastic_face_dof
                                                + m_a_compress_indexes.at(fc_id)*n_fbs;
    bool dirichlet = m_a_bnd.is_dirichlet_face(fc_id);
    for (size_t i = 0; i < n_fbs; i++)
      asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );
  }

  assert( asm_map.size() == lhs.rows() && asm_map.size() == lhs.cols() );
  
  for (size_t i = 0; i < lhs.rows(); i++) {
    if (!asm_map[i].assemble())
      continue;
    for (size_t j = 0; j < lhs.cols(); j++) {
      if ( asm_map[j].assemble() )
	m_triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], lhs(i,j)) );
    }
  }
  
  for (size_t i = 0; i < rhs.rows(); i++) {
    if (!asm_map[i].assemble())
      continue;
    RHS(int(asm_map[i])) += rhs(i);
  }
  
}
    

void scatter_e_mass_data(size_t e_cell_ind,
			 const Mesh& msh,
			 const typename Mesh::cell_type& cl,
			 const Matrix<T, Dynamic, Dynamic>& mass_matrix) {
        
  size_t n_cbs = get_e_cell_basis_data();
  std::vector<assembly_index> asm_map;
  asm_map.reserve(n_cbs);
  
  auto cell_LHS_offset    = e_cell_ind * n_cbs;
  
  for (size_t i = 0; i < n_cbs; i++)
    asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );
  
  assert( asm_map.size() == mass_matrix.rows() && asm_map.size() == mass_matrix.cols() );
  
  for (size_t i = 0; i < mass_matrix.rows(); i++) {
    if (!asm_map[i].assemble())
      continue;
    for (size_t j = 0; j < mass_matrix.cols(); j++) {
      if ( asm_map[j].assemble() )
	m_mass_triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], mass_matrix(i,j)) );
    }
  }
  
}


void scatter_a_mass_data(size_t a_cell_ind,
			 const Mesh& msh,
			 const typename Mesh::cell_type& cl,
			 const Matrix<T, Dynamic, Dynamic>& mass_matrix) {
  
  size_t n_cbs = get_a_cell_basis_data();
  std::vector<assembly_index> asm_map;
  asm_map.reserve(n_cbs);
  
  auto cell_LHS_offset    = a_cell_ind * n_cbs + m_n_elastic_cell_dof;
  
  for (size_t i = 0; i < n_cbs; i++)
    asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );
  
  assert( asm_map.size() == mass_matrix.rows() && asm_map.size() == mass_matrix.cols() );
  
  for (size_t i = 0; i < mass_matrix.rows(); i++) {
    if (!asm_map[i].assemble())
      continue;
    for (size_t j = 0; j < mass_matrix.cols(); j++) {
      if ( asm_map[j].assemble() )
	m_mass_triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], mass_matrix(i,j)) );
    }
  }
  
}


void scatter_ea_interface_data(const Mesh& msh,
			       const typename Mesh::face_type& face,
			       const Matrix<T, Dynamic, Dynamic>& interface_matrix) {

  auto vfbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
  auto sfbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
  
  std::vector<assembly_index> asm_map_e, asm_map_a;
  asm_map_e.reserve(vfbs);
  asm_map_a.reserve(sfbs);
  
  auto fc_id = msh.lookup(face);
  
  auto e_face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof
                                                + m_e_compress_indexes.at(fc_id)*vfbs;

  bool e_dirichlet = m_e_bnd.is_dirichlet_face(fc_id);
  for (size_t i = 0; i < vfbs; i++) {
    asm_map_e.push_back( assembly_index(e_face_LHS_offset+i, !e_dirichlet) );
  }
  
  auto a_face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_n_elastic_face_dof
                                                + m_a_compress_indexes.at(fc_id)*sfbs;
  
  bool a_dirichlet = m_a_bnd.is_dirichlet_face(fc_id);
  
  for (size_t i = 0; i < sfbs; i++) {
    asm_map_a.push_back( assembly_index(a_face_LHS_offset+i, !a_dirichlet) );
  }
  
  assert( asm_map_e.size() == interface_matrix.rows() && asm_map_a.size() == interface_matrix.cols() );
    
  for (size_t i = 0; i < interface_matrix.rows(); i++) {
    for (size_t j = 0; j < interface_matrix.cols(); j++) {
      m_c_triplets.push_back( Triplet<T>(asm_map_e[i], asm_map_a[j],   interface_matrix(i,j)) );
      m_c_triplets.push_back( Triplet<T>(asm_map_a[j], asm_map_e[i], - interface_matrix(i,j)) );
    }
  }
  
}


void scatter_e_rhs_data(size_t e_cell_ind,
			const Mesh& msh,
			const typename Mesh::cell_type& cl,
			const Matrix<T, Dynamic, 1>& rhs) {
  
  size_t n_cbs = get_e_cell_basis_data();
  std::vector<assembly_index> asm_map;
  asm_map.reserve(n_cbs);
  
  auto cell_LHS_offset = e_cell_ind * n_cbs;
  
  for (size_t i = 0; i < n_cbs; i++)
    asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );
  
  assert( asm_map.size() == rhs.rows() );
  
  for (size_t i = 0; i < rhs.rows(); i++) {
    if (!asm_map[i].assemble())
      continue;
    RHS(int(asm_map[i])) += rhs(i);
  }
  
}


void scatter_a_rhs_data(size_t a_cell_ind,
			const Mesh& msh,
			const typename Mesh::cell_type& cl,
			const Matrix<T, Dynamic, 1>& rhs) {

  size_t n_cbs = get_a_cell_basis_data();
  std::vector<assembly_index> asm_map;
  asm_map.reserve(n_cbs);
  
  auto cell_LHS_offset = a_cell_ind * n_cbs + m_n_elastic_cell_dof;
  
  for (size_t i = 0; i < n_cbs; i++)
    asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );
  
  assert( asm_map.size() == rhs.rows());
  
  for (size_t i = 0; i < rhs.rows(); i++) {
    if (!asm_map[i].assemble())
      continue;
    RHS(int(asm_map[i])) += rhs(i);
  }
  
}


void scatter_e_bc_data(size_t e_cell_ind,
		       const Mesh& msh,
		       const typename Mesh::cell_type& cl,
		       const Matrix<T, Dynamic, Dynamic>& lhs) {
    
  auto fcs = faces(msh, cl);
  size_t n_cbs = get_e_cell_basis_data();
  size_t n_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
  std::vector<assembly_index> asm_map;
  asm_map.reserve(n_cbs + n_fbs*fcs.size());
  
  auto cell_LHS_offset  = e_cell_ind * n_cbs;
  
  for (size_t i = 0; i < n_cbs; i++)
    asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );
  
  Matrix<T, Dynamic, 1> dirichlet_data = Matrix<T, Dynamic, 1>::Zero(n_cbs + fcs.size()*n_fbs);
  for (size_t face_i = 0; face_i < fcs.size(); face_i++) {
    auto fc = fcs[face_i];
    auto fc_id = msh.lookup(fc);
    auto face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof
                                                + m_e_compress_indexes.at(fc_id)*n_fbs;
    bool dirichlet = m_e_bnd.is_dirichlet_face(fc_id);
    
    for (size_t i = 0; i < n_fbs; i++)
      asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );
    
    if (dirichlet) {
      auto fb = make_vector_monomial_basis(msh, fc, m_hho_di.face_degree());
      auto dirichlet_fun  = m_e_bnd.dirichlet_boundary_func(fc_id);
      
      Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, fb);
      Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, fb, dirichlet_fun);
      dirichlet_data.block(n_cbs + face_i*n_fbs, 0, n_fbs, 1) = mass.llt().solve(rhs);
    }  
  }
  
  assert( asm_map.size() == lhs.rows() && asm_map.size() == lhs.cols() );
  
  for (size_t i = 0; i < lhs.rows(); i++) {
    if (!asm_map[i].assemble())
      continue;
    
    for (size_t j = 0; j < lhs.cols(); j++) {
      if ( !asm_map[j].assemble() )
	RHS(int(asm_map[i])) -= lhs(i,j) * dirichlet_data(j);
    }
  }
  
}


void scatter_a_bc_data(size_t a_cell_ind,
		       const Mesh& msh,
		       const typename Mesh::cell_type& cl,
		       const Matrix<T, Dynamic, Dynamic>& lhs) {
  
  auto fcs = faces(msh, cl);
  size_t n_cbs = get_a_cell_basis_data();
  size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
  std::vector<assembly_index> asm_map;
  asm_map.reserve(n_cbs + n_fbs*fcs.size());
  
  auto cell_LHS_offset = a_cell_ind * n_cbs + m_n_elastic_cell_dof;
  
  for (size_t i = 0; i < n_cbs; i++)
    asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );
  
  Matrix<T, Dynamic, 1> dirichlet_data = Matrix<T, Dynamic, 1>::Zero(n_cbs + fcs.size()*n_fbs);
  for (size_t face_i = 0; face_i < fcs.size(); face_i++) {
    auto fc = fcs[face_i];
    auto fc_id = msh.lookup(fc);
    auto face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_n_elastic_face_dof
                                                + m_a_compress_indexes.at(fc_id)*n_fbs;
            
    bool dirichlet = m_a_bnd.is_dirichlet_face(fc_id);
    
    for (size_t i = 0; i < n_fbs; i++)
      asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );
    
    if (dirichlet) {
      auto fb = make_scalar_monomial_basis(msh, fc, m_hho_di.face_degree());
      auto dirichlet_fun  = m_a_bnd.dirichlet_boundary_func(fc_id);
      Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, fb);
      Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, fb, dirichlet_fun);
      dirichlet_data.block(n_cbs + face_i*n_fbs, 0, n_fbs, 1) = mass.llt().solve(rhs);
    }  
  }

  assert( asm_map.size() == lhs.rows() && asm_map.size() == lhs.cols() );
  
  for (size_t i = 0; i < lhs.rows(); i++) {
    if (!asm_map[i].assemble())
      continue;
    
    for (size_t j = 0; j < lhs.cols(); j++) {
      if ( !asm_map[j].assemble() )
	RHS(int(asm_map[i])) -= lhs(i,j) * dirichlet_data(j);
    }
  }  
}


void scatter_e_cell_dof_data(size_t e_cell_ind,
			     const Mesh& msh,
			     const typename Mesh::cell_type& cell,
			     Matrix<T, Dynamic, 1>& x_glob,
			     Matrix<T, Dynamic, 1> x_proj_dof) const {

  size_t n_cbs = get_e_cell_basis_data();
  auto cell_ofs = e_cell_ind * n_cbs;
  x_glob.block(cell_ofs, 0, n_cbs, 1) = x_proj_dof;
}
    

void scatter_a_cell_dof_data(size_t a_cell_ind,
			     const Mesh& msh,
			     const typename Mesh::cell_type& cell,
			     Matrix<T, Dynamic, 1>& x_glob,
			     Matrix<T, Dynamic, 1> x_proj_dof) const {

  size_t n_cbs = get_a_cell_basis_data();
  auto cell_ofs = a_cell_ind * n_cbs + m_n_elastic_cell_dof;
  x_glob.block(cell_ofs, 0, n_cbs, 1) = x_proj_dof;
}
  
void scatter_e_face_dof_data(const Mesh& msh,
			     const typename Mesh::face_type& face,
			     Matrix<T, Dynamic, 1>& x_glob,
			     Matrix<T, Dynamic, 1> x_proj_dof) const {

  size_t n_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
  auto fc_id = msh.lookup(face);

  auto glob_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof
                                          + m_e_compress_indexes.at(fc_id)*n_fbs;
  x_glob.block(glob_offset, 0, n_fbs, 1) = x_proj_dof;
}


void scatter_a_face_dof_data(const Mesh& msh,
			     const typename Mesh::face_type& face,
			     Matrix<T, Dynamic, 1>& x_glob,
			     Matrix<T, Dynamic, 1> x_proj_dof) const {

  size_t n_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
  auto fc_id = msh.lookup(face);
  auto glob_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_n_elastic_face_dof
                                          + m_a_compress_indexes.at(fc_id)*n_fbs;
  x_glob.block(glob_offset, 0, n_fbs, 1) = x_proj_dof;  

}
