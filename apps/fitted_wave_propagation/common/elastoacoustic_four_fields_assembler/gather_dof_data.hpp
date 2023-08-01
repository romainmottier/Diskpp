
Matrix<T, Dynamic, 1> gather_e_dof_data(size_t cell_ind,
					const Mesh& msh,
					const typename Mesh::cell_type& cl,
					const Matrix<T, Dynamic, 1>& x_glob) const {

  // Cell index
  auto e_cell_ind = m_e_cell_index.find(cell_ind)->second;

  // Face index of the selected cell 
  auto num_faces = howmany_faces(msh, cl);

  // elastic basis 
  size_t n_cbs = get_e_cell_basis_data();
  size_t n_fbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);

  
  Matrix<T, Dynamic, 1> x_el(n_cbs + num_faces * n_fbs );
  x_el.block(0, 0, n_cbs, 1) = x_glob.block(e_cell_ind * n_cbs, 0, n_cbs, 1);

  auto fcs = faces(msh, cl);

  for (size_t i = 0; i < fcs.size(); i++) {

    auto fc  = fcs[i];
    auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
    if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
    const auto face_id = eid.second;
    
    if (m_e_bnd.is_dirichlet_face( face_id)) {
      auto fb = make_vector_monomial_basis(msh, fc, m_hho_di.face_degree());
      auto dirichlet_fun  = m_e_bnd.dirichlet_boundary_func(face_id);
      Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, fb);
      Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, fb, dirichlet_fun);
      x_el.block(n_cbs + i * n_fbs, 0, n_fbs, 1) = mass.llt().solve(rhs);
    }
    
    else {
      auto fc_id = msh.lookup(fc);
      auto face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof
	                                          + m_e_compress_indexes.at(fc_id)*n_fbs;
      x_el.block(n_cbs + i*n_fbs, 0, n_fbs, 1) = x_glob.block(face_LHS_offset, 0, n_fbs, 1);
    }
    
  }
  
  return x_el;
  
}


Matrix<T, Dynamic, 1> gather_a_dof_data(size_t cell_ind,
					const Mesh& msh,
					const typename Mesh::cell_type& cl,
					const Matrix<T, Dynamic, 1>& x_glob) const {
  
  auto a_cell_ind = m_a_cell_index.find(cell_ind)->second;
  auto num_faces  = howmany_faces(msh, cl);
  size_t n_cbs    = get_a_cell_basis_data();
  size_t n_fbs    = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
  
  Matrix<T, Dynamic, 1> x_local(n_cbs + num_faces * n_fbs );
  x_local.block(0,0,n_cbs,1) = x_glob.block(a_cell_ind*n_cbs + m_n_elastic_cell_dof, 0, n_cbs, 1);
  auto fcs = faces(msh, cl);
  
  for (size_t i = 0; i < fcs.size(); i++) {
    
    auto fc = fcs[i];
    auto eid = find_element_id(msh.faces_begin(), msh.faces_end(), fc);
    if (!eid.first) throw std::invalid_argument("This is a bug: face not found");
    const auto face_id = eid.second;
    
    if (m_a_bnd.is_dirichlet_face( face_id)) {
      auto fb = disk::make_scalar_monomial_basis(msh, fc, m_hho_di.face_degree());
      Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, fb, m_hho_di.face_degree());
      auto velocity = m_a_bnd.dirichlet_boundary_func(face_id);
      Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, fb, velocity, m_hho_di.face_degree());
      x_local.block(n_cbs + i * n_fbs, 0, n_fbs, 1) = mass.llt().solve(rhs);
    }
    
    else {
      auto fc_id = msh.lookup(fc);
      auto face_LHS_offset = m_n_elastic_cell_dof + m_n_acoustic_cell_dof + m_n_elastic_face_dof
	+ m_a_compress_indexes.at(fc_id)*n_fbs;
      x_local.block(n_cbs + i*n_fbs, 0, n_fbs, 1) = x_glob.block(face_LHS_offset, 0, n_fbs, 1);
    }
    
  }
  
  return x_local;
  
}
