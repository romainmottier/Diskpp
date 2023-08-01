
void assemble(const Mesh& msh,
	      std::function<static_vector<T, 2>(const typename Mesh::point_type& )> e_rhs_fun,
	      std::function<T(const typename Mesh::point_type& )> a_rhs_fun) {
  
  auto storage = msh.backend_storage();
  LHS.setZero();
  RHS.setZero();
  
  // elastic block
  for (auto e_chunk : m_e_material) {
    size_t e_cell_ind = m_e_cell_index[e_chunk.first];
    auto& cell = storage->surfaces[e_chunk.first];
    Matrix<T, Dynamic, Dynamic> mixed_operator_loc = e_mixed_operator(e_chunk.second,msh,cell);
    Matrix<T, Dynamic, 1> f_loc = e_mixed_rhs(msh, cell, e_rhs_fun);
    scatter_e_data(e_cell_ind, msh, cell, mixed_operator_loc, f_loc);
  }
  
  // acoustic block
  for (auto a_chunk : m_a_material) {
    size_t a_cell_ind = m_a_cell_index[a_chunk.first];
    auto& cell = storage->surfaces[a_chunk.first];
    Matrix<T, Dynamic, Dynamic> mixed_operator_loc = a_mixed_operator(a_chunk.second, msh, cell);
    Matrix<T, Dynamic, 1> f_loc = a_mixed_rhs(msh, cell, a_rhs_fun);
    scatter_a_data(a_cell_ind, msh, cell, mixed_operator_loc, f_loc);
  }
  finalize();
}

void finalize() {
  LHS.setFromTriplets(m_triplets.begin(),m_triplets.end());
  m_triplets.clear();
}

void assemble_coupling_terms(const Mesh& msh) {
  
  auto storage = msh.backend_storage();
  COUPLING.setZero();
  
  // Coupling blocks
  for (auto chunk : m_interface_cell_indexes) {
    auto& face = storage->edges[chunk.first];
    auto& e_cell = storage->surfaces[chunk.second.first];
    auto& a_cell = storage->surfaces[chunk.second.second];
    Matrix<T, Dynamic, Dynamic> interface_operator_loc = e_interface_operator(msh, face, e_cell, a_cell);
    scatter_ea_interface_data(msh, face, interface_operator_loc);
  }
  finalize_coupling();
}

void finalize_coupling() {
    COUPLING.setFromTriplets(m_c_triplets.begin(),m_c_triplets.end());
    m_c_triplets.clear();
}

void assemble_mass(const Mesh& msh, bool add_vector_mass_Q = true){
  
  auto storage = msh.backend_storage();
  MASS.setZero();
  
  // Elastic block
  for (auto e_chunk : m_e_material) {
    size_t e_cell_ind = m_e_cell_index[e_chunk.first];
    auto& cell = storage->surfaces[e_chunk.first];
    elastic_material_data<T> material = e_chunk.second;
    Matrix<T, Dynamic, Dynamic> mass_matrix = e_mass_operator(material, msh, cell, add_vector_mass_Q);
    scatter_e_mass_data(e_cell_ind,msh, cell, mass_matrix);
  }
  
  // Acoustic block
  for (auto a_chunk : m_a_material) {
    size_t a_cell_ind = m_a_cell_index[a_chunk.first];
    auto& cell = storage->surfaces[a_chunk.first];
    acoustic_material_data<T> material = a_chunk.second;
    Matrix<T, Dynamic, Dynamic> mass_matrix = a_mass_operator(material, msh, cell);
    scatter_a_mass_data(a_cell_ind, msh, cell, mass_matrix);
  }
  
  finalize_mass();
  
}

void finalize_mass() {
  MASS.setFromTriplets(m_mass_triplets.begin(),m_mass_triplets.end());
  m_mass_triplets.clear();
}

void assemble_rhs(const Mesh& msh,
		  std::function<static_vector<T, 2>(const typename Mesh::point_type& )> e_rhs_fun,
		  std::function<T(const typename Mesh::point_type& )> a_rhs_fun) {
  
  RHS.setZero();
  
#ifdef HAVE_INTEL_TBB2
  
  size_t n_cells = msh.cells_size();
  tbb::parallel_for(
		 size_t(0),
		 size_t(n_cells),
		 size_t(1),
		 [this,&msh,&rhs_fun] (size_t & cell_ind){
		   auto& cell = msh.backend_storage()->surfaces[cell_ind];
		   auto cell_basis   = make_vector_monomial_basis(msh, cell, m_hho_di.cell_degree());
		   Matrix<T, Dynamic, 1> f_loc = make_rhs(msh, cell, cell_basis, rhs_fun);
		   this->scatter_rhs_data(msh, cell, f_loc);
		 });
#else
  
  auto storage = msh.backend_storage();
  
  for (auto e_chunk : m_e_material) {
    size_t e_cell_ind = m_e_cell_index[e_chunk.first];
    auto& cell = storage->surfaces[e_chunk.first];
    Matrix<T, Dynamic, 1> f_loc = e_mixed_rhs(msh, cell, e_rhs_fun);
    scatter_e_rhs_data(e_cell_ind, msh, cell, f_loc);
  }
  
  for (auto a_chunk : m_a_material) {
    size_t a_cell_ind = m_a_cell_index[a_chunk.first];
    auto& cell = storage->surfaces[a_chunk.first];
    Matrix<T, Dynamic, 1> f_loc = a_mixed_rhs(msh, cell, a_rhs_fun);
    scatter_a_rhs_data(a_cell_ind, msh, cell, f_loc);
  }
  
#endif
  
  apply_bc(msh);

}

