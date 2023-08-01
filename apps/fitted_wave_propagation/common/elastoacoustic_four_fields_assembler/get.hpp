
/////////////////////////////////
///////////////////////////////// Boundary conditions
/////////////////////////////////

e_boundary_type & get_e_bc_conditions() {
  return m_e_bnd;
}
  
a_boundary_type & get_a_bc_conditions() {
  return m_a_bnd;
}

/////////////////////////////////
///////////////////////////////// Material data
/////////////////////////////////

std::map<size_t,elastic_material_data<T>> & get_e_material_data() {
  return m_e_material;
}

std::map<size_t,acoustic_material_data<T>> & get_a_material_data() {
  return m_a_material;
}

/////////////////////////////////
///////////////////////////////// Degree of freedom
/////////////////////////////////

size_t get_a_n_cells_dof() const {
  return m_n_acoustic_cell_dof;
}

size_t get_e_n_cells_dof() const {
  return m_n_elastic_cell_dof;
}

size_t get_n_face_dof() const {
  size_t n_face_dof = m_n_elastic_face_dof + m_n_acoustic_face_dof;
  return n_face_dof;
}

/////////////////////////////////
///////////////////////////////// Basis 
/////////////////////////////////

size_t get_e_cell_basis_data() const {
  size_t n_ten_cbs = disk::sym_matrix_basis_size(m_hho_di.grad_degree(),
						 Mesh::dimension,
						 Mesh::dimension);
  size_t n_vec_cbs = disk::vector_basis_size(m_hho_di.cell_degree(),Mesh::dimension, Mesh::dimension);
  size_t n_cbs = n_ten_cbs + n_vec_cbs;
  return n_cbs;
}
    
size_t get_a_cell_basis_data() const {
  size_t n_vel_scal_cbs = disk::scalar_basis_size(m_hho_di.reconstruction_degree(),Mesh::dimension) -1;
  size_t n_scal_cbs = disk::scalar_basis_size(m_hho_di.cell_degree(), Mesh::dimension);
  size_t n_cbs = n_vel_scal_cbs + n_scal_cbs;
  return n_cbs;
}
