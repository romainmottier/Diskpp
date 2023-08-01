
void apply_bc(const Mesh& msh){
  
#ifdef HAVE_INTEL_TBB2
  
  size_t n_cells = m_elements_with_bc_eges.size();
  tbb::parallel_for(size_t(0), size_t(n_cells), size_t(1), [this,&msh] (size_t & i) {
    size_t cell_ind = m_elements_with_bc_eges[i];
    auto& cell = msh.backend_storage()->surfaces[cell_ind];
    Matrix<T, Dynamic, Dynamic> laplacian_operator_loc = laplacian_operator(cell_ind, msh, cell);
    scatter_bc_data(msh, cell, laplacian_operator_loc);
  }
    );

#else
  
  auto storage = msh.backend_storage();
  for (auto& cell_ind : m_e_elements_with_bc_eges) {
    auto& cell = storage->surfaces[cell_ind];
    size_t e_cell_ind = m_e_cell_index[cell_ind];
    elastic_material_data<T> e_mat = m_e_material.find(cell_ind)->second;
    Matrix<T, Dynamic, Dynamic> mixed_operator_loc = e_mixed_operator(e_mat,msh,cell);
    scatter_e_bc_data(e_cell_ind, msh, cell, mixed_operator_loc);
  }
  
  for (auto& cell_ind : m_a_elements_with_bc_eges) {
    auto& cell = storage->surfaces[cell_ind];
    size_t a_cell_ind = m_a_cell_index[cell_ind];
    acoustic_material_data<T> a_mat = m_a_material.find(cell_ind)->second;
    Matrix<T, Dynamic, Dynamic> mixed_operator_loc = a_mixed_operator(a_mat, msh, cell);
    scatter_a_bc_data(a_cell_ind, msh, cell, mixed_operator_loc);
  }
  
#endif
  
}


Matrix<T, Dynamic, 1> e_neumman_bc_operator(
				      const Mesh& msh,
				      const typename Mesh::face_type& face,
				      const typename Mesh::cell_type& e_cell,
				      const typename Mesh::cell_type& a_cell,
				      std::function<T(const typename Mesh::point_type& )> a_vel_fun) {

  Matrix<T, Dynamic, Dynamic> neumann_operator;
  auto facdeg = m_hho_di.face_degree();
  auto vfbs = disk::vector_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1, Mesh::dimension);
  
  neumann_operator = Matrix<T, Dynamic, Dynamic>::Zero(vfbs, 1);
  
  auto vfb = make_vector_monomial_basis(msh, face, facdeg);
  auto sfb = make_scalar_monomial_basis(msh, face, facdeg);
  const auto qps = integrate(msh, face, facdeg);
  const auto n = disk::normal(msh, e_cell, face);
  
  for (auto& qp : qps) {
    const auto v_f_phi = vfb.eval_functions(qp.point());
    const auto s_f_phi = sfb.eval_functions(qp.point());
    
    assert(v_f_phi.rows() == vfbs);
    const auto n_dot_v_f_phi = disk::priv::inner_product(v_f_phi,
							 disk::priv::inner_product(qp.weight(), n));
    const auto result = disk::priv::outer_product(n_dot_v_f_phi, a_vel_fun(qp.point()));
    neumann_operator += result;
    }
  
    return neumann_operator;
  }
  

Matrix<T, Dynamic, Dynamic> a_neumman_bc_operator(
		    const Mesh& msh, const typename Mesh::face_type& face,
		    const typename Mesh::cell_type& e_cell,
		    const typename Mesh::cell_type& a_cell,
		    std::function<static_vector<T,2>(const typename Mesh::point_type& )> e_vel_fun) {

  Matrix<T, Dynamic, Dynamic> neumann_operator;
  auto facdeg = m_hho_di.face_degree();
  auto sfbs = disk::scalar_basis_size(facdeg, Mesh::dimension - 1);
  
  neumann_operator = Matrix<T, Dynamic, Dynamic>::Zero(sfbs, 1);
  
  auto vfb = make_vector_monomial_basis(msh, face, facdeg);
  auto sfb = make_scalar_monomial_basis(msh, face, facdeg);
  const auto qps = integrate(msh, face, facdeg);
  const auto n = disk::normal(msh, a_cell, face);
  
  for (auto& qp : qps) {
    const auto v_f_phi = vfb.eval_functions(qp.point());
    const auto s_f_phi = sfb.eval_functions(qp.point());
    
    assert(s_f_phi.rows() == sfbs);
    const auto n_dot_v_f = disk::priv::inner_product(e_vel_fun(qp.point()),
						     disk::priv::inner_product(qp.weight(), n));
    
    const auto result = disk::priv::inner_product(n_dot_v_f, s_f_phi);
    neumann_operator += result;
  }

  return neumann_operator;

}
