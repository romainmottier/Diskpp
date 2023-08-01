
void classify_cells(const Mesh& msh){
    
  m_e_elements_with_bc_eges.clear();
  for (auto& cell : msh) {
    auto cell_ind = msh.lookup(cell);
    auto face_list = faces(msh, cell);
    for (size_t face_i = 0; face_i < face_list.size(); face_i++) {
      auto fc = face_list[face_i];
      auto fc_id = msh.lookup(fc);
      bool is_dirichlet_Q = m_e_bnd.is_dirichlet_face(fc_id);
      if (is_dirichlet_Q) {
	m_e_elements_with_bc_eges.push_back(cell_ind);
	break;
      }
    }
  }
  
  m_a_elements_with_bc_eges.clear();
  for (auto& cell : msh) {
    typename Mesh::point_type bar = barycenter(msh, cell);
    if (bar.x() < 0) {
      continue;
    }
    
    auto cell_ind = msh.lookup(cell);
    auto face_list = faces(msh, cell);
    for (size_t face_i = 0; face_i < face_list.size(); face_i++) {
      auto fc = face_list[face_i];
      auto fc_id = msh.lookup(fc);
      bool is_dirichlet_Q = m_a_bnd.is_dirichlet_face(fc_id);
      if (is_dirichlet_Q) {
	m_a_elements_with_bc_eges.push_back(cell_ind);
	break;
      }
    }
  }
}


void build_cells_maps(){
  
  // elastic data
  size_t e_cell_ind = 0;
  for (auto chunk : m_e_material) {
    m_e_cell_index.insert(std::make_pair(chunk.first,e_cell_ind));
    e_cell_ind++;
  }
  
  // acoustic data
  size_t a_cell_ind = 0;
  for (auto chunk : m_a_material) {
    m_a_cell_index.insert(std::make_pair(chunk.first,a_cell_ind));
    a_cell_ind++;
  }
}
