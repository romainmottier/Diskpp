
///////////////////////////////////
/////////////////////////////////// 
/////////////////////////////////// 

void set_hdg_stabilization() {
  if(m_hho_di.cell_degree() > m_hho_di.face_degree()) {
    m_hho_stabilization_Q = false;
  }
  else {
    std::cout << "Proceeding with HHO stabilization cell and face degree are equal." << std::endl;
    std::cout << "cell degree = " << m_hho_di.cell_degree() << std::endl;
    std::cout << "face degree = " << m_hho_di.face_degree() << std::endl;
  }
}

///////////////////////////////////
/////////////////////////////////// Indexes des cellules possédant une face à l'interface
///////////////////////////////////

void set_interface_cell_indexes(std::map<size_t,std::pair<size_t,size_t>> & interface_cell_indexes) {
  m_interface_cell_indexes = interface_cell_indexes;
}


///////////////////////////////////
/////////////////////////////////// HHO stabilization : equal order
///////////////////////////////////

void set_hho_stabilization() {
  m_hho_stabilization_Q = true;
}


///////////////////////////////////
/////////////////////////////////// Stabilization : O(1) = true
///////////////////////////////////

void set_scaled_stabilization() {
  m_scaled_stabilization_Q = true;
}
