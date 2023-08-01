
#pragma once
#ifndef elastoacoustic_four_fields_assembler_hpp
#define elastoacoustic_four_fields_assembler_hpp

#include "bases/bases.hpp"
#include "methods/hho"
#include "../common/assembly_index.hpp"
#include "../common/acoustic_material_data.hpp"
#include "../common/elastic_material_data.hpp"
#include <map>

#ifdef HAVE_INTEL_TBB
#include <tbb/parallel_for.h>
#endif

template<typename Mesh>
class elastoacoustic_four_fields_assembler{
    
  typedef disk::BoundaryConditions<Mesh, false>    e_boundary_type;
  typedef disk::BoundaryConditions<Mesh, true>     a_boundary_type;
  using T = typename Mesh::coordinate_type;
  
  std::vector<size_t>                 m_e_compress_indexes;
  std::vector<size_t>                 m_e_expand_indexes;
  
  std::vector<size_t>                 m_a_compress_indexes;
  std::vector<size_t>                 m_a_expand_indexes;
  
  disk::hho_degree_info               m_hho_di;
  e_boundary_type                     m_e_bnd;
  a_boundary_type                     m_a_bnd;
  std::vector< Triplet<T> >           m_triplets;
  std::vector< Triplet<T> >           m_c_triplets;
  std::vector< Triplet<T> >           m_mass_triplets;
  std::map<size_t,elastic_material_data<T>> m_e_material;
  std::map<size_t,acoustic_material_data<T>> m_a_material;
  std::map<size_t,size_t> m_e_cell_index;
  std::map<size_t,size_t> m_a_cell_index;
  std::vector< size_t >               m_e_elements_with_bc_eges;
  std::vector< size_t >               m_a_elements_with_bc_eges;
  std::map<size_t,std::pair<size_t,size_t>>   m_interface_cell_indexes;
  
  size_t      m_n_edges;
  size_t      m_n_essential_edges;
  size_t      m_n_elastic_cell_dof;
  size_t      m_n_acoustic_cell_dof;
  size_t      m_n_elastic_face_dof;
  size_t      m_n_acoustic_face_dof;
  bool        m_hho_stabilization_Q;
  bool        m_scaled_stabilization_Q;
  
public:

  SparseMatrix<T>         LHS;
  Matrix<T, Dynamic, 1>   RHS;
  SparseMatrix<T>         MASS;
  SparseMatrix<T>         COUPLING;
  
  elastoacoustic_four_fields_assembler(const Mesh& msh,
				       const disk::hho_degree_info& hho_di,
				       const e_boundary_type& e_bnd,
				       const a_boundary_type& a_bnd,
				       std::map<size_t,elastic_material_data<T>> & e_material,
				       std::map<size_t,acoustic_material_data<T>> & a_material) :
                                       m_hho_di(hho_di),
				       m_e_bnd(e_bnd),
				       m_a_bnd(a_bnd),
				       m_e_material(e_material),
				       m_a_material(a_material),
				       m_hho_stabilization_Q(true),
				       m_scaled_stabilization_Q(false) {
     
    auto storage = msh.backend_storage();

    ////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////// Detection of elastic edges 
    ////////////////////////////////////////////////////////////

    auto is_e_dirichlet = [&](const typename Mesh::face& fc) -> bool {
      auto fc_id = msh.lookup(fc);
      return e_bnd.is_dirichlet_face(fc_id);
    };
    
    size_t n_e_essential_edges = std::count_if(msh.faces_begin(), msh.faces_end(), is_e_dirichlet);

    std::set<size_t> e_egdes;
    for (auto &chunk : m_e_material) {
      size_t cell_i = chunk.first;
      auto& cell = storage->surfaces[cell_i];
      auto cell_faces = faces(msh,cell);
      for (auto &face : cell_faces) {
	if (!is_e_dirichlet(face)) {
	  auto fc_id = msh.lookup(face);
	  e_egdes.insert(fc_id);
	}
      }
    }
    size_t n_e_edges = e_egdes.size();
    
    ////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////// Detection of acoustic edges 
    ////////////////////////////////////////////////////////////

    auto is_a_dirichlet = [&](const typename Mesh::face& fc) -> bool {
      auto fc_id = msh.lookup(fc);
      return a_bnd.is_dirichlet_face(fc_id);
    };

    size_t n_a_essential_edges = std::count_if(msh.faces_begin(), msh.faces_end(), is_a_dirichlet);    
       
    std::set<size_t> a_egdes;
    for (auto &chunk : m_a_material) {
      size_t cell_i = chunk.first;
      auto& cell = storage->surfaces[cell_i];
      auto cell_faces = faces(msh,cell);
      for (auto &face : cell_faces) {
	if (!is_a_dirichlet(face)) {
	  auto fc_id = msh.lookup(face);
	  a_egdes.insert(fc_id);
	}
      }
    }
    size_t n_a_edges = a_egdes.size();

    ////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////// 
    ////////////////////////////////////////////////////////////
    
    m_n_edges = msh.faces_size();
    m_n_essential_edges = n_e_essential_edges + n_a_essential_edges;

    m_e_compress_indexes.resize( m_n_edges );
    m_e_expand_indexes.resize( m_n_edges - m_n_essential_edges );
    
    m_a_compress_indexes.resize( m_n_edges );
    m_a_expand_indexes.resize( m_n_edges - m_n_essential_edges );    
        
    size_t e_compressed_offset = 0;
    for (auto face_id : e_egdes) {
      m_e_compress_indexes.at(face_id) = e_compressed_offset;
      m_e_expand_indexes.at(e_compressed_offset) = face_id;
      e_compressed_offset++;
    }
        
    size_t a_compressed_offset = 0;
    for (auto face_id : a_egdes) {
      m_a_compress_indexes.at(face_id) = a_compressed_offset;
      m_a_expand_indexes.at(a_compressed_offset) = face_id;
      a_compressed_offset++;
    }

    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////  Elastic and Acoustic basis 
    ////////////////////////////////////////////////////////////

    size_t n_cbs = get_e_cell_basis_data();
    size_t n_fbs = disk::vector_basis_size(m_hho_di.face_degree(),
					   Mesh::dimension - 1,
					   Mesh::dimension);
    
    size_t n_s_cbs = get_a_cell_basis_data();
    size_t n_s_fbs = disk::scalar_basis_size(m_hho_di.face_degree(), Mesh::dimension - 1);
    
    m_n_elastic_cell_dof = (n_cbs * m_e_material.size());
    m_n_acoustic_cell_dof = (n_s_cbs * m_a_material.size());
    
    m_n_elastic_face_dof = (n_fbs * n_e_edges);
    m_n_acoustic_face_dof = (n_s_fbs * n_a_edges);
    size_t system_size = m_n_elastic_cell_dof + m_n_acoustic_cell_dof
                                              + m_n_elastic_face_dof + m_n_acoustic_face_dof;
    
    LHS = SparseMatrix<T>( system_size, system_size );
    RHS = Matrix<T, Dynamic, 1>::Zero( system_size );
    MASS = SparseMatrix<T>( system_size, system_size );
    COUPLING = SparseMatrix<T>( system_size, system_size );
    
    classify_cells(msh);
    build_cells_maps();
  }

  #include "elastoacoustic_four_fields_assembler/get.hpp"
  #include "elastoacoustic_four_fields_assembler/setting.hpp"
  #include "elastoacoustic_four_fields_assembler/assemble.hpp"
  #include "elastoacoustic_four_fields_assembler/projection.hpp"
  #include "elastoacoustic_four_fields_assembler/scatter_data.hpp"
  #include "elastoacoustic_four_fields_assembler/hho_operators.hpp"
  #include "elastoacoustic_four_fields_assembler/gather_dof_data.hpp"
  #include "elastoacoustic_four_fields_assembler/boundary_conditions.hpp"
  #include "elastoacoustic_four_fields_assembler/material_properties.hpp" 
    
};
  
#endif /* elastoacoustic_four_fields_assembler_hpp */
