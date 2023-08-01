#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <cmath>
#include <memory>
#include <sstream>
#include <fstream>
#include <list>
#include <getopt.h>


#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
using namespace Eigen;

#include "timecounter.h"
#include "methods/hho"
#include "geometry/geometry.hpp"
#include "boundary_conditions/boundary_conditions.hpp"
#include "output/silo.hpp"

// application common sources
#include "../common/display_settings.hpp"
#include "../common/fitted_geometry_builders.hpp"
#include "../common/linear_solver.hpp"
#include "../common/acoustic_material_data.hpp"
#include "../common/elastic_material_data.hpp"
//#include "../common/scal_vec_analytic_functions.hpp"
#include "../common/scal_vec_analytic_functions_original.hpp"
#include "../common/preprocessor.hpp"
#include "../common/postprocessor.hpp"

// implicit RK schemes
#include "../common/dirk_hho_scheme.hpp"
#include "../common/dirk_butcher_tableau.hpp"

// Prototypes
#include "Prototypes/IHHOFirstOrder.hpp"
#include "Prototypes/IHHOSecondOrder.hpp"
#include "Prototypes/HeterogeneousGar6more2DIHHOFirstOrder.hpp"
#include "Prototypes/HeterogeneousGar6more2DIHHOSecondOrder.hpp"
#include "Prototypes/HeterogeneousIHHOFirstOrder.hpp"



int main(int argc, char **argv){

  //////////////////////////////////
  ////////////////////////////////// Convergence tests
  //////////////////////////////////
  
  //IHHOFirstOrder(argc, argv);
  //IHHOSecondOrder(argc, argv);

  //////////////////////////////////
  ////////////////////////////////// Test case to compare with Gar6more 
  //////////////////////////////////

  //HeterogeneousGar6more2DIHHOFirstOrder(argc, argv);
  //HeterogeneousGar6more2DIHHOSecondOrder(argc, argv);

  //////////////////////////////////
  ////////////////////////////////// Heterogeneous coupling
  //////////////////////////////////

  HeterogeneousIHHOFirstOrder(argc, argv);

}





