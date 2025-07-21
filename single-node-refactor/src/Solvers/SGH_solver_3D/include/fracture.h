#ifndef Fracture_H
#define Fracture_H
#include "matar.h"
#include "mesh.h"
#include "state.h"
#include "simulation_parameters.h"

// struct for fracture cohesive zones

struct cohesive_zones_t {
    size_t gid; // global node id
    CArrayKokkos <size_t> node_pairs; // will need to size this inside of a function in the source file 
    // member functions defined in this header file and sized inside of the source file
    void initialize(Mesh_t& mesh, State_t State); // in fracture.cpp can go in and say what initialize does
    // would look something like void node_pairs_t::initialize(const Mesh_t &mesh, ...)
    // this is where the algorithim to find the unique node pairs (boundary nodes) will go
    // only thing that should be in sgh_setup.cpp is calling this function

    cohesive_zones_t(); 

};



#endif