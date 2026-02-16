/**********************************************************************************************
© 2020. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and
to permit others to do so.
This program is open source under the BSD-3 License.
Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:
1.  Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer.
2.  Redistributions in binary form must reproduce the above copyright notice, this list of
conditions and the following disclaimer in the documentation and/or other materials
provided with the distribution.
3.  Neither the name of the copyright holder nor the names of its contributors may be used
to endorse or promote products derived from this software without specific prior
written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**********************************************************************************************/
#ifndef Fracture_H
#define Fracture_H
#include "matar.h"
#include "mesh.h"
#include "state.h"
#include "simulation_parameters.h"

// struct for fracture cohesive zones

struct cohesive_zones_t {
    // member functions defined in this header file and sized inside of the source file
    
    size_t gid; // global node id
    size_t lid; // local node id
    size_t nvcz; // number of actual cohesive zone pairs
    CArrayKokkos <size_t> nodes_gid; // global ids of the nodes in the cohesive zone
    CArrayKokkos <size_t> overlapping_node_gids; // node pairs with overlapping coordinates ; // will need to size this inside of a function in the source file 
    CArrayKokkos <size_t> vczconn; // 2D [num_pairs][2] definition for function finds the max number of elements that any cohesive zone node is part of
    //CArrayKokkos<size_t> cohesive_zone_info;
    CArrayKokkos<int> cz_info;
    size_t max_elem_in_cohesive_zone;
    CArrayKokkos<double> internal_vars_n0; // 1/28/2026 addition: storage for internal vars at t_n


    // Iso-parametric coordinates of the patch nodes (1D array of size mesh.num_nodes_in_surf)
    // For a standard linear hex, xi = [-1.0, 1.0, 1.0, -1.0], eta = [-1.0, -1.0, 1.0, 1.0]
    // For now, these are the same for all surface objects, but should they be different, then remove static and look to
    // cohesive_zones_t::initialize for how to set these values
    CArrayKokkos<double> xi;  // xi coordinates
    CArrayKokkos<double> eta;  // eta coordinates
    static size_t num_nodes_in_surf;  // number of nodes on the surface
    static constexpr size_t max_nodes = 4;  // max number of nodes on the surface; for allocating memory at compile time


    void initialize(Mesh_t& mesh, State_t& State); // in fracture.cpp can go in and say what initialize does
    // would look something like void node_pairs_t::initialize(const Mesh_t &mesh, ...)
    // this is where the algorithim to find the unique node pairs (boundary nodes) will go
    // only thing that should be in sgh_setup.cpp is calling this function

    cohesive_zones_t(); 

    cohesive_zones_t(const ViewCArrayKokkos<double> &points, const ViewCArrayKokkos<double> &vel_points,
                    const ViewCArrayKokkos<double> &internal_force_points,
                    const ViewCArrayKokkos<double> &fracture_force_points, const ViewCArrayKokkos<double> &mass_points_);


    // START OF FRACTURE FUNCTION AND ARRAY DECLARATIONS

    size_t cohesive_zone_elem_count(const CArrayKokkos<size_t>& overlapping_node_gids, const RaggedRightArrayKokkos<size_t>& elems_in_node, const Mesh_t& mesh);

    KOKKOS_FUNCTION
    void compute_face_geometry(const DCArrayKokkos<double> &nodes,
                                const Mesh_t &mesh,
                                const DCArrayKokkos<double> &node_coords,
                                const DCArrayKokkos<size_t> &conn,
                                const size_t surf,
                                const size_t elem,
                                ViewCArrayKokkos<double> &n,
                                ViewCArrayKokkos<double> &r,
                                ViewCArrayKokkos<double> &s,
                                ViewCArrayKokkos<double> &cenface
                            ) const;


    CArrayKokkos<int> build_cohesive_zone_info(
    const Mesh_t& mesh,
    const State_t& state,
    const CArrayKokkos<size_t>& overlapping_node_gids,   
    const size_t max_elem_in_cohesive_zone,              // from cohesive_zone_elem_count()
    const double tol                                     // centroid coincidence tolerance
    );

    CArrayKokkos<int> cohesive_zone_faces(
       const CArrayKokkos<size_t>& overlapping_node_gids,
       const size_t max_elem_in_cohesive_zone
    );


    KOKKOS_FUNCTION
    void oriented(
        Mesh_t& mesh,
        DCArrayKokkos<double>& pos,      // reference  coords (num_nodes x 3) 
        CArrayKokkos<size_t>& overlapping_node_gids, 
        CArrayKokkos<int>& cz_info,      // from build_cohesive_zone_info()
        size_t max_elem_in_cohesive_zone,
        double tol,                 // centroid coincidence tolerance (ABS distance)
        CArrayKokkos<double>& cohesive_zone_orientation       // (nvcz x 6): [nx_t,ny_t,nz_t, nx_tdt,ny_tdt,nz_tdt]
    ); 

    KOKKOS_FUNCTION
    void ucmap(
        const DCArrayKokkos<double>& pos,
        const DCArrayKokkos<double>& vel,
        const CArrayKokkos<double>& cohesive_zone_orientation,
        const CArrayKokkos<size_t>& overlapping_node_gids,
        const double dt_stage, 
        CArrayKokkos<double>& local_opening    // (overlapping_node_gids.dims(0) x 4): [un_t, utan_t, un_tdt, utan_tdt]
    );

    KOKKOS_FUNCTION
    void cohesive_zone_var_update(
        const CArrayKokkos<double>& local_opening,
        const double dt_stage,
        const double time_value, // ADDED IN FOR DEBUGGING
        const CArrayKokkos<size_t>& overlapping_node_gids,
        const RaggedRightArrayKokkos<double>& stress_bc_global_vars,
        const int bdy_set,
        const ViewCArrayKokkos<double>& internal_vars,      // (overlapping_node_gids.dims(0), 4 + num_prony_terms)
        const ViewCArrayKokkos<double>& delta_internal_vars // (overlapping_node_gids.dims(0), 4 + num_prony_terms)
    );

    CArrayKokkos<double> internal_vars;
    CArrayKokkos<double> delta_internal_vars;

    KOKKOS_FUNCTION
    void cohesive_zone_loads(
        Mesh_t &mesh,
        const DCArrayKokkos<double> &pos,
        const CArrayKokkos<size_t> &overlapping_node_gids,
        const CArrayKokkos<double> &cohesive_zone_orientation,
        const CArrayKokkos<int> &cz_info,
        const size_t max_elem_in_cohesive_zone,
        const ViewCArrayKokkos<double> &internal_vars,
        const ViewCArrayKokkos<double> &delta_internal_vars,
        CArrayKokkos<double> &pair_area,
        const ViewCArrayKokkos<double> &F_cz
    ); 

    // END OF FRACTURE FUNCTION AND ARRAY DECLARATIONS
};

#endif