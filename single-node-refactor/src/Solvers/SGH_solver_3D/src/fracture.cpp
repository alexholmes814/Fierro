#include <stdio.h>
#include </home/alexholmes814/MATAR/src/include/matar.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <iomanip>
#include "mesh.h"
#include "state.h"

using namespace mtr; // matar namespace

// This code runs a 3D Total Lagrangian Method Linear Elastic Finite Element Analysis with Viscoelastic Cohezive Zones between elements for modeling fracturebased upon a user provided input file
// The analysis as of Oct 21st, 2024 will consider an isotropic homogenous linear elastic bulk material and neglects body forces
// Interpolation functions are that of the Lagrange Family for an 8-node 1st order brick element
// Written by Gavin Whetstone

// This function reads the intialization lines of the input file based upon input filename and stores
// the necessary values to then read the rest of the input file
// NNODE: number of nodes
// NEL: number of elements
// NDBC: number of dirichlet boundary conditions
// NPL: number of points loads
// NTL: number of traction loads
// BCFLAG: decides type of BC application
// NLS: number of load steps
// E: Young's Modulus
// nu: Poisson's Ratio
// t: thickness
// tol: convergence tolerance
// NUP: number of unique node pairs with cohesive zones between them
// a1: alpha1 parameter in damage evolution law
// n: n parameter in damage evolution law
// Einf: constant term in the prony series
// delt: delta_t value
// NPT: number of prony series terms after Einf
// uns: u_n^* is the characteristic length for the lambda calculation for the VCZ local normal direction
// urt: u_t^* is the characteristic length for the lambda calculation for the VCZ local tangent direction
void readFirstLines(const std::string& filename, int& NNODE, int& NEL, int& NDBC, int& NPL, int& NTL, int& BCFLAG, int& NLS, double& E, double& nu, double& t, double& tol, int& NUP, double& a1, double& n, double& Einf, double& delt, int& NPT, double& uns, double& uts) {
    // Create an input file stream
    std::ifstream inputFile(filename);

    std::string line;

    // Read the first line
    std::getline(inputFile, line);

    // Extract the values from the line into the initialized integers
    std::istringstream iss(line);
    iss >> NNODE >> NEL >> NDBC >> NPL >> NTL >> BCFLAG;

    // Read and extract from second line
    std::getline(inputFile, line);
    iss.clear();
    iss.str(line);
    iss >> NLS >> E >> nu >> t >> tol;

    // Reand extract from third line
    std::getline(inputFile, line);
    iss.clear();
    iss.str(line);
    iss >> NUP >> a1 >> n >> Einf >> delt >> NPT >> uns >> uts;

    inputFile.close();

}
// This function reads the rest of the input file and stores the prony series parameters, nodal coordinates, element connectivity, cohesive zobe unique node pair connectivity, and boundary condition values
// Outputs are stored in Eandrhom, NODES, CONN, DBCS, PLS, UPs, and TLS
void readTheRest(const std::string& filename, int NUP, int NPT, int NNODE, int NEL, int NDBC, int NPL, int NTL, int BCFLAG, int NLS, CArray <double> NODES, CArray <int> CONN, CArray <double> DBCS, CArray <double> PLS, CArray <double> TLS, CArray <double> Eandrhom, CArray <int> UPs) {
    // Create an input file stream
    std::ifstream inputFile(filename);

    // Skip the first second, and third line
    std::string firstLines;
    std::getline(inputFile, firstLines);
    std::getline(inputFile, firstLines);
    std::getline(inputFile, firstLines);

    // Read and process the remaining lines
    std::string line;
    for (int i = 0; i < NPT; i++) {
        std::getline(inputFile, line);
        std::istringstream iss(line);
        iss >> Eandrhom(i,0) >> Eandrhom(i,1);
    }
    for (int i = 0; i < NNODE; i++) {
        std::getline(inputFile, line);
        std::istringstream iss(line);
        iss >> NODES(i,0) >> NODES(i,1) >> NODES(i,2);
    }
    for (int i = 0; i < NEL; i++) {
        std::getline(inputFile, line);
        std::istringstream iss(line);
        iss >> CONN(i,0) >> CONN(i,1) >> CONN(i,2) >> CONN(i,3) >> CONN(i,4) >> CONN(i,5) >> CONN(i,6) >> CONN(i,7);
    }
    for (int i = 0; i < NUP; i++) {
        std::getline(inputFile, line);
        std::istringstream iss(line);
        iss >> UPs(i,0) >> UPs(i,1);
    }
    if (BCFLAG == 0) {
        for (int i = 0; i < NDBC; i++) {
            std::getline(inputFile, line);
            std::istringstream iss(line);
            iss >> DBCS(i,0) >> DBCS(i,1) >> DBCS(i,2);
        }
        for (int i = 0; i < NPL; i++) {
            std::getline(inputFile, line);
            std::istringstream iss(line);
            iss >> PLS(i,0) >> PLS(i,1) >> PLS(i,2);
        }
        for (int i = 0; i < NTL; i++) {
            std::getline(inputFile, line);
            std::istringstream iss(line);
            iss >> TLS(i,0) >> TLS(i,1) >> TLS(i,2) >> TLS(i,3) >> TLS(i,4);
        }
    }
    else {
        for (int i = 0; i < NDBC; i++) {
            std::getline(inputFile, line);
            std::istringstream iss(line);
            for (int j = 0; j < 2 + NLS; j++) {
                iss >> DBCS(i,j);
            }
        }
        for (int i = 0; i < NPL; i++) {
            std::getline(inputFile, line);
            std::istringstream iss(line);
            for (int j = 0; j < 2 + NLS; j++) {
                iss >> PLS(i,j);
            }
        }
        for (int i = 0; i < NTL; i++) {
            std::getline(inputFile, line);
            std::istringstream iss(line);
            iss >> TLS(i,0) >> TLS(i,1);
            for (int j = 0; j < NLS; j++) {
                for (int k = 0; k < 3; k++) {
                    iss >> TLS(i,2 + 3 * j + k);
                }
            }
        }
    }

}