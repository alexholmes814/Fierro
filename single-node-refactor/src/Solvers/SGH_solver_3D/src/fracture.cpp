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
#include "fracture.h"

using namespace mtr; // matar namespace
/* 
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

// This function calculates shape function values and their derivatives at point (xi1, xi2, xi3)
// Values are stored in the input array psi
void MasterShapes(FArray <double> psi, double xi1, double xi2, double xi3) {
    // psi
    psi(0,0) = 0.125*(1-xi1)*(1-xi2)*(1-xi3);
    psi(1,0) = 0.125*(1-xi1)*(1-xi2)*(1+xi3);
    psi(2,0) = 0.125*(1+xi1)*(1-xi2)*(1+xi3);
    psi(3,0) = 0.125*(1+xi1)*(1-xi2)*(1-xi3);
    psi(4,0) = 0.125*(1-xi1)*(1+xi2)*(1-xi3);
    psi(5,0) = 0.125*(1-xi1)*(1+xi2)*(1+xi3);
    psi(6,0) = 0.125*(1+xi1)*(1+xi2)*(1+xi3);
    psi(7,0) = 0.125*(1+xi1)*(1+xi2)*(1-xi3);
    // dpsidxi1
    psi(0,1) = -0.125*(1-xi2)*(1-xi3);
    psi(1,1) = -0.125*(1-xi2)*(1+xi3);
    psi(2,1) = 0.125*(1-xi2)*(1+xi3);
    psi(3,1) = 0.125*(1-xi2)*(1-xi3);
    psi(4,1) = -0.125*(1+xi2)*(1-xi3);
    psi(5,1) = -0.125*(1+xi2)*(1+xi3);
    psi(6,1) = 0.125*(1+xi2)*(1+xi3);
    psi(7,1) = 0.125*(1+xi2)*(1-xi3);
    // dpsidxi2
    psi(0,2) = -0.125*(1-xi1)*(1-xi3);
    psi(1,2) = -0.125*(1-xi1)*(1+xi3);
    psi(2,2) = -0.125*(1+xi1)*(1+xi3);
    psi(3,2) = -0.125*(1+xi1)*(1-xi3);
    psi(4,2) = 0.125*(1-xi1)*(1-xi3);
    psi(5,2) = 0.125*(1-xi1)*(1+xi3);
    psi(6,2) = 0.125*(1+xi1)*(1+xi3);
    psi(7,2) = 0.125*(1+xi1)*(1-xi3);
    // dpsidxi3
    psi(0,3) = -0.125*(1-xi1)*(1-xi2);
    psi(1,3) = 0.125*(1-xi1)*(1-xi2);
    psi(2,3) = 0.125*(1+xi1)*(1-xi2);
    psi(3,3) = -0.125*(1+xi1)*(1-xi2);
    psi(4,3) = -0.125*(1-xi1)*(1+xi2);
    psi(5,3) = 0.125*(1-xi1)*(1+xi2);
    psi(6,3) = 0.125*(1+xi1)*(1+xi2);
    psi(7,3) = -0.125*(1+xi1)*(1+xi2);
}

// This function calcules the material matrix for isotropic linear elasticity based upon E and nu
// Output is C Matrix for isotropic linear elastic solid strain formulation
CArray <double> CMaterial(double E, double nu) {
    auto Cmat = CArray <double> (6,6);
    Cmat.set_values(0);
    int coef = E/((1+nu)*(1-2*nu));
    Cmat(0,0) = coef*(1-nu);
    Cmat(0,1) = coef*nu;
    Cmat(0,2) = coef*nu;
    Cmat(1,0) = coef*nu;
    Cmat(1,1) = coef*(1-nu);
    Cmat(1,2) = coef*nu;
    Cmat(2,0) = coef*nu;
    Cmat(2,1) = coef*nu;
    Cmat(2,2) = coef*(1-nu);
    Cmat(3,3) = coef*(1-2*nu)/2;
    Cmat(4,4) = coef*(1-2*nu)/2;
    Cmat(5,5) = coef*(1-2*nu)/2;
    return Cmat;
}

// This function pulls the current displacements and nodal coordinates based upon
// inputs of element number, nodal coordinate matrix, total displacement vector, and connectivity
// The current location of the nodes for the element are calculated
// Outputs elcoords with cols [x1, x2, x3] and uel with cols [u1, u2, u3]
void ElemCoords(CArray <double> elcoords, CArray <double> uel, int elnum, CArray <double> NODES, CArray <int> CONN, CArray <double> ut, CArray <double> us) {
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 3; j++) {
            elcoords(i,j) = NODES(CONN(elnum,i),j);
            uel(i,j) = ut(3 * CONN(elnum,i) + j) + us(3 * CONN(elnum,i) + j);
        }
    }
}

// This function calculates the jacobian, its determinant, its inverse, displacement gradient, global derivatives, and second PK stress
// for element calculations at a point based on inputs of master shape function derivative values, element displacement vector,
// material matrix, and global nodal coordinates
// Outputs are stored in dpsig, gradu, and Jdet (Jinv is only used for calculating global derivatives and are therefore
// unnecessary to return from this function)
void Gradients(CArray <double> dpsig, CArray <double> gradu, double& detJ, FArray <double> dpsiloc, CArray <double> elcoords, CArray <double> uel, CArray <double> S01, CArray <double> E01, CArray <double> C) {
    // creating and initializing Jacobian array, inverse Jacobian, adjoint matrix (transpose of cofactor matrix), green-lagrange strain matrix, lambda and mu conversion
    auto J = CArray <double> (3,3);
    J.set_values(0);
    auto Jinv = CArray <double> (3,3);
    auto adj = CArray <double> (3,3);
    
    // 1) elcoords and dpsiloc -> J with indices of [ [dx1dxi1 dx2dxi1 dx3dxi1] [dx1dxi2 dx2dxi2 dx3dxi2] [dx1dxi3 dx2dxi3 dx3dxi3] ]
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 8; k++) {
                J(i,j) += elcoords(k,j) * dpsiloc(k,i+1);
            }
        }
    }
    
    // 2) J -> detJ and Jinv
    detJ = J(0,0)*(J(1,1)*J(2,2) - J(2,1)*J(1,2)) - J(0,1)*(J(1,0)*J(2,2) - J(2,0)*J(1,2)) + J(0,2)*(J(1,0)*J(2,1) - J(2,0)*J(1,1));
    adj(0,0) = J(1,1)*J(2,2) - J(2,1)*J(1,2);
    adj(0,1) = -(J(0,1)*J(2,2) - J(2,1)*J(0,2));
    adj(0,2) = J(0,1)*J(1,2) - J(1,1)*J(0,2);
    adj(1,0) = -(J(1,0)*J(2,2) - J(2,0)*J(1,2));
    adj(1,1) = J(0,0)*J(2,2) - J(2,0)*J(0,2);
    adj(1,2) = -(J(0,0)*J(1,2) - J(1,0)*J(0,2));
    adj(2,0) = J(1,0)*J(2,1) - J(2,0)*J(1,1);
    adj(2,1) = -(J(0,0)*J(2,1) - J(2,0)*J(0,1));
    adj(2,2) = J(0,0)*J(1,1) - J(1,0)*J(0,1);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            Jinv(i,j) = adj(i,j) / detJ;
        }
    }

    // 3) Jinv and dpsiloc -> dpsig with cols [dpsidx1 dpsidx2 dpsidx3]
    for (int i = 0; i < 8; i++) {
        dpsig(i,0) = Jinv(0,0)*dpsiloc(i,1) + Jinv(0,1)*dpsiloc(i,2) + Jinv(0,2)*dpsiloc(i,3);
        dpsig(i,1) = Jinv(1,0)*dpsiloc(i,1) + Jinv(1,1)*dpsiloc(i,2) + Jinv(1,2)*dpsiloc(i,3);
        dpsig(i,2) = Jinv(2,0)*dpsiloc(i,1) + Jinv(2,1)*dpsiloc(i,2) + Jinv(2,2)*dpsiloc(i,3);
    }

    // 4) dpsig and uel -> gradu [ [dudx dvdx dwdx] [dudy dvdy dwdy] [dudz dvdz dwdz] ]
    gradu.set_values(0);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 8; k++) {
                gradu(i,j) += uel(k,j) * dpsig(k,i);
            }
        }
    }

    // 5) gradu -> S01 (in Voigt notation [S11 S22 S33 S23 S31 S12])
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                E01(i,j) = 0.5*(gradu(i,j) + gradu(j,i) + gradu(k,i)*gradu(k,j));
            }
        }
    }

    S01(0) = C(0,0) * E01(0,0) + C(0,1) * E01(1,1) + C(0,2) * E01(2,2);
    S01(1) = C(1,0) * E01(0,0) + C(1,1) * E01(1,1) + C(1,2) * E01(2,2);
    S01(2) = C(2,0) * E01(0,0) + C(2,1) * E01(1,1) + C(2,2) * E01(2,2);
    S01(3) = C(3,3) * 2 * E01(1,2);
    S01(4) = C(4,4) * 2 * E01(2,0);
    S01(5) = C(5,5) * 2 * E01(0,1);
}

// This function calculates the element matrices at a single gauss point based upon inputs for Gradients() and volume gauss point
// Outputs are stored in Kel and F01el
void ElemMats(CArray <double> Kel, CArray <double> F01el, CArray <double> dpsig, CArray <double> gradu, double& detJ, CArray <double> elcoords, CArray <double> uel, CArray <double> S01, CArray <double> C) {
    // Define intermediate arrays
    auto B1 = CArray <double> (6,24);
    auto B2 = CArray <double> (9,24);
    auto K1 = CArray <double> (24,24);
    auto K2 = CArray <double> (24,24);
    B1.set_values(0);
    B2.set_values(0);
    K1.set_values(0);
    K2.set_values(0);
    
    // Definine a few variables just for ease of writing equations
    double ux = gradu(0,0);
    double uy = gradu(1,0);
    double uz = gradu(2,0);
    double vx = gradu(0,1);
    double vy = gradu(1,1);
    double vz = gradu(2,1);
    double wx = gradu(0,2);
    double wy = gradu(1,2);
    double wz = gradu(2,2);

    // Calculating B1 and B2
    for (int i = 0; i < 8; i++) {
        B1(0,3*i) = dpsig(i,0)*(1+ux);
        B1(1,3*i) = dpsig(i,1)*uy;
        B1(2,3*i) = dpsig(i,2)*uz;
        B1(3,3*i) = dpsig(i,1)*uz+dpsig(i,2)*uy;
        B1(4,3*i) = dpsig(i,2)*(1+ux) + dpsig(i,0)*uz;
        B1(5,3*i) = dpsig(i,1)*(1+ux) + dpsig(i,0)*uy;

        B1(0,3*i+1) = dpsig(i,0)*vx;
        B1(1,3*i+1) = dpsig(i,1)*(1+vy);
        B1(2,3*i+1) = dpsig(i,2)*vz;
        B1(3,3*i+1) = dpsig(i,1)*vz + dpsig(i,2)*(1+vy);
        B1(4,3*i+1) = dpsig(i,0)*vz + dpsig(i,2)*vx;
        B1(5,3*i+1) = dpsig(i,0)*(1+vy) + dpsig(i,1)*vx;

        B1(0,3*i+2) = dpsig(i,0)*wx;
        B1(1,3*i+2) = dpsig(i,1)*wy;
        B1(2,3*i+2) = dpsig(i,2)*(1+wz);
        B1(3,3*i+2) = dpsig(i,1)*(1+wz) + dpsig(i,2)*wy;
        B1(4,3*i+2) = dpsig(i,0)*(1+wz) + dpsig(i,2)*wx;
        B1(5,3*i+2) = dpsig(i,0)*wy + dpsig(i,1)*wx;

        B2(0,3*i) = dpsig(i,0);
        B2(1,3*i) = dpsig(i,1);
        B2(2,3*i) = dpsig(i,2);
        
        B2(3,3*i+1) = dpsig(i,0);
        B2(4,3*i+1) = dpsig(i,1);
        B2(5,3*i+1) = dpsig(i,2);

        B2(6,3*i+2) = dpsig(i,0);
        B2(7,3*i+2) = dpsig(i,1);
        B2(8,3*i+2) = dpsig(i,2);
    }
    
    // Calculate K1 = B1^T * C * B1
    // i = rows(K1) = 24   j = rows(K1) = 24   k = rows(C) = 6   m = cols(C) = 6
    for (int i = 0; i < 24; i++) {
        for (int j = 0; j < 24; j++) {
            for (int k = 0; k < 6; k++) {
                for (int m = 0; m < 6; m++) {
                    K1(i,j) += B1(k,i) * C(k,m) * B1(m,j);
                }
            }
        }
    }
    
    // Calculate K2 = B2^T * S * B2
    // i = rows(K2) = 24   j = rows(K2) = 24   k = rows(S) = 9   m = cols(S) = 9
    // See equation 9.4.7 in Reddy nonlinear fem book for why S is 9x9 repeating its 3x3
    // ***********************************************************************************
    // ORDERING OF SHEAR COMPONENTS UNCLEAR TO ME IN 3D SO CHECK HERE IF THINGS BREAK
    // ***********************************************************************************
    auto Smat = CArray <double> (9,9);
    Smat.set_values(0);
    Smat(0,0) = S01(0);
    Smat(0,1) = S01(5);
    Smat(0,2) = S01(4);
    Smat(1,0) = S01(5);
    Smat(1,1) = S01(1);
    Smat(1,2) = S01(3);
    Smat(2,0) = S01(4);
    Smat(2,1) = S01(3);
    Smat(2,2) = S01(2);

    Smat(3,3) = S01(0);
    Smat(3,4) = S01(5);
    Smat(3,5) = S01(4);
    Smat(4,3) = S01(5);
    Smat(4,4) = S01(1);
    Smat(4,5) = S01(3);
    Smat(5,3) = S01(4);
    Smat(5,4) = S01(3);
    Smat(5,5) = S01(2);

    Smat(6,6) = S01(0);
    Smat(6,7) = S01(5);
    Smat(6,8) = S01(4);
    Smat(7,6) = S01(5);
    Smat(7,7) = S01(1);
    Smat(7,8) = S01(3);
    Smat(8,6) = S01(4);
    Smat(8,7) = S01(3);
    Smat(8,8) = S01(2);

    for (int i = 0; i < 24; i++) {
        for (int j = 0; j < 24; j++) {
            for (int k = 0; k < 9; k++) {
                for (int m = 0; m < 9; m++) {
                    K2(i,j) += B2(k,i) * Smat(k,m) * B2(m,j);
                }
            }
        }
    }

    // Calculate Kel
    for (int i = 0; i < 24; i++) {
        for (int j = 0; j < 24; j++) {
            Kel(i,j) += detJ*(K1(i,j) + K2(i,j));
        }
    }

    // calculate F01el
    // ***********************************************************************************
    // ORDERING OF SHEAR COMPONENTS UNCLEAR TO ME IN 3D SO CHECK HERE IF THINGS BREAK
    // ***********************************************************************************
    for (int i = 0; i < 24; i++) {
        for (int j = 0; j < 6; j++) {
            F01el(i) += detJ*B1(j,i)*S01(j);
        }
    }
}

// Applies Dirichlet boundary conditions to the global stiffness matrix and force vector
void Dirichlet(CArray <double> Kg, CArray <double> Fg, CArray <double> DBCS, int NDBC, int NNODE) {
    // Adding pseudo force vector to global force vector
    for (int i = 0; i < NDBC; i++) {
        for (int j = 0; j < 3 * NNODE; j++) {
            Fg(j) += -Kg(j,3 * static_cast<int>(DBCS(i,0)) + static_cast<int>(DBCS(i,1))) * DBCS(i,2);
        }
    }
    // Adjusting stiffness matrix and global force vector to reflect the input value
    for (int i = 0; i < NDBC; i++) {
        for (int j = 0; j < 3 * NNODE; j++) {
            Kg(j,3 * static_cast<int>(DBCS(i,0)) + static_cast<int>(DBCS(i,1))) = 0;
            Kg(3 * static_cast<int>(DBCS(i,0)) + static_cast<int>(DBCS(i,1)),j) = 0;
        }
        Kg(3 * static_cast<int>(DBCS(i,0)) + static_cast<int>(DBCS(i,1)),3 * static_cast<int>(DBCS(i,0)) + static_cast<int>(DBCS(i,1))) = 1;
        Fg(3 * static_cast<int>(DBCS(i,0)) + static_cast<int>(DBCS(i,1))) = DBCS(i,2);
    }
}

// Applies point loads to the global force vector
// Only input ONE SINGLE point load for each degree of freedom
void PointLoad(CArray <double> Fg, CArray <double> PLS, int NPL) {
    for (int i = 0; i < NPL; i++) {
        Fg(3 * static_cast<int>(PLS(i,0)) + static_cast<int>(PLS(i,1))) += PLS(i,2);
    }
}

// Applies traction boundary conditions to the global force vector (only includes uniform tractions as of 10/29/2024)
void Traction(CArray <double> Fg, int elem, int patch, CArray <int> CONN, CArray <double> elcoords, double tx1, double tx2, double tx3, FArray <double> gp0, FArray <double> gp1, FArray <double> gp2, FArray <double> gp3) {
    // initializing jacobian arrays for area calculations and area
    auto J0 = CArray <double> (3,3);
    auto J1 = CArray <double> (3,3);
    auto J2 = CArray <double> (3,3);
    auto J3 = CArray <double> (3,3);
    J0.set_values(0);
    J1.set_values(0);
    J2.set_values(0);
    J3.set_values(0);
    auto Jsur = CArray <double> (4);

    // calculating jacobian matrix for each gauss point
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 8; k++) {
                J0(i,j) += elcoords(k,j) * gp0(k,i+1);
                J1(i,j) += elcoords(k,j) * gp1(k,i+1);
                J2(i,j) += elcoords(k,j) * gp2(k,i+1);
                J3(i,j) += elcoords(k,j) * gp3(k,i+1);
            }
        }
    }

    // calculating surface jacobian magnitude
    int i;
    int j;
    switch(patch) {
        case 0:
            i = 1;
            j = 2;
            Jsur(0) = sqrt(pow((J0(i,1)*J0(j,2) - J0(j,1)*J0(i,2)),2) + pow((J0(i,0)*J0(j,2) - J0(j,0)*J0(i,2)),2) + pow((J0(i,0)*J0(j,1) - J0(j,0)*J0(i,1)),2));
            Jsur(1) = sqrt(pow((J1(i,1)*J1(j,2) - J1(j,1)*J1(i,2)),2) + pow((J1(i,0)*J1(j,2) - J1(j,0)*J1(i,2)),2) + pow((J1(i,0)*J1(j,1) - J1(j,0)*J1(i,1)),2));
            Jsur(2) = sqrt(pow((J2(i,1)*J2(j,2) - J2(j,1)*J2(i,2)),2) + pow((J2(i,0)*J2(j,2) - J2(j,0)*J2(i,2)),2) + pow((J2(i,0)*J2(j,1) - J2(j,0)*J2(i,1)),2));
            Jsur(3) = sqrt(pow((J3(i,1)*J3(j,2) - J3(j,1)*J3(i,2)),2) + pow((J3(i,0)*J3(j,2) - J3(j,0)*J3(i,2)),2) + pow((J3(i,0)*J3(j,1) - J3(j,0)*J3(i,1)),2));
            break;
        case 1:
            i = 1;
            j = 2;
            Jsur(0) = sqrt(pow((J0(i,1)*J0(j,2) - J0(j,1)*J0(i,2)),2) + pow((J0(i,0)*J0(j,2) - J0(j,0)*J0(i,2)),2) + pow((J0(i,0)*J0(j,1) - J0(j,0)*J0(i,1)),2));
            Jsur(1) = sqrt(pow((J1(i,1)*J1(j,2) - J1(j,1)*J1(i,2)),2) + pow((J1(i,0)*J1(j,2) - J1(j,0)*J1(i,2)),2) + pow((J1(i,0)*J1(j,1) - J1(j,0)*J1(i,1)),2));
            Jsur(2) = sqrt(pow((J2(i,1)*J2(j,2) - J2(j,1)*J2(i,2)),2) + pow((J2(i,0)*J2(j,2) - J2(j,0)*J2(i,2)),2) + pow((J2(i,0)*J2(j,1) - J2(j,0)*J2(i,1)),2));
            Jsur(3) = sqrt(pow((J3(i,1)*J3(j,2) - J3(j,1)*J3(i,2)),2) + pow((J3(i,0)*J3(j,2) - J3(j,0)*J3(i,2)),2) + pow((J3(i,0)*J3(j,1) - J3(j,0)*J3(i,1)),2));
            break;
        case 2:
            i = 0;
            j = 2;
            Jsur(0) = sqrt(pow((J0(i,1)*J0(j,2) - J0(j,1)*J0(i,2)),2) + pow((J0(i,0)*J0(j,2) - J0(j,0)*J0(i,2)),2) + pow((J0(i,0)*J0(j,1) - J0(j,0)*J0(i,1)),2));
            Jsur(1) = sqrt(pow((J1(i,1)*J1(j,2) - J1(j,1)*J1(i,2)),2) + pow((J1(i,0)*J1(j,2) - J1(j,0)*J1(i,2)),2) + pow((J1(i,0)*J1(j,1) - J1(j,0)*J1(i,1)),2));
            Jsur(2) = sqrt(pow((J2(i,1)*J2(j,2) - J2(j,1)*J2(i,2)),2) + pow((J2(i,0)*J2(j,2) - J2(j,0)*J2(i,2)),2) + pow((J2(i,0)*J2(j,1) - J2(j,0)*J2(i,1)),2));
            Jsur(3) = sqrt(pow((J3(i,1)*J3(j,2) - J3(j,1)*J3(i,2)),2) + pow((J3(i,0)*J3(j,2) - J3(j,0)*J3(i,2)),2) + pow((J3(i,0)*J3(j,1) - J3(j,0)*J3(i,1)),2));
            break;
        case 3:
            i = 0;
            j = 2;
            Jsur(0) = sqrt(pow((J0(i,1)*J0(j,2) - J0(j,1)*J0(i,2)),2) + pow((J0(i,0)*J0(j,2) - J0(j,0)*J0(i,2)),2) + pow((J0(i,0)*J0(j,1) - J0(j,0)*J0(i,1)),2));
            Jsur(1) = sqrt(pow((J1(i,1)*J1(j,2) - J1(j,1)*J1(i,2)),2) + pow((J1(i,0)*J1(j,2) - J1(j,0)*J1(i,2)),2) + pow((J1(i,0)*J1(j,1) - J1(j,0)*J1(i,1)),2));
            Jsur(2) = sqrt(pow((J2(i,1)*J2(j,2) - J2(j,1)*J2(i,2)),2) + pow((J2(i,0)*J2(j,2) - J2(j,0)*J2(i,2)),2) + pow((J2(i,0)*J2(j,1) - J2(j,0)*J2(i,1)),2));
            Jsur(3) = sqrt(pow((J3(i,1)*J3(j,2) - J3(j,1)*J3(i,2)),2) + pow((J3(i,0)*J3(j,2) - J3(j,0)*J3(i,2)),2) + pow((J3(i,0)*J3(j,1) - J3(j,0)*J3(i,1)),2));
            break;
        case 4:
            i = 0;
            j = 1;
            Jsur(0) = sqrt(pow((J0(i,1)*J0(j,2) - J0(j,1)*J0(i,2)),2) + pow((J0(i,0)*J0(j,2) - J0(j,0)*J0(i,2)),2) + pow((J0(i,0)*J0(j,1) - J0(j,0)*J0(i,1)),2));
            Jsur(1) = sqrt(pow((J1(i,1)*J1(j,2) - J1(j,1)*J1(i,2)),2) + pow((J1(i,0)*J1(j,2) - J1(j,0)*J1(i,2)),2) + pow((J1(i,0)*J1(j,1) - J1(j,0)*J1(i,1)),2));
            Jsur(2) = sqrt(pow((J2(i,1)*J2(j,2) - J2(j,1)*J2(i,2)),2) + pow((J2(i,0)*J2(j,2) - J2(j,0)*J2(i,2)),2) + pow((J2(i,0)*J2(j,1) - J2(j,0)*J2(i,1)),2));
            Jsur(3) = sqrt(pow((J3(i,1)*J3(j,2) - J3(j,1)*J3(i,2)),2) + pow((J3(i,0)*J3(j,2) - J3(j,0)*J3(i,2)),2) + pow((J3(i,0)*J3(j,1) - J3(j,0)*J3(i,1)),2));
            break;
        case 5:
            i = 0;
            j = 1;
            Jsur(0) = sqrt(pow((J0(i,1)*J0(j,2) - J0(j,1)*J0(i,2)),2) + pow((J0(i,0)*J0(j,2) - J0(j,0)*J0(i,2)),2) + pow((J0(i,0)*J0(j,1) - J0(j,0)*J0(i,1)),2));
            Jsur(1) = sqrt(pow((J1(i,1)*J1(j,2) - J1(j,1)*J1(i,2)),2) + pow((J1(i,0)*J1(j,2) - J1(j,0)*J1(i,2)),2) + pow((J1(i,0)*J1(j,1) - J1(j,0)*J1(i,1)),2));
            Jsur(2) = sqrt(pow((J2(i,1)*J2(j,2) - J2(j,1)*J2(i,2)),2) + pow((J2(i,0)*J2(j,2) - J2(j,0)*J2(i,2)),2) + pow((J2(i,0)*J2(j,1) - J2(j,0)*J2(i,1)),2));
            Jsur(3) = sqrt(pow((J3(i,1)*J3(j,2) - J3(j,1)*J3(i,2)),2) + pow((J3(i,0)*J3(j,2) - J3(j,0)*J3(i,2)),2) + pow((J3(i,0)*J3(j,1) - J3(j,0)*J3(i,1)),2));
            break;
    }

    // calculating area
    double Area = Jsur(0) + Jsur(1) + Jsur(2) + Jsur(3);
    
    // adding force contribution based on patch number and element number
    switch(patch) {
        case 0:
            // patch: xi=-1
            Fg(3 * CONN(elem,0)) += tx1*Area/4;
            Fg(3 * CONN(elem,1)) += tx1*Area/4;
            Fg(3 * CONN(elem,4)) += tx1*Area/4;
            Fg(3 * CONN(elem,5)) += tx1*Area/4;
            Fg(3 * CONN(elem,0) + 1) += tx2*Area/4;
            Fg(3 * CONN(elem,1) + 1) += tx2*Area/4;
            Fg(3 * CONN(elem,4) + 1) += tx2*Area/4;
            Fg(3 * CONN(elem,5) + 1) += tx2*Area/4;
            Fg(3 * CONN(elem,0) + 2) += tx3*Area/4;
            Fg(3 * CONN(elem,1) + 2) += tx3*Area/4;
            Fg(3 * CONN(elem,4) + 2) += tx3*Area/4;
            Fg(3 * CONN(elem,5) + 2) += tx3*Area/4;
            break;
        case 1:
            // patch: xi=+1
            Fg(3 * CONN(elem,2)) += tx1*Area/4;
            Fg(3 * CONN(elem,3)) += tx1*Area/4;
            Fg(3 * CONN(elem,6)) += tx1*Area/4;
            Fg(3 * CONN(elem,7)) += tx1*Area/4;
            Fg(3 * CONN(elem,2) + 1) += tx2*Area/4;
            Fg(3 * CONN(elem,3) + 1) += tx2*Area/4;
            Fg(3 * CONN(elem,6) + 1) += tx2*Area/4;
            Fg(3 * CONN(elem,7) + 1) += tx2*Area/4;
            Fg(3 * CONN(elem,2) + 2) += tx3*Area/4;
            Fg(3 * CONN(elem,3) + 2) += tx3*Area/4;
            Fg(3 * CONN(elem,6) + 2) += tx3*Area/4;
            Fg(3 * CONN(elem,7) + 2) += tx3*Area/4;
            break;
        case 2:
            // patch: eta=-1
            Fg(3 * CONN(elem,0)) += tx1*Area/4;
            Fg(3 * CONN(elem,1)) += tx1*Area/4;
            Fg(3 * CONN(elem,2)) += tx1*Area/4;
            Fg(3 * CONN(elem,3)) += tx1*Area/4;
            Fg(3 * CONN(elem,0) + 1) += tx2*Area/4;
            Fg(3 * CONN(elem,1) + 1) += tx2*Area/4;
            Fg(3 * CONN(elem,2) + 1) += tx2*Area/4;
            Fg(3 * CONN(elem,3) + 1) += tx2*Area/4;
            Fg(3 * CONN(elem,0) + 2) += tx3*Area/4;
            Fg(3 * CONN(elem,1) + 2) += tx3*Area/4;
            Fg(3 * CONN(elem,2) + 2) += tx3*Area/4;
            Fg(3 * CONN(elem,3) + 2) += tx3*Area/4;
            break;
        case 3:
            // patch: eta=+1
            Fg(3 * CONN(elem,6)) += tx1*Area/4;
            Fg(3 * CONN(elem,7)) += tx1*Area/4;
            Fg(3 * CONN(elem,4)) += tx1*Area/4;
            Fg(3 * CONN(elem,5)) += tx1*Area/4;
            Fg(3 * CONN(elem,6) + 1) += tx2*Area/4;
            Fg(3 * CONN(elem,7) + 1) += tx2*Area/4;
            Fg(3 * CONN(elem,4) + 1) += tx2*Area/4;
            Fg(3 * CONN(elem,5) + 1) += tx2*Area/4;
            Fg(3 * CONN(elem,6) + 2) += tx3*Area/4;
            Fg(3 * CONN(elem,7) + 2) += tx3*Area/4;
            Fg(3 * CONN(elem,4) + 2) += tx3*Area/4;
            Fg(3 * CONN(elem,5) + 2) += tx3*Area/4;
            break;
        case 4:
            // patch: zeta=-1
            Fg(3 * CONN(elem,0)) += tx1*Area/4;
            Fg(3 * CONN(elem,3)) += tx1*Area/4;
            Fg(3 * CONN(elem,4)) += tx1*Area/4;
            Fg(3 * CONN(elem,7)) += tx1*Area/4;
            Fg(3 * CONN(elem,0) + 1) += tx2*Area/4;
            Fg(3 * CONN(elem,3) + 1) += tx2*Area/4;
            Fg(3 * CONN(elem,4) + 1) += tx2*Area/4;
            Fg(3 * CONN(elem,7) + 1) += tx2*Area/4;
            Fg(3 * CONN(elem,0) + 2) += tx3*Area/4;
            Fg(3 * CONN(elem,3) + 2) += tx3*Area/4;
            Fg(3 * CONN(elem,4) + 2) += tx3*Area/4;
            Fg(3 * CONN(elem,7) + 2) += tx3*Area/4;
            break;
        case 5:
            // patch: zeta=+1
            Fg(3 * CONN(elem,2)) += tx1*Area/4;
            Fg(3 * CONN(elem,1)) += tx1*Area/4;
            Fg(3 * CONN(elem,6)) += tx1*Area/4;
            Fg(3 * CONN(elem,5)) += tx1*Area/4;
            Fg(3 * CONN(elem,2) + 1) += tx2*Area/4;
            Fg(3 * CONN(elem,1) + 1) += tx2*Area/4;
            Fg(3 * CONN(elem,6) + 1) += tx2*Area/4;
            Fg(3 * CONN(elem,5) + 1) += tx2*Area/4;
            Fg(3 * CONN(elem,2) + 2) += tx3*Area/4;
            Fg(3 * CONN(elem,1) + 2) += tx3*Area/4;
            Fg(3 * CONN(elem,6) + 2) += tx3*Area/4;
            Fg(3 * CONN(elem,5) + 2) += tx3*Area/4;
            break;
    }
    //prarr(Fg);
}

// Perform Gaussian Elimination on the augmented matrix [A|b], original b remains unchanged
void GaussElim(CArray <double> A, CArray <double> b, CArray <double> x) {
    int n = b.size();
    
    // Forward elimination
    for (int k = 0; k < n - 1; k++) {
        for (int i = k + 1; i < n; i++) {
            double factor = A(i,k) / A(k,k);
            for (int j = k; j < n; j++) {
                A(i,j) -= factor * A(k,j);
            }
            b(i) -= factor * b(k);
        }
    }
    
    // Backward substitution
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += A(i,j) * x(j);
        }
        x(i) = (b(i) - sum) / A(i,i);
    }
}

// This function calculates the stress at a gauss point
void postprocess(CArray <double> elcoords, CArray <double> uel, int k, CArray <double> NODES, CArray <int> CONN, CArray <double> ut, CArray <double> us, CArray <double> dpsig, CArray <double> gradu, double detJ, FArray <double> gp, CArray <double> S01, CArray <double> E01, CArray <double> C) {
    // Running necessary functions to get values to calculate strain and stress in a given element at a point
    ElemCoords(elcoords, uel, k, NODES, CONN, ut, us);
    Gradients(dpsig, gradu, detJ, gp, elcoords, uel, S01, E01, C);
    
}

// This function calculates the global location of gauss points in each element
CArray <double> gpglob(int NEL, double gp, CArray <double> elcoords, CArray <double> uel, CArray <double> NODES, CArray <int> CONN, CArray <double> ut, CArray <double> us, FArray <double> gp0, FArray <double> gp1, FArray <double> gp2, FArray <double> gp3, FArray <double> gp4, FArray <double> gp5, FArray <double> gp6, FArray <double> gp7) {
    // initializing output array
    auto gpcoords = CArray <double> (NEL,24);
    
    // initializing intermediate variables for calculating jacobian
    auto psi = FArray <double> (8,4);
    auto J0 = CArray <double> (3,3);
    auto J1 = CArray <double> (3,3);
    auto J2 = CArray <double> (3,3);
    auto J3 = CArray <double> (3,3);
    auto J4 = CArray <double> (3,3);
    auto J5 = CArray <double> (3,3);
    auto J6 = CArray <double> (3,3);
    auto J7 = CArray <double> (3,3);
    auto dxi = CArray <double> (3);
    auto dx = CArray <double> (3);

    // looping over each element
    for (int i = 0; i < NEL; i++) {
        J0.set_values(0);
        J1.set_values(0);
        J2.set_values(0);
        J3.set_values(0);
        J4.set_values(0);
        J5.set_values(0);
        J6.set_values(0);
        J7.set_values(0);
        ElemCoords(elcoords,uel,i,NODES,CONN,ut,us);
        // calculating jacobians
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                for (int m = 0; m < 8; m++) {
                    J0(j,k) += elcoords(m,k) * gp0(m,j+1);
                    J1(j,k) += elcoords(m,k) * gp1(m,j+1);
                    J2(j,k) += elcoords(m,k) * gp2(m,j+1);
                    J3(j,k) += elcoords(m,k) * gp3(m,j+1);
                    J4(j,k) += elcoords(m,k) * gp4(m,j+1);
                    J5(j,k) += elcoords(m,k) * gp5(m,j+1);
                    J6(j,k) += elcoords(m,k) * gp6(m,j+1);
                    J7(j,k) += elcoords(m,k) * gp7(m,j+1);
                }
            }
        }
        for (int j = 0; j < 8; j++) {
            dx.set_values(0);
            // calculating dxi1, dxi2, dxi3 based upon which gauss point in the element wrt [xi1,xi2,xi3] = [1,1,1]
            // then calculating dx1, dx2, dx3 with matrix multiplication dx = [J]{dxi}
            switch(j) {
                case 0:
                    dxi(0) = 1+gp;
                    dxi(1) = 1+gp;
                    dxi(2) = 1+gp;
                    for (int k = 0; k < 3; k++) {
                        for (int m = 0; m < 3; m++) {
                            dx(k) += J0(k,m)*dxi(m);
                        }
                    }
                    break;
                case 1:
                    dxi(0) = 1+gp;
                    dxi(1) = 1+gp;
                    dxi(2) = 1-gp;
                    for (int k = 0; k < 3; k++) {
                        for (int m = 0; m < 3; m++) {
                            dx(k) += J1(k,m)*dxi(m);
                        }
                    }
                    break;
                case 2:
                    dxi(0) = 1-gp;
                    dxi(1) = 1+gp;
                    dxi(2) = 1-gp;
                    for (int k = 0; k < 3; k++) {
                        for (int m = 0; m < 3; m++) {
                            dx(k) += J2(k,m)*dxi(m);
                        }
                    }
                    break;
                case 3:
                    dxi(0) = 1-gp;
                    dxi(1) = 1+gp;
                    dxi(2) = 1+gp;
                    for (int k = 0; k < 3; k++) {
                        for (int m = 0; m < 3; m++) {
                            dx(k) += J3(k,m)*dxi(m);
                        }
                    }
                    break;
                case 4:
                    dxi(0) = 1+gp;
                    dxi(1) = 1-gp;
                    dxi(2) = 1+gp;
                    for (int k = 0; k < 3; k++) {
                        for (int m = 0; m < 3; m++) {
                            dx(k) += J4(k,m)*dxi(m);
                        }
                    }
                    break;
                case 5:
                    dxi(0) = 1+gp;
                    dxi(1) = 1-gp;
                    dxi(2) = 1-gp;
                    for (int k = 0; k < 3; k++) {
                        for (int m = 0; m < 3; m++) {
                            dx(k) += J5(k,m)*dxi(m);
                        }
                    }
                    break;
                case 6:
                    dxi(0) = 1-gp;
                    dxi(1) = 1-gp;
                    dxi(2) = 1-gp;
                    for (int k = 0; k < 3; k++) {
                        for (int m = 0; m < 3; m++) {
                            dx(k) += J6(k,m)*dxi(m);
                        }
                    }
                    break;
                case 7:
                    dxi(0) = 1-gp;
                    dxi(1) = 1-gp;
                    dxi(2) = 1+gp;
                    for (int k = 0; k < 3; k++) {
                        for (int m = 0; m < 3; m++) {
                            dx(k) += J7(k,m)*dxi(m);
                        }
                    }
                    break;    
            }
            //printf("%.5f  %.5f  %.5f\n", dx(0),dx(1),dx(2));
            for (int k = 0; k < 3; k++) {
                gpcoords(i,3*j+k) = elcoords(6,k) - dx(k);
            }
        }
    }
    return gpcoords;
}
 */
// START OF VISCOELASTIC COHESIVE ZONE FUNCTIONS

cohesive_zones_t::cohesive_zones_t() {
    // constructor for cohesive zones
}

// initialize the identification of cohesive zones
// this is an algorithim for identifying cohesive zones in a mesh
// it loops over all nodal coordinates and identifies overlapping nodal coordinate pairs
void cohesive_zones_t::initialize(Mesh_t& mesh, State_t& State){
    // the following code counts the number of boundary nodes and checks for node overlaps (2 nodes with the same coordinates)
    // this is the beginning step to setting up cohesive zones for fracture
                   
    // counting the number of boundary nodes
    size_t num_bdy_nodes = mesh.num_bdy_nodes;
    std::cout << "Number of boundary nodes: " << num_bdy_nodes << std::endl;
    std::cout << "Total boundary nodes: " << mesh.num_bdy_nodes << std::endl;

    // total number of boundary nodes across all sets
    size_t total_bdy_nodes = 0;
    for (size_t i = 0; i < mesh.num_bdy_sets; ++i) {
        std::cout << "Boundary nodes in set " << i << ": " << mesh.num_bdy_nodes_in_set(i) << std::endl;
    }
    

    size_t overlap_count = 0;

    const double tol = 0.000001; //0.000001; //e-3; // adjust as needed; added just in case coordinate pairs are close but not exactly equal

    for (size_t i = 0; i < mesh.num_bdy_nodes; ++i) { // outer loop goes through each boundary node
        size_t node_i = mesh.bdy_nodes(i);
        for (size_t j = i + 1; j < mesh.num_bdy_nodes; ++j) { // Avoid duplicate pairs and i == j (starts from i +1), for each pair, compares coordinates k = 0, 1, 2
            size_t node_j = mesh.bdy_nodes(j);
            bool overlap = true;

            for (size_t k = 0; k < 3; ++k) {
                double diff = State.node.coords(node_i, k) - State.node.coords(node_j, k);
                if (std::abs(diff) > tol) {
                 overlap = false;
                break;
                }
            }
            if (overlap) {
            std::cout << "Overlap found between node " << node_i << " and node " << node_j << std::endl;
            ++overlap_count;
            }
        }
    }

    std::cout << "Total overlapping node pairs: " << overlap_count << std::endl; //prints overlapping node pairs

    // printing nodal coordinates to check for overlaps 
    for (size_t i = 0; i < mesh.num_bdy_nodes; ++i) {
        size_t node_id = mesh.bdy_nodes(i);
        std::cout << "Node " << node_id << " coordinates: ";
        for (size_t k = 0; k < 3; ++k) { // 3D
            std::cout << State.node.coords(node_id, k) << " ";
        }
        std::cout << std::endl;
    }   


}

/* 

// returns max number of elements any VCZ node is part of connectivity for in order to define VCZ array sizes properly
int elcount(CArray<double> nodes, CArray<int> conn, int ne, int nvcz, CArray<int> vczconn) {
    // initializing variables for intermediate calculations:
    int count0 = 0;
    int count1 = 0;
    int maxel = 0;
    // looping over all VCZ elements
    for (int i = 0; i < nvcz; i++) {
        count0 = 0;
        count1 = 0;
        // find how many elements the nodes are part of the connectivity for
        for (int j = 0; j < ne; j++) {
            for (int k = 0; k < 8; k++) {
                if (vczconn(i,0) == conn(j,k)) {
                    count0 += 1;
                }
                if (vczconn(i,1) == conn(j,k)) {
                    count1 += 1;
                }
            }
        }

        // updating maxel
        if (count0 > maxel) {
            maxel = count0;
        }
        if (count1 > maxel) {
            maxel = count1;
        }
    }
    return maxel;
}

// inputs: nodes, conn, ne, nvcz, vczconn
// outputs: vczelem
CArray<int> elems(CArray<double> nodes, CArray<int> conn, int ne, int nvcz, CArray<int> vczconn, int maxel) {
    // initializing output array
    auto vczelem = CArray <int> (nvcz,2*maxel);
    for (int i = 0; i < nvcz; i++) {
        for (int j = 0; j < 2*maxel; j++) {
            vczelem(i,j) = -1;
        }
    }

    // initialize intermediate variables 
    int count0 = 0;
    int count1 = 0;
    
    // looping over all VCZ elements
    for (int i = 0; i < nvcz; i++) {
        // find how many elements the nodes are part of the connectivity for
        count0 = 0;
        count1 = 0;
        for (int j = 0; j < ne; j++) {
            for (int k = 0; k < 8; k++) {
                if (vczconn(i,0) == conn(j,k)) {
                    vczelem(i,count0) = j;
                    count0 += 1;
                }
                if (vczconn(i,1) == conn(j,k)) {
                    vczelem(i,count1 + maxel) = j;
                    count1 += 1;
                }
            }
        }
    }
    return vczelem;
}

// inputs: NODES, conn, interpvals, patch number, element number
// output: n vector, r vector, s vector, and center of face in physical coordinates
void paravecs(CArray<double> nodes, CArray<int> conn, FArray<double> interpvals, int pn, int en, CArray<double> n, CArray<double> r, CArray<double> s, CArray<double> cenface) {
    // reset cenface, r, and s vectors
    for (int i = 0; i < 3; i++) {
        r(i,0) = 0;
        s(i,0) = 0;
        cenface(i,0) = 0;
    }
    
    // populate interpvals based on patch number, note patch has constant normal vector so any point on the patch is fine
    // pull node numbers for the patch
    auto pnodes = CArray <int> (4,1);
    // pull which master dof are tangent
    auto mastan = CArray <int> (2,1);
    switch (pn) {
    case -1:
        break;
    case 0:
        MasterShapes(interpvals,-1, 0, 0);
        pnodes(0,0) = 0;
        pnodes(1,0) = 1;
        pnodes(2,0) = 4;
        pnodes(3,0) = 5;
        mastan(0,0) = 1;
        mastan(1,0) = 2;
        break;
    case 1:
        MasterShapes(interpvals,1, 0, 0);
        pnodes(0,0) = 2;
        pnodes(1,0) = 3;
        pnodes(2,0) = 6;
        pnodes(3,0) = 7;
        mastan(0,0) = 1;
        mastan(1,0) = 2;
        break;
    case 2:
        MasterShapes(interpvals,0, -1, 0);
        pnodes(0,0) = 0;
        pnodes(1,0) = 1;
        pnodes(2,0) = 2;
        pnodes(3,0) = 3;
        mastan(0,0) = 0;
        mastan(1,0) = 2;
        break;
    case 3:
        MasterShapes(interpvals,0, 1, 0);
        pnodes(0,0) = 4;
        pnodes(1,0) = 5;
        pnodes(2,0) = 6;
        pnodes(3,0) = 7;
        mastan(0,0) = 0;
        mastan(1,0) = 2;
        break;
    case 4:
        MasterShapes(interpvals,0, 0, -1);
        pnodes(0,0) = 0;
        pnodes(1,0) = 3;
        pnodes(2,0) = 4;
        pnodes(3,0) = 7;
        mastan(0,0) = 0;
        mastan(1,0) = 1;
        break;
    case 5:
        MasterShapes(interpvals,0, 0, 1);
        pnodes(0,0) = 1;
        pnodes(1,0) = 2;
        pnodes(2,0) = 5;
        pnodes(3,0) = 6;
        mastan(0,0) = 0;
        mastan(1,0) = 1;
        break;
    }
    
    // calculating tangent vectors and cenface values
    if (pn != -1) {
        for (int i = 0; i < 4; i++) {
            cenface(0,0) += nodes(conn(en,pnodes(i,0)),0) * interpvals(pnodes(i,0),0);
            cenface(1,0) += nodes(conn(en,pnodes(i,0)),1) * interpvals(pnodes(i,0),0);
            cenface(2,0) += nodes(conn(en,pnodes(i,0)),2) * interpvals(pnodes(i,0),0);

            r(0,0) += nodes(conn(en,pnodes(i,0)),0) * interpvals(pnodes(i,0),mastan(0,0) + 1);
            r(1,0) += nodes(conn(en,pnodes(i,0)),1) * interpvals(pnodes(i,0),mastan(0,0) + 1);
            r(2,0) += nodes(conn(en,pnodes(i,0)),2) * interpvals(pnodes(i,0),mastan(0,0) + 1);

            s(0,0) += nodes(conn(en,pnodes(i,0)),0) * interpvals(pnodes(i,0),mastan(1,0) + 1);
            s(1,0) += nodes(conn(en,pnodes(i,0)),1) * interpvals(pnodes(i,0),mastan(1,0) + 1);
            s(2,0) += nodes(conn(en,pnodes(i,0)),2) * interpvals(pnodes(i,0),mastan(1,0) + 1);
        }
    }

    // cross product for normal vector
    if (pn == 1 || pn == 2 || pn == 5) {
        // s cross r
        n(0,0) = s(1,0) * r(2,0) - r(1,0) * s(2,0);
        n(1,0) = s(0,0) * r(2,0) - r(0,0) * s(2,0);
        n(2,0) = s(0,0) * r(1,0) - r(0,0) * s(1,0);
    }
    else {
        n(0,0) = r(1,0) * s(2,0) - s(1,0) * r(2,0);
        n(1,0) = r(0,0) * s(2,0) - s(0,0) * r(2,0);
        n(2,0) = r(0,0) * s(1,0) - s(0,0) * r(1,0);
    }

    // normalizing to direction vectors
    double lenr = sqrt(r(0,0) * r(0,0) + r(1,0) * r(1,0) + r(2,0) * r(2,0));
    double lens = sqrt(s(0,0) * s(0,0) + s(1,0) * s(1,0) + s(2,0) * s(2,0));
    double lenn = sqrt(n(0,0) * n(0,0) + n(1,0) * n(1,0) + n(2,0) * n(2,0));

    for (int i = 0; i < 3; i++) {
        n(i,0) /= lenn;
        r(i,0) /= lenr;
        s(i,0) /= lens;
    }

}
 */

// // inputs: NODES, conn, NE, NVCZ, VCZconn, VCZelem, interpvals
// // outputs: VCZinfo
// CArray<int> faces(CArray<double> nodes, CArray<int> conn, int ne, int nvcz, CArray<int> vczconn, CArray<int> vczelem, FArray<double> interpvals, int maxel) {
//     // initializing output array
//     auto vczinfo = CArray <int> (nvcz,6*maxel);
//     for (int i = 0; i < nvcz; i++) {
//         for (int j = 0; j < 2 * maxel; j++) {
//             vczinfo(i,j) = vczelem(i,j);
//             vczinfo(i,j + 2 * maxel) = -1;
//             vczinfo(i,j + 4 * maxel) = -1;
//         }
//     }
    
//     // initializing intermediate array
//     auto vczfaces = CArray <int> (nvcz,6*maxel);
//     for (int i = 0; i < nvcz; i++) {
//         for (int j = 0; j < 6 * maxel; j++) {
//             vczfaces(i,j) = -1;
//         }
//     }
    
//     // intermediate value
//     int cnt;
    
//     // looping over all VCZ elements to populate VCZfaces based upon VCZelem and connectivities
//     for (int i = 0; i < nvcz; i++) {
//         cnt = 0;
//         for (int j = 0; j < maxel; j++) {
//             if (vczelem(i,j) != -1) {
//                 for (int k = 0; k < 8; k++) {
//                     if (conn(vczelem(i,j),k) == vczconn(i,0)) {
//                         vczinfo(i,4*maxel+cnt) = k;
//                         switch (k) {
//                         case 0:
//                             vczfaces(i,3 * j) = 0;
//                             vczfaces(i,3 * j + 1) = 2;
//                             vczfaces(i,3 * j + 2) = 4;
//                             break;
//                         case 1:
//                             vczfaces(i,3 * j) = 0;
//                             vczfaces(i,3 * j + 1) = 2;
//                             vczfaces(i,3 * j + 2) = 5;
//                             break;
//                         case 2:
//                             vczfaces(i,3 * j) = 1;
//                             vczfaces(i,3 * j + 1) = 2;
//                             vczfaces(i,3 * j + 2) = 5;
//                             break;
//                         case 3:
//                             vczfaces(i,3 * j) = 1;
//                             vczfaces(i,3 * j + 1) = 2;
//                             vczfaces(i,3 * j + 2) = 4;
//                             break;
//                         case 4:
//                             vczfaces(i,3 * j) = 0;
//                             vczfaces(i,3 * j + 1) = 3;
//                             vczfaces(i,3 * j + 2) = 4;
//                             break;
//                         case 5:
//                             vczfaces(i,3 * j) = 0;
//                             vczfaces(i,3 * j + 1) = 3;
//                             vczfaces(i,3 * j + 2) = 5;
//                             break;
//                         case 6:
//                             vczfaces(i,3 * j) = 1;
//                             vczfaces(i,3 * j + 1) = 3;
//                             vczfaces(i,3 * j + 2) = 5;
//                             break;
//                         case 7:
//                             vczfaces(i,3 * j) = 1;
//                             vczfaces(i,3 * j + 1) = 3;
//                             vczfaces(i,3 * j + 2) = 4;
//                             break;
//                         }
//                     }
//                 }
//             }
//             if (vczelem(i,j+maxel) != -1) {
//                 for (int k = 0; k < 8; k++) {
//                     if (conn(vczelem(i,j+maxel),k) == vczconn(i,1)) {
//                         vczinfo(i,5*maxel+cnt) = k;
//                         switch (k) {
//                         case 0:
//                             vczfaces(i,3 * maxel + 3 * j) = 0;
//                             vczfaces(i,3 * maxel + 3 * j + 1) = 2;
//                             vczfaces(i,3 * maxel + 3 * j + 2) = 4;
//                             break;
//                         case 1:
//                             vczfaces(i,3 * maxel + 3 * j) = 0;
//                             vczfaces(i,3 * maxel + 3 * j + 1) = 2;
//                             vczfaces(i,3 * maxel + 3 * j + 2) = 5;
//                             break;
//                         case 2:
//                             vczfaces(i,3 * maxel + 3 * j) = 1;
//                             vczfaces(i,3 * maxel + 3 * j + 1) = 2;
//                             vczfaces(i,3 * maxel + 3 * j + 2) = 5;
//                             break;
//                         case 3:
//                             vczfaces(i,3 * maxel + 3 * j) = 1;
//                             vczfaces(i,3 * maxel + 3 * j + 1) = 2;
//                             vczfaces(i,3 * maxel + 3 * j + 2) = 4;
//                             break;
//                         case 4:
//                             vczfaces(i,3 * maxel + 3 * j) = 0;
//                             vczfaces(i,3 * maxel + 3 * j + 1) = 3;
//                             vczfaces(i,3 * maxel + 3 * j + 2) = 4;
//                             break;
//                         case 5:
//                             vczfaces(i,3 * maxel + 3 * j) = 0;
//                             vczfaces(i,3 * maxel + 3 * j + 1) = 3;
//                             vczfaces(i,3 * maxel + 3 * j + 2) = 5;
//                             break;
//                         case 6:
//                             vczfaces(i,3 * maxel + 3 * j) = 1;
//                             vczfaces(i,3 * maxel + 3 * j + 1) = 3;
//                             vczfaces(i,3 * maxel + 3 * j + 2) = 5;
//                             break;
//                         case 7:
//                             vczfaces(i,3 * maxel + 3 * j) = 1;
//                             vczfaces(i,3 * maxel + 3 * j + 1) = 3;
//                             vczfaces(i,3 * maxel + 3 * j + 2) = 4;
//                             break;
//                         }
//                     }
//                 }
//             }
//             cnt += 1;
//         }
//     }
    
//     // initialize intermediate variables
//     auto nj = CArray <double> (3,1);
//     auto rj = CArray <double> (3,1);
//     auto sj = CArray <double> (3,1);
//     auto nk = CArray <double> (3,1);
//     auto rk = CArray <double> (3,1);
//     auto sk = CArray <double> (3,1);
//     auto cenfacej = CArray <double> (3,1);
//     auto cenfacek = CArray <double> (3,1);
//     for (int i = 0; i < 3; i++) {
//         nj(i,0) = 0;
//         rj(i,0) = 0;
//         sj(i,0) = 0;
//         nk(i,0) = 0;
//         rk(i,0) = 0;
//         sk(i,0) = 0;
//         cenfacej(i,0) = 0;
//         cenfacek(i,0) = 0;
//     }
//     int count = 0;
    
//     // looping over all VCZ elements to populate VCZinfo based upon VCZelem, VCZfaces, and connectivities
//     for (int i = 0; i < nvcz; i++) {
//         count = 0;
//         for (int j = 0; j < 3 * maxel; j++) {
//             // calculating normal vector for j face
//             paravecs(nodes, conn, interpvals, vczfaces(i,j), vczelem(i,static_cast<int>(floor(j / 3))), nj, rj, sj, cenfacej);
//             for (int k = 0; k < 3 * maxel; k++) {
//                 // calculating normal vector for k face
//                 paravecs(nodes, conn, interpvals, vczfaces(i,k + 3 * maxel), vczelem(i,static_cast<int>(floor(k / 3) + maxel)), nk, rk, sk, cenfacek);
//                 // checking that normal vectors are opposing and that the faces are in contact with eachother (center of face overlaps)
//                 if (nj(0,0) == -nk(0,0) && nj(1,0) == -nk(1,0) && nj(2,0) == -nk(2,0) && cenfacej(0,0) == cenfacek(0,0) && cenfacej(1,0) == cenfacek(1,0) && cenfacej(2,0) == cenfacek(2,0)) {
//                     vczinfo(i,count + 2 * maxel) = vczfaces(i,j);
//                     vczinfo(i,3 * maxel + count) = vczfaces(i,k + 3 * maxel);
//                     count += 1;
//                     break;
//                 }
//             }
//             if (vczinfo(i,count + 2 * maxel) != -1) {
//                 break;
//             }
//         }
//     }
    
//     return vczinfo;

// }

// // sequentially calling preprocessing functions
// CArray<int> VCZmeshpreprocess(CArray<double> nodes, CArray<int> conn, int ne, int nvcz, CArray<int> vczconn, FArray<double> interpvals) {
//     int maxel = elcount(nodes, conn, ne, nvcz, vczconn);
//     CArray<int> vczelem = elems(nodes, conn, ne, nvcz, vczconn, maxel);
//     CArray<int> vczinfo = faces(nodes, conn, ne, nvcz, vczconn, vczelem, interpvals, maxel);
//     /* vczinfo columns are each by number of maxel for total size of [nvcz by 6*maxel]
//     [elems first glob node     elems second glob node     patch first glob node     patch second glob node
//     local conn index first glob node     local conn index second glob node] */
//     return vczinfo;
    
// }

// // This function takes reference configuration and total displacements then returns a set of local direction vectors for the VCZ
// // based upon the average orientation of the elements the VCZ nodes are part of
// // inputs: nodes, conn, ne, ut, us, vczconn, vczinfo, nvcz, interpvals
// // outputs: vczorient (updated orientation of VCZs based on current configuration)
// void oriented(CArray<double> nodes, CArray<int> conn, int ne, int nvcz, CArray<double> ut, CArray<double> us, CArray<int> vczinfo, FArray<double> interpvals, CArray<double> vczorient) {
//     // resetting output array
//     vczorient.set_values(0);
    
//     // initializing intermediate variables
//     int ln0;
//     int ln1;
//     int lp0;
//     int lp1;
//     int el0;
//     int el1;
//     int cnt;
//     double n0magt;
//     double n0magtdt;
//     auto currel0t = CArray <double> (8,3);
//     auto currel1t = CArray <double> (8,3);
//     auto currel0tdt = CArray <double> (8,3);
//     auto currel1tdt = CArray <double> (8,3);
//     auto pnodes = CArray <int> (4,1);
//     auto mastan = CArray <int> (2,1);
//     auto r0t = CArray <double> (3);
//     auto s0t = CArray <double> (3);
//     auto r0tdt = CArray <double> (3);
//     auto s0tdt = CArray <double> (3);
//     auto r1t = CArray <double> (3);
//     auto s1t = CArray <double> (3);
//     auto r1tdt = CArray <double> (3);
//     auto s1tdt = CArray <double> (3);
//     auto n0t = CArray <double> (3);
//     auto n0tdt = CArray <double> (3);
//     auto n1t = CArray <double> (3);
//     auto n1tdt = CArray <double> (3);
//     n0t.set_values(0);

//     // looping through each vcz
//     for (int i = 0; i < nvcz; i++) {
//         cnt = 0;
//         // looping over possible number of elements to average the orientation of
//         for (int j = 0; j < vczinfo.dims(1)/6; j++) {
//             // resetting tangent vectors
//             r0t.set_values(0);
//             s0t.set_values(0);
//             r0tdt.set_values(0);
//             s0tdt.set_values(0);

//             // checking that the index is populated with an element reference
//             if (vczinfo(i,j+2*vczinfo.dims(1)/6) != -1) {
//                 // storing element reference information, local node on the elements, and the patch of the element
//                 // only using 0 index for average normal direction of first global node as of 12/20/2024
//                 el0 = vczinfo(i,j);
//                 el1 = vczinfo(i,j+vczinfo.dims(1)/6);
//                 lp0 = vczinfo(i,j+2*vczinfo.dims(1)/6);
//                 lp1 = vczinfo(i,j+3*vczinfo.dims(1)/6);
//                 ln0 = vczinfo(i,j+4*vczinfo.dims(1)/6);
//                 ln1 = vczinfo(i,j+5*vczinfo.dims(1)/6);

//                 // finding current configuration of the elements
//                 for (int k = 0; k < 8; k++) {
//                     for (int m = 0; m < 3; m++) {
//                         currel0t(k,m) = nodes(conn(el0,k),m) + ut(3*conn(el0,k)+m);
//                         currel1t(k,m) = nodes(conn(el1,k),m) + ut(3*conn(el0,k)+m);
//                         currel0tdt(k,m) = nodes(conn(el0,k),m) + ut(3*conn(el0,k)+m) + us(3*conn(el0,k)+m);
//                         currel1tdt(k,m) = nodes(conn(el1,k),m) + ut(3*conn(el0,k)+m) + us(3*conn(el0,k)+m);
//                     }
//                 }

//                 // finding tangent vectors
//                 switch (lp0) {
//                 case 0:
//                     MasterShapes(interpvals,-1,0,0);
//                     pnodes(0,0) = 0;
//                     pnodes(1,0) = 1;
//                     pnodes(2,0) = 4;
//                     pnodes(3,0) = 5;
//                     mastan(0,0) = 1;
//                     mastan(1,0) = 2;
//                     break;
//                 case 1:
//                     MasterShapes(interpvals,1,0,0);
//                     pnodes(0,0) = 2;
//                     pnodes(1,0) = 3;
//                     pnodes(2,0) = 6;
//                     pnodes(3,0) = 7;
//                     mastan(0,0) = 1;
//                     mastan(1,0) = 2;
//                     break;
//                 case 2:
//                     MasterShapes(interpvals,0,-1,0);
//                     pnodes(0,0) = 0;
//                     pnodes(1,0) = 1;
//                     pnodes(2,0) = 2;
//                     pnodes(3,0) = 3;
//                     mastan(0,0) = 0;
//                     mastan(1,0) = 2;
//                     break;
//                 case 3:
//                     MasterShapes(interpvals,0,1,0);
//                     pnodes(0,0) = 4;
//                     pnodes(1,0) = 5;
//                     pnodes(2,0) = 6;
//                     pnodes(3,0) = 7;
//                     mastan(0,0) = 0;
//                     mastan(1,0) = 2;
//                     break;
//                 case 4:
//                     MasterShapes(interpvals,0,0,-1);
//                     pnodes(0,0) = 0;
//                     pnodes(1,0) = 3;
//                     pnodes(2,0) = 4;
//                     pnodes(3,0) = 7;
//                     mastan(0,0) = 0;
//                     mastan(1,0) = 1;
//                     break;
//                 case 5:
//                     MasterShapes(interpvals,0,0,1);
//                     pnodes(0,0) = 1;
//                     pnodes(1,0) = 2;
//                     pnodes(2,0) = 5;
//                     pnodes(3,0) = 6;
//                     mastan(0,0) = 0;
//                     mastan(1,0) = 1;
//                     break;
//                 }
                
//                 for (int k = 0; k < 4; k++) {
//                     // r vector
//                     r0t(0) += currel0t(pnodes(k,0),0) * interpvals(pnodes(k,0),mastan(0,0) + 1);
//                     r0t(1) += currel0t(pnodes(k,0),1) * interpvals(pnodes(k,0),mastan(0,0) + 1);
//                     r0t(2) += currel0t(pnodes(k,0),2) * interpvals(pnodes(k,0),mastan(0,0) + 1);
//                     r0tdt(0) += currel0tdt(pnodes(k,0),0) * interpvals(pnodes(k,0),mastan(0,0) + 1);
//                     r0tdt(1) += currel0tdt(pnodes(k,0),1) * interpvals(pnodes(k,0),mastan(0,0) + 1);
//                     r0tdt(2) += currel0tdt(pnodes(k,0),2) * interpvals(pnodes(k,0),mastan(0,0) + 1);
//                     // s vector
//                     s0t(0) += currel0t(pnodes(k,0),0) * interpvals(pnodes(k,0),mastan(1,0) + 1);
//                     s0t(1) += currel0t(pnodes(k,0),1) * interpvals(pnodes(k,0),mastan(1,0) + 1);
//                     s0t(2) += currel0t(pnodes(k,0),2) * interpvals(pnodes(k,0),mastan(1,0) + 1);
//                     s0tdt(0) += currel0tdt(pnodes(k,0),0) * interpvals(pnodes(k,0),mastan(1,0) + 1);
//                     s0tdt(1) += currel0tdt(pnodes(k,0),1) * interpvals(pnodes(k,0),mastan(1,0) + 1);
//                     s0tdt(2) += currel0tdt(pnodes(k,0),2) * interpvals(pnodes(k,0),mastan(1,0) + 1);
//                 }

//                 // cross product for normal vector
//                 if (lp0 == 0 || lp0 == 3 || lp0 == 4) {
//                     // s cross r
//                     n0t(0) += s0t(1) * r0t(2) - r0t(1) * s0t(2);
//                     n0t(1) += -(s0t(0) * r0t(2) - r0t(0) * s0t(2));
//                     n0t(2) += s0t(0) * r0t(1) - r0t(0) * s0t(1);
//                     n0tdt(0) += s0tdt(1) * r0tdt(2) - r0tdt(1) * s0tdt(2);
//                     n0tdt(1) += -(s0tdt(0) * r0tdt(2) - r0tdt(0) * s0tdt(2));
//                     n0tdt(2) += s0tdt(0) * r0tdt(1) - r0tdt(0) * s0tdt(1);
//                 }
//                 else {
//                     // r cross s
//                     n0t(0) += r0t(1) * s0t(2) - s0t(1) * r0t(2);
//                     n0t(1) += -(r0t(0) * s0t(2) - s0t(0) * r0t(2));
//                     n0t(2) += r0t(0) * s0t(1) - s0t(0) * r0t(1);
//                     n0tdt(0) += r0tdt(1) * s0tdt(2) - s0tdt(1) * r0tdt(2);
//                     n0tdt(1) += -(r0tdt(0) * s0tdt(2) - s0tdt(0) * r0tdt(2));
//                     n0tdt(2) += r0tdt(0) * s0tdt(1) - s0tdt(0) * r0tdt(1);
//                 }
//                 cnt += 1;
//             }
//         }

//         // averaging
//         n0t(0) /= cnt;
//         n0t(1) /= cnt;
//         n0t(2) /= cnt;
//         n0tdt(0) /= cnt;
//         n0tdt(1) /= cnt;
//         n0tdt(2) /= cnt;

//         // normalizing
//         n0magt = sqrt(n0t(0)*n0t(0) + n0t(1)*n0t(1) + n0t(2)*n0t(2));
//         n0magtdt = sqrt(n0tdt(0)*n0tdt(0) + n0tdt(1)*n0tdt(1) + n0tdt(2)*n0tdt(2));
//         //n0magtdt = sqrt(n0tdt(i,0)*n0tdt(i,0) + n0tdt(i,1)*n0tdt(i,1) + n0tdt(i,2)*n0tdt(i,2));
//         for (int j = 0; j < 3; j++) {
//             vczorient(i,j) = n0t(j) / n0magt;
//             vczorient(i,j+3) = n0tdt(j) / n0magtdt;
//         }
//     }
// }



// // calculates the opening of the crack in global coordinates and then maps it to local coordinates using normal vectors in
// // vczorient array
// // inputs: nodes, ut, us, vczorient, vczconn, nvcz
// // output: ulocvcz (4 columns: unt untdt utt utdt)
// void ucmap(CArray<double>& nodes, CArray<double>& ut, CArray<double>& us, CArray<double>& vczorient, CArray<int>& vczconn, int nvcz, CArray<double>& ulocvcz) {
//     // intermediate variables
//     auto uglobt = CArray <double> (3);
//     auto uglobtdt = CArray <double> (3);
//     double umagt = 0;
//     double umagtdt = 0;
    
//     // looping through each vcz
//     for (int i = 0; i < nvcz; i++) {
//         // calculating displacement between vcz nodes in global frame
//         uglobt(0) = (nodes(vczconn(i,1),0) + ut(3*vczconn(i,1))) - (nodes(vczconn(i,0),0) + ut(3*vczconn(i,0)));
//         uglobt(1) = (nodes(vczconn(i,1),1) + ut(3*vczconn(i,1)+1)) - (nodes(vczconn(i,0),1) + ut(3*vczconn(i,0)+1));
//         uglobt(2) = (nodes(vczconn(i,1),2) + ut(3*vczconn(i,1)+2)) - (nodes(vczconn(i,0),2) + ut(3*vczconn(i,0)+2));
//         uglobtdt(0) = (nodes(vczconn(i,1),0) + ut(3*vczconn(i,1)) + us(3*vczconn(i,1))) - (nodes(vczconn(i,0),0) + ut(3*vczconn(i,0)) + us(3*vczconn(i,0)));
//         uglobtdt(1) = (nodes(vczconn(i,1),1) + ut(3*vczconn(i,1)+1) + us(3*vczconn(i,1)+1)) - (nodes(vczconn(i,0),1) + ut(3*vczconn(i,0)+1) + us(3*vczconn(i,0)+1));
//         uglobtdt(2) = (nodes(vczconn(i,1),2) + ut(3*vczconn(i,1)+2) + us(3*vczconn(i,1)+2)) - (nodes(vczconn(i,0),2) + ut(3*vczconn(i,0)+2) + us(3*vczconn(i,0)+2));
        
//         // dotting with normal vector to get normal component of displacement
//         ulocvcz(i,0) = uglobt(0)*vczorient(i,0) + uglobt(1)*vczorient(i,1) + uglobt(2)*vczorient(i,2);
//         ulocvcz(i,2) = uglobtdt(0)*vczorient(i,3) + uglobtdt(1)*vczorient(i,4) + uglobtdt(2)*vczorient(i,5);

//         // calcuting magnitude of global u vectors
//         umagt = sqrt(uglobt(0)*uglobt(0) + uglobt(1)*uglobt(1) + uglobt(2)*uglobt(2));
//         umagtdt = sqrt(uglobtdt(0)*uglobtdt(0) + uglobtdt(1)*uglobtdt(1) + uglobtdt(2)*uglobtdt(2));
        
//         // calculating tangent component of displacement assuming that us* == ur*
//         ulocvcz(i,1) = sqrt(abs(umagt*umagt - ulocvcz(i,0)*ulocvcz(i,0)));
//         ulocvcz(i,3) = sqrt(abs(umagtdt*umagtdt - ulocvcz(i,2)*ulocvcz(i,2)));
//     }
    
// }

// // This function calculates delta_delta_internal variables (ddinvars)
// // Inputs: cohezive zone local displacments (uvcz [un(t), ut(t), un(t)+dun, ut(t)+dut]x3), time step size (delt), cohesive zone properties (Einf, a1, n, Eandrhom, uns, uts),
// // NPT (number of prony terms), nvcz (number of cohesive zones), internal variables [lambdadot(t-delt) , alpha , <tn,tt> , sig_1 through sig_m] (invars),
// // Output: delta_internal variables [lambdadot(t), d_alpha, d_<tn,tt>, d_sig_1 through d_sig_m]xNCZ (dinvars)
// void internalvars(CArray<double> uvcz, double delt, double Einf, double a1, double n, CArray<double> Eandrhom, double uns, double uts, int npt, int nvcz, CArray<double> invars, CArray<double> dinvars) {
    
//     // Initializing variables for calculations
//     double lambdadt = 0;
//     double lambdat = 0;
//     double dadt = 0;
//     double Edelt = 0;
//     double sumsig = 0;
//     double sumsigexp = 0;

//     // Calculating E(delt) NOTE: This is a redundant calculation so long as delt is uniform across time steps that ideally should be moved elsewhere to avoid repeating it unless delt is made to vary across steps
//     Edelt += Einf;
//     for (int i = 0; i < npt; i++) {
//         Edelt += Eandrhom(i,0) * Eandrhom(i,1) * (1 - exp(-delt / Eandrhom(i,1))) / delt;
//     }

//     // Looping over each unique pair
//     for (int i = 0; i < nvcz; i++) {

//         // Calculating the lambda(t)+dlambda and lambda(t) values
//         lambdadt = sqrt((uvcz(i,2) / uns) * (uvcz(i,2) / uns) + (uvcz(i,3) / uts) * (uvcz(i,3) / uts));
//         lambdat = sqrt((uvcz(i,0) / uns) * (uvcz(i,0) / uns) + (uvcz(i,1) / uts) * (uvcz(i,1) / uts));

//         // Initializing and calculating lambda rate for current load step (lambdadot(t))
//         dinvars(i,0) = (lambdadt - lambdat) / delt;

//         // Calculating dalpha/dt for current load step
//         if (dinvars(i,0) > 0) {
//             dadt = a1 * pow(((lambdadt + lambdat) / 2), n);
//         }
//         else {
//             dadt = 0;
//         }

//         // Preventing divide by zero errors in final calculation: the value is set equal to the minimum possible positive value that can be represented by a double in C++
//         if (lambdat == 0) {
//             lambdat = std::numeric_limits<double>::min();
//         }
//         if (lambdadt == 0) {
//             lambdadt = std::numeric_limits<double>::min();
//         }
        
//         // Updating dalpha
//         dinvars(i,1) = dadt * delt;

//         // Updating dsigma values for prony terms
//         for (int j = 0; j < npt; j++) {
//             dinvars(i,4 + j) = exp(-delt / Eandrhom(j,1)) * invars(i,4 + j) + Eandrhom(j,0) * Eandrhom(j,1) * invars(i,0) * (1 - exp(-delt / Eandrhom(j,1)));
//         }

//         // Resetting sumsig and sumsigexp on each iteration of the loop prior to calculating them
//         sumsig = 0;
//         sumsigexp = 0;

//         // Calculating sigma sums and sigma product sums (the deltaE term in the residual traction) in the residual traction terms
//         for (int j = 0; j < npt; j++) {
//             sumsig += dinvars(i,4 + j);
//             sumsigexp += (1 - exp(-delt / Eandrhom(j,1))) * dinvars(i,4 + j);
//         }

//         // Enforcing alpha domain limitations
//         if ((invars(i,1) + dinvars(i,1)) > 1) {
//             dinvars(i,1) = 1 - invars(i,1);
//         }

//         // Updating dtraction values for each cohesive zone along the edge
//         dinvars(i,2) = uvcz(i,2) / (uns * lambdadt) * (1 - (invars(i,1) + dinvars(i,1))) * Edelt * dinvars(i,0) * delt
//             + uvcz(i,2) / (uns * lambdadt) * (1 - (invars(i,1) + dinvars(i,1))) * (Einf * lambdat + sumsig)
//             - uvcz(i,0) / (uns * lambdat) * (1 - (invars(i,1))) * (Einf * lambdat + sumsig)
//             + uvcz(i,2) / (uns * lambdadt) * (1 - (invars(i,1) + dinvars(i,1))) * -sumsigexp;
//         dinvars(i,3) = uvcz(i,3) / (uts * lambdadt) * (1 - (invars(i,1) + dinvars(i,1))) * Edelt * dinvars(i,0) * delt
//             + uvcz(i,3) / (uts * lambdadt) * (1 - (invars(i,1) + dinvars(i,1))) * (Einf * lambdat + sumsig)
//             - uvcz(i,1) / (uts * lambdat) * (1 - (invars(i,1))) * (Einf * lambdat + sumsig)
//             + uvcz(i,3) / (uts * lambdadt) * (1 - (invars(i,1) + dinvars(i,1))) * -sumsigexp;
//     }
// }

// // This function calculates the loads due to the cohesive zone tractions in the global frame of reference
// // inputs: vczorient, vczconn, nvcz, invars, dinvars, vczinfo, nodes, conn, ut, us
// // output: Fvcz (global load vector from cohesive zones), KVCZ (global stiffness contributions from cohesive zones)
// void vczloads(CArray<double> vczorient, CArray<int> vczconn, double nvcz, CArray<double> invars, CArray<double> dinvars, CArray<int> vczinfo, CArray<double> nodes, CArray<int> conn, CArray<double> ut, CArray<double> us, CArray<double> Eandrhom, int npt, double delt, double Einf, double uns, double uts, CArray<double> Fvcz,CArray<double> Kvcz) {
//     // reset Fvcz and Kvcz
//     Fvcz.set_values(0);
//     Kvcz.set_values(0);
    
//     // intermediate variables
//     auto curr = CArray <double> (8,3); // current configuration of an element
//     double gp = 0.5773502691896257;
//     auto gps = CArray <double> (4,3);
//     auto J = CArray <double> (3,3);
//     double area;
//     auto interpvals = FArray <double> (8,4);
//     auto Fn = CArray <double> (3);
//     auto Ft = CArray <double> (3);
//     double tanmag;
//     double udotn;
//     auto uglobtdt = CArray <double> (3);
//     auto tanvec = CArray <double> (3);
//     double kn;
//     double kt;
//     double alpha;
//     double beta;
//     auto kvcz = CArray <double> (6,6);

//     // Edelt calculator for stiffness of truss element values
//     double Edelt = 0;
//     Edelt += Einf;
//     for (int i = 0; i < npt; i++) {
//         Edelt += Eandrhom(i,0) * Eandrhom(i,1) * (1 - exp(-delt / Eandrhom(i,1))) / delt;
//     }

//     // looping over each cohesive zone
//     for (int i = 0; i < nvcz; i++) {
//         // looping over each element for the cohesive zones
//         for (int j = 0; j < vczinfo.dims(1)/6; j++) {
//             // checking that the element index in vczinfo is populated with an element number
//             if (vczinfo(i,j) != -1) {

//                 // calculate current configuration of the element in the index
//                 for (int k = 0; k < 8; k++) {
//                     for (int m = 0; m < 3; m++) {
//                         curr(k,m) = nodes(conn(vczinfo(i,j),k),m) + ut(3 * conn(vczinfo(i,j),k) + m) + us(3 * conn(vczinfo(i,j),k) + m);
//                     }
//                 }
                
//                 // finding gauss points to loop through based on patch number
//                 switch(vczinfo(i,2*vczinfo.dims(1)/6+j)){
//                 case 0:
//                     gps(0,0) = -1;
//                     gps(1,0) = -1;
//                     gps(2,0) = -1;
//                     gps(3,0) = -1;
//                     gps(0,1) = -gp;
//                     gps(1,1) = gp;
//                     gps(2,1) = -gp;
//                     gps(3,1) = gp;
//                     gps(0,2) = -gp;
//                     gps(1,2) = -gp;
//                     gps(2,2) = gp;
//                     gps(3,2) = gp;
//                     break;
//                 case 1:
//                     gps(0,0) = 1;
//                     gps(1,0) = 1;
//                     gps(2,0) = 1;
//                     gps(3,0) = 1;
//                     gps(0,1) = -gp;
//                     gps(1,1) = gp;
//                     gps(2,1) = -gp;
//                     gps(3,1) = gp;
//                     gps(0,2) = -gp;
//                     gps(1,2) = -gp;
//                     gps(2,2) = gp;
//                     gps(3,2) = gp;
//                     break;
//                 case 2:
//                     gps(0,1) = -1;
//                     gps(1,1) = -1;
//                     gps(2,1) = -1;
//                     gps(3,1) = -1;
//                     gps(0,0) = -gp;
//                     gps(1,0) = gp;
//                     gps(2,0) = -gp;
//                     gps(3,0) = gp;
//                     gps(0,2) = -gp;
//                     gps(1,2) = -gp;
//                     gps(2,2) = gp;
//                     gps(3,2) = gp;
//                     break;
//                 case 3:
//                     gps(0,1) = 1;
//                     gps(1,1) = 1;
//                     gps(2,1) = 1;
//                     gps(3,1) = 1;
//                     gps(0,0) = -gp;
//                     gps(1,0) = gp;
//                     gps(2,0) = -gp;
//                     gps(3,0) = gp;
//                     gps(0,2) = -gp;
//                     gps(1,2) = -gp;
//                     gps(2,2) = gp;
//                     gps(3,2) = gp;
//                     break;
//                 case 4:
//                     gps(0,2) = -1;
//                     gps(1,2) = -1;
//                     gps(2,2) = -1;
//                     gps(3,2) = -1;
//                     gps(0,0) = -gp;
//                     gps(1,0) = gp;
//                     gps(2,0) = -gp;
//                     gps(3,0) = gp;
//                     gps(0,1) = -gp;
//                     gps(1,1) = -gp;
//                     gps(2,1) = gp;
//                     gps(3,1) = gp;
//                     break;
//                 case 5:
//                     gps(0,2) = 1;
//                     gps(1,2) = 1;
//                     gps(2,2) = 1;
//                     gps(3,2) = 1;
//                     gps(0,0) = -gp;
//                     gps(1,0) = gp;
//                     gps(2,0) = -gp;
//                     gps(3,0) = gp;
//                     gps(0,1) = -gp;
//                     gps(1,1) = -gp;
//                     gps(2,1) = gp;
//                     gps(3,1) = gp;
//                     break;
//                 }
                
//                 // jacobian calculations for the patch area
//                 area = 0;
//                 for (int k = 0; k < 4; k++) {

//                     // calculating jacobian
//                     MasterShapes(interpvals, gps(k,0),gps(k,1),gps(k,2));
//                     J.set_values(0);
//                     for (int m = 0; m < 3; m++) {
//                         for (int o = 0; o < 3; o++) {
//                             for (int p = 0; p < 8; p++) {
//                                 J(m,o) +=  curr(p,o)*interpvals(p,m+1);
//                             }
//                         }
//                     }
                    
//                     // cross product of jacobian on surface to find area
//                     switch(2*vczinfo.dims(1)/6+j) {
//                     case 0:
//                         area += sqrt(pow((J(1,1)*J(2,2) - J(2,1)*J(1,2)),2) + pow((J(1,0)*J(2,2) - J(2,0)*J(1,2)),2) + pow((J(1,0)*J(2,1) - J(2,0)*J(1,1)),2));
//                         break;
//                     case 1:
//                         area += sqrt(pow((J(1,1)*J(2,2) - J(2,1)*J(1,2)),2) + pow((J(1,0)*J(2,2) - J(2,0)*J(1,2)),2) + pow((J(1,0)*J(2,1) - J(2,0)*J(1,1)),2));
//                         break;
//                     case 2:
//                         area += sqrt(pow((J(0,1)*J(2,2) - J(2,1)*J(0,2)),2) + pow((J(0,0)*J(2,2) - J(2,0)*J(0,2)),2) + pow((J(0,0)*J(2,1) - J(2,0)*J(0,1)),2));
//                         break;
//                     case 3:
//                         area += sqrt(pow((J(0,1)*J(2,2) - J(2,1)*J(0,2)),2) + pow((J(0,0)*J(2,2) - J(2,0)*J(0,2)),2) + pow((J(0,0)*J(2,1) - J(2,0)*J(0,1)),2));
//                         break;
//                     case 4:
//                         area += sqrt(pow((J(0,1)*J(1,2) - J(1,1)*J(0,2)),2) + pow((J(0,0)*J(1,2) - J(1,0)*J(0,2)),2) + pow((J(0,0)*J(1,1) - J(1,0)*J(0,1)),2));
//                         break;
//                     case 5:
//                         area += sqrt(pow((J(0,1)*J(1,2) - J(1,1)*J(0,2)),2) + pow((J(0,0)*J(1,2) - J(1,0)*J(0,2)),2) + pow((J(0,0)*J(1,1) - J(1,0)*J(0,1)),2));
//                         break;
//                     }
//                 }
                
//                 // calculating normal force vector
//                 Fn(0) = (invars(i,2)+dinvars(i,2))*(area/4)*vczorient(i,3);
//                 Fn(1) = (invars(i,2)+dinvars(i,2))*(area/4)*vczorient(i,4);
//                 Fn(2) = (invars(i,2)+dinvars(i,2))*(area/4)*vczorient(i,5);
                
//                 // calculating tangent vector     u_t = u - (u dot n)n
//                 uglobtdt(0) = (nodes(vczconn(i,1),0) + ut(3*vczconn(i,1)) + us(3*vczconn(i,1))) - (nodes(vczconn(i,0),0) + ut(3*vczconn(i,0)) + us(3*vczconn(i,0)));
//                 uglobtdt(1) = (nodes(vczconn(i,1),1) + ut(3*vczconn(i,1)+1) + us(3*vczconn(i,1)+1)) - (nodes(vczconn(i,0),1) + ut(3*vczconn(i,0)+1) + us(3*vczconn(i,0)+1));
//                 uglobtdt(2) = (nodes(vczconn(i,1),2) + ut(3*vczconn(i,1)+2) + us(3*vczconn(i,1)+2)) - (nodes(vczconn(i,0),2) + ut(3*vczconn(i,0)+2) + us(3*vczconn(i,0)+2));
//                 udotn = uglobtdt(0)*vczorient(i,3) + uglobtdt(1)*vczorient(i,4) + uglobtdt(2)*vczorient(i,5);
//                 tanvec(0) = uglobtdt(0) - udotn*vczorient(i,3);
//                 tanvec(1) = uglobtdt(1) - udotn*vczorient(i,4);
//                 tanvec(2) = uglobtdt(2) - udotn*vczorient(i,5);
//                 tanmag = sqrt(tanvec(0)*tanvec(0) + tanvec(1)*tanvec(1) + tanvec(2)*tanvec(2));
//                 if (tanmag != 0) {
//                     tanvec(0) /= tanmag;
//                     tanvec(1) /= tanmag;
//                     tanvec(2) /= tanmag;
//                 }

//                 // calculating tangent force vector
//                 Ft(0) = (invars(i,3)+dinvars(i,3))*(area/4)*tanvec(0);
//                 Ft(1) = (invars(i,3)+dinvars(i,3))*(area/4)*tanvec(1);
//                 Ft(2) = (invars(i,3)+dinvars(i,3))*(area/4)*tanvec(2);
//                 //prarr(Ft);
//                 // adding to global force vector
//                 Fvcz(3*vczconn(i,0)) += Fn(0) + Ft(0);
//                 Fvcz(3*vczconn(i,0)+1) += Fn(1) + Ft(1);
//                 Fvcz(3*vczconn(i,0)+2) += Fn(2) + Ft(2);
//                 Fvcz(3*vczconn(i,1)) -= Fn(0) + Ft(0);
//                 Fvcz(3*vczconn(i,1)+1) -= Fn(1) + Ft(1);
//                 Fvcz(3*vczconn(i,1)+2) -= Fn(2) + Ft(2);

//                 // calculating rotation angles of truss element
//                 // if u_x is zero it causes divide by zero errors
//                 if (uglobtdt(0) == 0) {
                    
//                     if (uglobtdt(1) == 0) {
//                         // u_x and u_y both zero then no rotation wrt z axis
//                         alpha = 0;
//                     } else if (uglobtdt(1) > 0) {
//                         // u_x is zero and u_y is positive, truss orientation wrt z axis is 90 degrees (straight up)
//                         alpha = M_PI/2;
//                     } else {
//                         // u_x is zero and u_y is negative, truss orientation wrt z axis is 270 degrees (straight down)
//                         alpha = 3*M_PI/2;
//                     }

//                     if (uglobtdt(2) == 0) {
//                         // u_x and u_z both zero then no rotation wrt y axis
//                         beta = 0;
//                     } else if (uglobtdt(2) > 0) {
//                         // u_x is zero and u_y is positive, truss orientation wrt y axis is 90 degrees (straight up)
//                         beta = 3*M_PI/2;
//                     } else {
//                         // u_x is zero and u_y is negative, truss orientation wrt y axis is 270 degrees (straight down)
//                         beta = M_PI/2;
//                     }

//                 } else if (uglobtdt(0) > 0) {
//                     // if u_x is positive then the orientation wrt z and y axes is in either quadrant 1 or 4 of unit circle
//                     alpha = atan(uglobtdt(1)/uglobtdt(0));
//                     beta = atan(-uglobtdt(2)/uglobtdt(0));
//                 } else {
//                     // if u_x is negative then the orientation wrt z and y axes is in either quadrant 2 or 3 of unit circle
//                     // this case requires a correction to account for full 360 rotation due to range limits of arctan
//                     alpha = atan(uglobtdt(1)/uglobtdt(0))+M_PI;
//                     beta = atan(-uglobtdt(2)/uglobtdt(0))+M_PI;
//                 }
//                 /* if (uglobtdt(0) != 0) {
//                     alpha = atan(uglobtdt(1)/uglobtdt(0));
//                     beta = atan(-uglobtdt(2)/uglobtdt(0));
//                 } else {
//                     alpha = 0;
//                     beta = 0;
//                 } */
//                 //printf("a = %f    b = %f\n\n",alpha*180/3.14159,beta*180/3.14159);

//                 // calculating directional stiffness in local frame
//                 // note: only one tangent value needed bc of assumption urs == uss
//                 kn = (1 - (invars(i,1) + dinvars(i,1))) * Edelt * (area/4) / uns;
//                 kt = (1 - (invars(i,1) + dinvars(i,1))) * Edelt * (area/4) / uts;

//                 // populating local stiffness matrix based on +ccw turns
//                 // first rotation about x2 axis (beta), second rotation about x3 axis (alpha)
//                 kvcz(0,0) = cos(alpha)*cos(alpha)*cos(beta)*cos(beta)*kn + (cos(alpha)*cos(alpha)*sin(beta)*sin(beta) + sin(alpha)*sin(alpha))*kt;
//                 kvcz(0,1) = cos(alpha)*cos(beta)*cos(beta)*sin(alpha)*kn + (cos(alpha)*sin(alpha)*sin(beta)*sin(beta) - cos(alpha)*sin(alpha))*kt;
//                 kvcz(0,2) = (-cos(alpha)*cos(beta)*sin(beta))*kn + cos(alpha)*cos(beta)*sin(beta)*kt;
//                 kvcz(0,3) = (-cos(alpha)*cos(alpha)*cos(beta)*cos(beta))*kn + (- cos(alpha)*cos(alpha)*sin(beta)*sin(beta) - sin(alpha)*sin(alpha))*kt;
//                 kvcz(0,4) = (-cos(alpha)*cos(beta)*cos(beta)*sin(alpha))*kn + (cos(alpha)*sin(alpha) - cos(alpha)*sin(alpha)*sin(beta)*sin(beta))*kt;
//                 kvcz(0,5) = cos(alpha)*cos(beta)*sin(beta)*kn + (-cos(alpha)*cos(beta)*sin(beta))*kt;
//                 kvcz(1,1) = cos(beta)*cos(beta)*sin(alpha)*sin(alpha)*kn + (cos(alpha)*cos(alpha) + sin(alpha)*sin(alpha)*sin(beta)*sin(beta))*kt;
//                 kvcz(1,2) = (-cos(beta)*sin(alpha)*sin(beta))*kn + cos(beta)*sin(alpha)*sin(beta)*kt;
//                 kvcz(1,3) = (-cos(alpha)*cos(beta)*cos(beta)*sin(alpha))*kn + (cos(alpha)*sin(alpha) - cos(alpha)*sin(alpha)*sin(beta)*sin(beta))*kt;
//                 kvcz(1,4) = (-cos(beta)*cos(beta)*sin(alpha)*sin(alpha))*kn + (- cos(alpha)*cos(alpha) - sin(alpha)*sin(alpha)*sin(beta)*sin(beta))*kt;
//                 kvcz(1,5) = cos(beta)*sin(alpha)*sin(beta)*kn + (-cos(beta)*sin(alpha)*sin(beta))*kt;
//                 kvcz(2,2) = sin(beta)*sin(beta)*kn + cos(beta)*cos(beta)*kt;
//                 kvcz(2,3) = cos(alpha)*cos(beta)*sin(beta)*kn + (-cos(alpha)*cos(beta)*sin(beta))*kt;
//                 kvcz(2,4) = cos(beta)*sin(alpha)*sin(beta)*kn + (-cos(beta)*sin(alpha)*sin(beta))*kt;
//                 kvcz(2,5) = (-sin(beta)*sin(beta))*kn + (-cos(beta)*cos(beta))*kt;
//                 kvcz(3,3) = cos(alpha)*cos(alpha)*cos(beta)*cos(beta)*kn + (cos(alpha)*cos(alpha)*sin(beta)*sin(beta) + sin(alpha)*sin(alpha))*kt;
//                 kvcz(3,4) = cos(alpha)*cos(beta)*cos(beta)*sin(alpha)*kn + (cos(alpha)*sin(alpha)*sin(beta)*sin(beta) - cos(alpha)*sin(alpha))*kt;
//                 kvcz(3,5) = (-cos(alpha)*cos(beta)*sin(beta))*kn + cos(alpha)*cos(beta)*sin(beta)*kt;
//                 kvcz(4,4) = cos(beta)*cos(beta)*sin(alpha)*sin(alpha)*kn + (cos(alpha)*cos(alpha) + sin(alpha)*sin(alpha)*sin(beta)*sin(beta))*kt;
//                 kvcz(4,5) = (-cos(beta)*sin(alpha)*sin(beta))*kn + cos(beta)*sin(alpha)*sin(beta)*kt;
//                 kvcz(5,5) = sin(beta)*sin(beta)*kn + cos(beta)*cos(beta)*kt;

//                 // symmetric fill
//                 kvcz(1,0) = kvcz(0,1);
//                 kvcz(2,0) = kvcz(0,2);
//                 kvcz(2,1) = kvcz(1,2);
//                 kvcz(3,0) = kvcz(0,3);
//                 kvcz(3,1) = kvcz(1,3);
//                 kvcz(3,2) = kvcz(2,3);
//                 kvcz(4,0) = kvcz(0,4);
//                 kvcz(4,1) = kvcz(1,4);
//                 kvcz(4,2) = kvcz(2,4);
//                 kvcz(4,3) = kvcz(3,4);
//                 kvcz(5,0) = kvcz(0,5);
//                 kvcz(5,1) = kvcz(1,5);
//                 kvcz(5,2) = kvcz(2,5);
//                 kvcz(5,3) = kvcz(3,5);
//                 kvcz(5,4) = kvcz(4,5);

//                 // populating Kvcz based on vcz connectivity
//                 for (int m = 0; m < 2; m++) {
//                     for (int o = 0; o < 2; o++) {
//                         for (int p = 0; p < 3; p++) {
//                             for (int q = 0; q < 3; q++) {
//                                 Kvcz(3*vczconn(i,m)+p,3*vczconn(i,o)+q) += kvcz(3*m + p, 3*o + q); 
//                             }
//                         }
//                     }
//                 }

//             }
//         }
//     }

// }

// // This function calls all necessary functions to update the vcz behavior based on a displacement field
// // inputs: nodes, conn, ne, ut, us, vczconn, vczinfo, nvcz, interpvals, vczorient, ulocvcz, npt, vcz material parameters (Einf, a1, n, Eandrhom, uns, uts), delt, invars, dinvars
// // outputs: Fvcz (global load vector due to vczs) and Kvcz (global stiffness matrix to characterize vczs)
// void vczupdate(CArray<double> nodes, CArray<int> conn, int ne, CArray<double> ut, CArray<double> us, CArray<int> vczconn, CArray<int> vczinfo, int nvcz, FArray<double> interpvals, CArray<double> vczorient, CArray<double> ulocvcz, int npt, double Einf, double a1, double n, CArray<double> Eandrhom, double uns, double uts, double delt, CArray<double> invars, CArray<double> dinvars, CArray<double> Fvcz, CArray<double> Kvcz) {
//     // find global orientation normal vectors
//     oriented(nodes,conn,ne,nvcz,ut,us,vczinfo,interpvals,vczorient);
//     // map displacement into local normal and tangent frame of reference;
//     ucmap(nodes,ut,us,vczorient,vczconn,nvcz,ulocvcz);
//     printf(" ");
//     // updating interal variables of cohesive zones
//     internalvars(ulocvcz,delt,Einf,a1,n,Eandrhom,uns,uts,npt,nvcz,invars,dinvars);
//     // populating global force vector and stiffness matrix
//     vczloads(vczorient,vczconn,nvcz,invars,dinvars,vczinfo,nodes,conn,ut,us,Eandrhom,npt,delt,Einf,uns,uts,Fvcz,Kvcz);
// }

// int main()
// {
//     // Initializing input variables from first 2 lines of input file then populating them with readFirstLines()
//     int NNODE;
//     int NEL;
//     int NDBC;
//     int NPL;
//     int NTL;
//     int BCFLAG;
//     int NLS;
//     double E;
//     double nu;
//     double t;
//     double tol;
//     int NUP;
//     double a1;
//     double n;
//     double Einf;
//     double delt;
//     int NPT;
//     double uns;
//     double uts;
//     std::string filename = "Input.txt";
//     readFirstLines(filename, NNODE, NEL, NDBC, NPL, NTL, BCFLAG, NLS, E, nu, t, tol, NUP, a1, n, Einf, delt, NPT, uns, uts);

//     // Initializing input variables based upon sizes from readFirstLines output then populating them with readTheRest()
//     auto NODES = CArray <double> (NNODE,3);
//     auto CONN = CArray <int> (NEL,8);
//     auto Eandrhom = CArray <double> (NPT,2);
//     auto UPs = CArray <int> (NUP,2);
//     int DBCcols;
//     int PLScols;
//     int TLScols;
//     if (BCFLAG == 0) {
//         DBCcols = 3;
//         PLScols = 3;
//         TLScols = 5;
//     }
//     else {
//         DBCcols = 2+NLS;
//         PLScols = 2+NLS;
//         TLScols = 2+(3*NLS);
//     }
//     auto DBCS = CArray <double> (NDBC, DBCcols);
//     auto PLS = CArray <double> (NPL, PLScols);
//     auto TLS = CArray <double> (NTL, TLScols);
//     readTheRest(filename, NUP, NPT, NNODE, NEL, NDBC, NPL, NTL, BCFLAG, NLS, NODES, CONN, DBCS, PLS, TLS, Eandrhom, UPs);
    
//     // Defining gauss point interpolation function values
//     // (gp__ where V means volume, p# means patch number, and the final number is that subset gauss point)
//     // 2 pt gauss quadrature (linear element) so V spans 0-7, p# spans 0-5, each patch spans 0-3, and weights are always 1
//     double gp = 0.5773502691896257;
//     auto gpV0 = FArray <double> (8,4);
//     MasterShapes(gpV0,-gp,-gp,-gp);
//     auto gpV1 = FArray <double> (8,4);
//     MasterShapes(gpV1,-gp,-gp,gp);
//     auto gpV2 = FArray <double> (8,4);
//     MasterShapes(gpV2,gp,-gp,gp);
//     auto gpV3 = FArray <double> (8,4);
//     MasterShapes(gpV3,gp,-gp,-gp);
//     auto gpV4 = FArray <double> (8,4);
//     MasterShapes(gpV4,-gp,gp,-gp);
//     auto gpV5 = FArray <double> (8,4);
//     MasterShapes(gpV5,-gp,gp,gp);
//     auto gpV6 = FArray <double> (8,4);
//     MasterShapes(gpV6,gp,gp,gp);
//     auto gpV7 = FArray <double> (8,4);
//     MasterShapes(gpV7,gp,gp,-gp);
    
//     // patch definition: [0, 1, 2, 3, 4, 5] = [xi1=-1, xi1=+1, xi2=-1, xi2=+1, xi3=-1, xi3=+1]
//     auto gpp00 = FArray <double> (8,4);
//     MasterShapes(gpp00,-1,-gp,-gp);
//     auto gpp01 = FArray <double> (8,4);
//     MasterShapes(gpp01,-1,-gp,gp);
//     auto gpp02 = FArray <double> (8,4);
//     MasterShapes(gpp02,-1,gp,-gp);
//     auto gpp03 = FArray <double> (8,4);
//     MasterShapes(gpp03,-1,gp,gp);

//     auto gpp10 = FArray <double> (8,4);
//     MasterShapes(gpp10,1,-gp,-gp);
//     auto gpp11 = FArray <double> (8,4);
//     MasterShapes(gpp11,1,-gp,gp);
//     auto gpp12 = FArray <double> (8,4);
//     MasterShapes(gpp12,1,gp,-gp);
//     auto gpp13 = FArray <double> (8,4);
//     MasterShapes(gpp13,1,gp,gp);

//     auto gpp20 = FArray <double> (8,4);
//     MasterShapes(gpp20,-gp,-1,-gp);
//     auto gpp21 = FArray <double> (8,4);
//     MasterShapes(gpp21,-gp,-1,gp);
//     auto gpp22 = FArray <double> (8,4);
//     MasterShapes(gpp22,gp,-1,-gp);
//     auto gpp23 = FArray <double> (8,4);
//     MasterShapes(gpp23,gp,-1,gp);

//     auto gpp30 = FArray <double> (8,4);
//     MasterShapes(gpp30,-gp,1,-gp);
//     auto gpp31 = FArray <double> (8,4);
//     MasterShapes(gpp31,-gp,1,gp);
//     auto gpp32 = FArray <double> (8,4);
//     MasterShapes(gpp32,gp,1,-gp);
//     auto gpp33 = FArray <double> (8,4);
//     MasterShapes(gpp33,gp,1,gp);

//     auto gpp40 = FArray <double> (8,4);
//     MasterShapes(gpp40,-gp,-gp,-1);
//     auto gpp41 = FArray <double> (8,4);
//     MasterShapes(gpp41,gp,-gp,-1);
//     auto gpp42 = FArray <double> (8,4);
//     MasterShapes(gpp42,-gp,gp,-1);
//     auto gpp43 = FArray <double> (8,4);
//     MasterShapes(gpp43,gp,gp,-1);

//     auto gpp50 = FArray <double> (8,4);
//     MasterShapes(gpp50,-gp,-gp,1);
//     auto gpp51 = FArray <double> (8,4);
//     MasterShapes(gpp51,gp,-gp,1);
//     auto gpp52 = FArray <double> (8,4);
//     MasterShapes(gpp52,-gp,gp,1);
//     auto gpp53 = FArray <double> (8,4);
//     MasterShapes(gpp53,gp,gp,1);

//     // calculating C matrix for linear elastic bulk elements
//     CArray <double> CMat = CMaterial(E,nu);
    
//     // Initializing necessary arrays and variables to be used during load step -> elements nested loops:
//     // elcoords 8x3, Kg 3NNx3NN, Fg 3NNx1, F02g 3NNx1, F01g 3NNx1, ut 3NNx1, us 3NNx1, dus 3NNx1, Kel 24x24, F01el 24x1,
//     // uel 8x3, S01 6x1, E01 3x3, gradu 3x3, detJ, dpsig 8x3
//     auto elcoords = CArray <double> (8,3);
//     auto Kg = CArray <double> (3*NNODE,3*NNODE);
//     auto Fg = CArray <double> (3*NNODE);
//     auto F02g = CArray <double> (3*NNODE);
//     auto F01g = CArray <double> (3*NNODE);
//     auto ut = CArray <double> (3*NNODE);
//     auto us = CArray <double> (3*NNODE);
//     auto dus = CArray <double> (3*NNODE);
//     auto Kel = CArray <double> (24,24);
//     auto F01el = CArray <double> (24);
//     auto uel = CArray <double> (8,3);
//     auto S01 = CArray <double> (6);
//     auto E01 = CArray <double> (3,3);
//     auto gradu = CArray <double> (3,3);
//     double detJ;
//     auto dpsig = CArray <double> (8,3);

//     // Initialize matrices that see +=
//     Kg.set_values(0);
//     Fg.set_values(0);
//     F02g.set_values(0);
//     F01g.set_values(0);
//     ut.set_values(0);
//     us.set_values(0);
//     dus.set_values(0);
//     Kel.set_values(0);
//     F01el.set_values(0);

//     // calling for mesh preprocessing to handle viscoelastic cohesive zones
//     auto vczconn = CArray <int> (NUP,2);
//     auto interpvals = FArray <double> (8,4);
//     CArray <int> vczinfo = VCZmeshpreprocess(NODES,CONN,NEL,NUP,UPs,interpvals);

//     // allocations for arrays needed to handle cohesive zone behavior
//     auto invars = CArray<double> (NUP,4+NPT);
//     auto dinvars = CArray<double> (NUP,4+NPT);
//     invars.set_values(0);
//     dinvars.set_values(0);
//     auto vczorient = CArray <double> (NUP,6);
//     auto ulocvcz = CArray <double> (NUP,4);
//     auto Fvcz = CArray<double> (3*NNODE);
//     auto Kvcz = CArray<double> (3*NNODE,3*NNODE);

//     // Defining boundary condition stepping arrays
//     int DBsize;
//     int PLsize;
//     int TLsize;
//     // BCFLAG: 0 for uniform stepping, 1 for custom stepping
//     // defining number of columns in step arrays
//     if (BCFLAG == 0) {
//         DBsize = 1;
//         PLsize = 1;
//         TLsize = 3;
//     }
//     else {
//         DBsize = NLS;
//         PLsize = NLS;
//         TLsize = 3*NLS;
//     }
//     // initializing step arrays
//     auto DBstep = CArray <double> (NDBC,DBsize);
//     auto PLstep = CArray <double> (NPL,PLsize);
//     auto TLstep = CArray <double> (NTL,TLsize);
//     // populating step arrays
//     if (BCFLAG == 0) {
//         for (int i = 0; i < NDBC; i++) {
//             DBstep(i,0) = DBCS(i,2) / NLS;
//         }

//         for (int i = 0; i < NPL; i++) {
//             PLstep(i,0) = PLS(i,2) / NLS;
//         }

//         for (int i = 0; i < NTL; i++) {
//             TLstep(i,0) = TLS(i,2) / NLS;
//             TLstep(i,1) = TLS(i,3) / NLS;
//             TLstep(i,2) = TLS(i,4) / NLS;
//         }
//     }
//     else {
//         for (int i = 0; i < NDBC; i++) {
//             for (int j = 0; j < NLS; j++) {
//                 if (j == 0) {
//                     // accounting for first load step assuming that the body is initially undeformed
//                     DBstep(i,j) = DBCS(i,2);
//                 }
//                 else {
//                     DBstep(i,j) = DBCS(i,j + 2) - DBCS(i,j + 1);
//                 }
//             }
//         }

//         for (int i = 0; i < NPL; i++) {
//             for (int j = 0; j < NLS; j++) {
//                 if (j == 0) {
//                     // accounting for first load step assuming that the body is initially undeformed
//                     PLstep(i,j) = PLS(i,2);
//                 }
//                 else {
//                     PLstep(i,j) = PLS(i,j + 2) - PLS(i,j + 1);
//                 }
//             }
//         }

//         for (int i = 0; i < NTL; i++) {
//             for (int j = 0; j < NLS; j++) {
//                 for (int k = 0; k < 3; k++) {
//                     if (j == 0) {
//                         // accounting for first load step assuming that the body is initially undeformed
//                         TLstep(i,k) = TLS(i,2 + k);
//                     }
//                     else {
//                         TLstep(i,3 * j + k) = TLS(i,3 * j + k + 2) - TLS(i,3 * (j - 1) + k + 2);
//                     }
//                 }
//             }
//         }
//     }

//     // Defining arrays for applying boundary conditions on the current load step
//     // PL02 carries current point loads, TL02 carries current traction loads, DB02 carries current prescribed displacement
//     // DBiter carries zeros for maintaining that the displacement step vector satifies DB02 during iteration
//     auto PL02 = CArray <double> (NPL,3);
//     auto TL02 = CArray <double> (NTL,5);
//     auto DB12 = CArray <double> (NDBC,3);
//     auto DBiter = CArray <double> (NDBC,3);
//     PL02.set_values(0);
//     TL02.set_values(0);
//     DB12.set_values(0);
//     DBiter.set_values(0);
//     // pulling node numbers and dofs for point loads and dirichlet, pulling element number and patch number for traction loads
//     for (int i = 0; i < NPL; i++) {
//         for (int j = 0; j < 2; j++) {
//             PL02(i,j) = PLS(i,j);
//         }
//     }
//     for (int i = 0; i < NTL; i++) {
//         for (int j = 0; j < 2; j++) {
//             TL02(i,j) = TLS(i,j);
//         }
//     }
//     for (int i = 0; i < NDBC; i++) {
//         for (int j = 0; j < 2; j++) {
//             DB12(i,j) = DBCS(i,j);
//             DBiter(i,j) = DBCS(i,j);
//         }
//     }

//     // finding global location of gauss points
//     CArray <double> gpcoords = gpglob(NEL,gp,elcoords,uel,NODES,CONN,ut,us,gpV0,gpV1,gpV2,gpV3,gpV4,gpV5,gpV6,gpV7);

//     // Opening Output file
//     std::ofstream outputFile("output.txt");

//     // Writing initial configuration to output file
//     outputFile << "Initial Configuration:" << "\n" << "Node     x     y     z" << "\n";
//     for (int i = 0; i < NNODE; i++) {
//         outputFile << i << "     " << NODES(i,0) << "     " << NODES(i,1) << "     " << NODES(i,2) << "\n";
//     }
//     outputFile << "\n";

//     // Initializing numerator and denominator of enorm equation for simpler calculation
//     double nnorm = 0;
//     double dnorm = 0;

//     // Defining max allowed iterations for a single load step and euclidean norm
//     // Euclidean norm of displacement is initialized at greater than tolerance
//     // so that the first check does not fail and the convergence loop will run
//     int maxiter = 100;
//     double enorm = 2 * tol;

//     // Beginning of load step loop
//     for (int i = 0; i < NLS; i++) {
//         // Updating boundary conditions for the current step from 01 to 02
//         if (BCFLAG == 0) {
//             for (int j = 0; j < NDBC; j++) {
//                 DB12(j,2) = DBstep(j,0);
//             }
//             for (int j = 0; j < NPL; j++) {
//                 PL02(j,2) += PLstep(j,0);
//             }
//             for (int j = 0; j < NTL; j++) {
//                 TL02(j,2) += TLstep(j,0);
//                 TL02(j,3) += TLstep(j,1);
//                 TL02(j,4) += TLstep(j,2);
//             }
//         }
//         else {
//             for (int j = 0; j < NDBC; j++) {
//                 DB12(j,2) = DBstep(j,i);
//             }
//             for (int j = 0; j < NPL; j++) {
//                 PL02(j,2) += PLstep(j,i);
//             }
//             for (int j = 0; j < NTL; j++) {
//                 for (int k = 0; k < 3; k++) {
//                     TL02(j,2 + k) += TLstep(j,3 * i + k);
//                 }
//             }
//         }

//         // Applying the loads for this load step to define F02g for the load step
//         PointLoad(F02g, PL02, NPL);
//         for (int j = 0; j < NTL; j++) {
//             ElemCoords(elcoords, uel, static_cast<int>(TL02(j,0)), NODES, CONN, ut, us);
//             switch(static_cast<int>(TL02(j,1))) {
//             case 0:
//                 Traction(F02g, static_cast<int>(TL02(j,0)), static_cast<int>(TL02(j,1)), CONN, elcoords, TL02(j,2), TL02(j,3), TL02(j,4), gpp00, gpp01, gpp02, gpp03);
//                 break;
//             case 1:
//                 Traction(F02g, static_cast<int>(TL02(j,0)), static_cast<int>(TL02(j,1)), CONN, elcoords, TL02(j,2), TL02(j,3), TL02(j,4), gpp10, gpp11, gpp12, gpp13);
//                 break;
//             case 2:
//                 Traction(F02g, static_cast<int>(TL02(j,0)), static_cast<int>(TL02(j,1)), CONN, elcoords, TL02(j,2), TL02(j,3), TL02(j,4), gpp20, gpp21, gpp22, gpp23);
//                 break;
//             case 3:
//                 Traction(F02g, static_cast<int>(TL02(j,0)), static_cast<int>(TL02(j,1)), CONN, elcoords, TL02(j,2), TL02(j,3), TL02(j,4), gpp30, gpp31, gpp32, gpp33);
//                 break;
//             case 4:
//                 Traction(F02g, static_cast<int>(TL02(j,0)), static_cast<int>(TL02(j,1)), CONN, elcoords, TL02(j,2), TL02(j,3), TL02(j,4), gpp40, gpp41, gpp42, gpp43);
//                 break;
//             case 5:
//                 Traction(F02g, static_cast<int>(TL02(j,0)), static_cast<int>(TL02(j,1)), CONN, elcoords, TL02(j,2), TL02(j,3), TL02(j,4), gpp50, gpp51, gpp52, gpp53);
//                 break;
//             }
//         }

//         // resetting dinvars for the load step
//         dinvars.set_values(0);
        
//         // initializing cohesive zone vector and array for first iteration
//         vczupdate(NODES,CONN,NEL,ut,us,UPs,vczinfo,NUP,interpvals,vczorient,ulocvcz,NPT,Einf,a1,n,Eandrhom,uns,uts,delt,invars,dinvars,Fvcz,Kvcz);
        
//         // enforcing dirichlet conditions on us
//         for (int k = 0; k < NDBC; k++) {
//             us(3 * static_cast<int>(DB12(k,0)) + static_cast<int>(DB12(k,1))) = DB12(k,2);
//         }

//         // Resetting Euclidean norm before initiating convergence check loop so that
//         // the subsequent load steps after the first do not fail to run
//         enorm = 2 * tol;

//         // Beginning of convergence check loop
//         for (int j = 0; j < maxiter; j++) {
//             // Checking for convergence
//             if (enorm < tol) {
//                 // uncomment to print iteration number to console
//                 //printf("j is %i\n", j);
//                 break;
//             }
            
//             //printf("LS is %i\nIT is %i\n\n",i,j);

//             // Resetting global arrays
//             Kg.set_values(0);
//             F01g.set_values(0);
//             Fg.set_values(0);

//             // Resetting enorm variables
//             nnorm = 0;
//             dnorm = 0;

//             // Beginning of element loop
//             for (int k = 0; k < NEL; k++) {
//                 // Pulling element geometry
//                 ElemCoords(elcoords, uel, k, NODES, CONN, ut, us);

//                 // Resetting element arrays
//                 Kel.set_values(0);
//                 F01el.set_values(0);

//                 // Calculating element matrices
//                 // Each gauss point is added to overall element matrices and so doesn't need to be reset between function calls
//                 Gradients(dpsig, gradu, detJ, gpV0, elcoords, uel, S01, E01, CMat);
//                 ElemMats(Kel,F01el,dpsig,gradu,detJ,elcoords,uel,S01,CMat);
//                 Gradients(dpsig, gradu, detJ, gpV1, elcoords, uel, S01, E01, CMat);
//                 ElemMats(Kel,F01el,dpsig,gradu,detJ,elcoords,uel,S01,CMat);
//                 Gradients(dpsig, gradu, detJ, gpV2, elcoords, uel, S01, E01, CMat);
//                 ElemMats(Kel,F01el,dpsig,gradu,detJ,elcoords,uel,S01,CMat);
//                 Gradients(dpsig, gradu, detJ, gpV3, elcoords, uel, S01, E01, CMat);
//                 ElemMats(Kel,F01el,dpsig,gradu,detJ,elcoords,uel,S01,CMat);
//                 Gradients(dpsig, gradu, detJ, gpV4, elcoords, uel, S01, E01, CMat);
//                 ElemMats(Kel,F01el,dpsig,gradu,detJ,elcoords,uel,S01,CMat);
//                 Gradients(dpsig, gradu, detJ, gpV5, elcoords, uel, S01, E01, CMat);
//                 ElemMats(Kel,F01el,dpsig,gradu,detJ,elcoords,uel,S01,CMat);
//                 Gradients(dpsig, gradu, detJ, gpV6, elcoords, uel, S01, E01, CMat);
//                 ElemMats(Kel,F01el,dpsig,gradu,detJ,elcoords,uel,S01,CMat);
//                 Gradients(dpsig, gradu, detJ, gpV7, elcoords, uel, S01, E01, CMat);
//                 ElemMats(Kel,F01el,dpsig,gradu,detJ,elcoords,uel,S01,CMat);
                
//                 // Assembly to global matrices
//                 for (int m = 0; m < 8; m++) {
//                     for (int o = 0; o < 8; o++) {
//                         Kg(3*CONN(k,m),3*CONN(k,o)) += Kel(3*m,3*o);
//                         Kg(3*CONN(k,m),3*CONN(k,o)+1) += Kel(3*m,3*o + 1);
//                         Kg(3*CONN(k,m),3*CONN(k,o)+2) += Kel(3*m,3*o + 2);
//                         Kg(3*CONN(k,m)+1,3*CONN(k,o)) += Kel(3*m+1,3*o);
//                         Kg(3*CONN(k,m)+1,3*CONN(k,o)+1) += Kel(3*m+1,3*o + 1);
//                         Kg(3*CONN(k,m)+1,3*CONN(k,o)+2) += Kel(3*m+1,3*o + 2);
//                         Kg(3*CONN(k,m)+2,3*CONN(k,o)) += Kel(3*m+2,3*o);
//                         Kg(3*CONN(k,m)+2,3*CONN(k,o)+1) += Kel(3*m+2,3*o + 1);
//                         Kg(3*CONN(k,m)+2,3*CONN(k,o)+2) += Kel(3*m+2,3*o + 2);
//                     }
//                 }
//                 for (int m = 0; m < 8; m++) {
//                     F01g(3*CONN(k,m)) += F01el(3*m);
//                     F01g(3*CONN(k,m) + 1) += F01el(3*m + 1);
//                     F01g(3*CONN(k,m) + 2) += F01el(3*m + 2);
//                 }
//             }

//             // adding vcz stiffness contributions
//             for (int k = 0; k < 3*NNODE; k++) {
//                 for (int m = 0; m < 3*NNODE; m++) {
//                     Kg(k,m) += Kvcz(k,m);
//                 }
//             }
            
//             // Forming full Fg
//             for (int k = 0; k < 3 * NNODE; k++) {
//                 Fg(k) = F02g(k) - F01g(k) + Fvcz(k);
//             }
            
//             // Applying Dirichlet boundary conditions
//             Dirichlet(Kg,Fg,DBiter,NDBC,NNODE);

//             // Solving system of equations
//             GaussElim(Kg,Fg,dus);
            
//             // updating us
//             for (int k = 0; k < 3 * NNODE; k++) {
//                 us(k) += dus(k);
//             }

//             // updating vcz for next iteration
//             vczupdate(NODES,CONN,NEL,ut,us,UPs,vczinfo,NUP,interpvals,vczorient,ulocvcz,NPT,Einf,a1,n,Eandrhom,uns,uts,delt,invars,dinvars,Fvcz,Kvcz);

//             // solving for error norm
//             for (int k = 0; k < 3 * NNODE; k++) {
//                 nnorm += dus(k)*dus(k);
//                 dnorm += us(k)*us(k);
//             }
//             enorm = sqrt(nnorm/dnorm);

//         }
//         // updating total displacement
//         for (int k = 0; k < 3 * NNODE; k++) {
//             ut(k) += us(k);
//         }

//         // reset displacement step and force vector that holds tractions and point loads
//         us.set_values(0);
//         F02g.set_values(0);

//         // Updating internal variables post load step convergence
//         for (int j = 0; j < NUP; j++) {
//             invars(j,0) = dinvars(j,0);
//             invars(j,1) += dinvars(j,1);
//             invars(j,2) += dinvars(j,2);
//             invars(j,3) += dinvars(j,3);
//             for (int k = 0; k < NPT; k++) {
//                 invars(j,4 + k) = dinvars(j,4 + k);
//             }
//         }

//         // printing cohesive zone values for output checking (probably should add to output file in so way at some point?)
//         /* for (int j = 0; j < 1; j++) {
//             printf("%f     %f\n",invars(j,1),invars(j,2));
//         } */
//         printf("%f     %f\n",invars(0,1),invars(0,2));
//         //std::cout << invars(0,1) << "     " << invars(0,2) << std::endl;
//         //std::cout << invars.dims(0) << "     " << invars.dims(1) << std::endl;
//         //printf("\n");
//         //prarr(invars);


//         // Post processing and writing to output
//         // Outputting displacement field
//         outputFile << "Load Step " << i << ":" << "\n";
//         outputFile << "Total Displacement:" << "\n" << "Node     x           y              z" << "\n";
//         for (int j = 0; j < NNODE; j++) {
//             outputFile << j << "     " << std::setprecision(2) << std::scientific << ut(3 * j) << "     " << ut(3 * j + 1) << "     " << ut(3 * j + 2) << "\n";
//         }
//         // Outputting strain field
//         outputFile << "\n" << "x          y          z          Exx           Eyy            Ezz            Exy          Exz         Eyz" << "\n";
//         for (int j = 0; j < NEL; j++) {
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV0,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,0) << "     " << gpcoords(j,1) << "     " << gpcoords(j,2) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << E01(0,0) << "     " << E01(1,1) << "     " << E01(2,2) << "     " << E01(0,1) << "     " << E01(0,2) << "     " << E01(1,2) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV1,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,3) << "     " << gpcoords(j,4) << "     " << gpcoords(j,5) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << E01(0,0) << "     " << E01(1,1) << "     " << E01(2,2) << "     " << E01(0,1) << "     " << E01(0,2) << "     " << E01(1,2) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV2,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,6) << "     " << gpcoords(j,7) << "     " << gpcoords(j,8) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << E01(0,0) << "     " << E01(1,1) << "     " << E01(2,2) << "     " << E01(0,1) << "     " << E01(0,2) << "     " << E01(1,2) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV3,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,9) << "     " << gpcoords(j,10) << "     " << gpcoords(j,11) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << E01(0,0) << "     " << E01(1,1) << "     " << E01(2,2) << "     " << E01(0,1) << "     " << E01(0,2) << "     " << E01(1,2) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV4,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,12) << "     " << gpcoords(j,13) << "     " << gpcoords(j,14) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << E01(0,0) << "     " << E01(1,1) << "     " << E01(2,2) << "     " << E01(0,1) << "     " << E01(0,2) << "     " << E01(1,2) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV5,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,15) << "     " << gpcoords(j,16) << "     " << gpcoords(j,17) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << E01(0,0) << "     " << E01(1,1) << "     " << E01(2,2) << "     " << E01(0,1) << "     " << E01(0,2) << "     " << E01(1,2) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV6,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,18) << "     " << gpcoords(j,19) << "     " << gpcoords(j,20) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << E01(0,0) << "     " << E01(1,1) << "     " << E01(2,2) << "     " << E01(0,1) << "     " << E01(0,2) << "     " << E01(1,2) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV7,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,21) << "     " << gpcoords(j,22) << "     " << gpcoords(j,23) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << E01(0,0) << "     " << E01(1,1) << "     " << E01(2,2) << "     " << E01(0,1) << "     " << E01(0,2) << "     " << E01(1,2) << "\n";
//         }
//         // Outputting stress field
//         outputFile << "\n" << "x          y          z          Sxx           Syy            Szz             Sxy           Sxz          Syz" << "\n";
//         for (int j = 0; j < NEL; j++) {
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV0,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,0) << "     " << gpcoords(j,1) << "     " << gpcoords(j,2) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << S01(0) << "     " << S01(1) << "     " << S01(2) << "     " << S01(5) << "     " << S01(4) << "     " << S01(3) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV1,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,3) << "     " << gpcoords(j,4) << "     " << gpcoords(j,5) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << S01(0) << "     " << S01(1) << "     " << S01(2) << "     " << S01(5) << "     " << S01(4) << "     " << S01(3) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV2,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,6) << "     " << gpcoords(j,7) << "     " << gpcoords(j,8) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << S01(0) << "     " << S01(1) << "     " << S01(2) << "     " << S01(5) << "     " << S01(4) << "     " << S01(3) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV3,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,9) << "     " << gpcoords(j,10) << "     " << gpcoords(j,11) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << S01(0) << "     " << S01(1) << "     " << S01(2) << "     " << S01(5) << "     " << S01(4) << "     " << S01(3) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV4,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,12) << "     " << gpcoords(j,13) << "     " << gpcoords(j,14) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << S01(0) << "     " << S01(1) << "     " << S01(2) << "     " << S01(5) << "     " << S01(4) << "     " << S01(3) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV5,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,15) << "     " << gpcoords(j,16) << "     " << gpcoords(j,17) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << S01(0) << "     " << S01(1) << "     " << S01(2) << "     " << S01(5) << "     " << S01(4) << "     " << S01(3) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV6,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,18) << "     " << gpcoords(j,19) << "     " << gpcoords(j,20) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << S01(0) << "     " << S01(1) << "     " << S01(2) << "     " << S01(5) << "     " << S01(4) << "     " << S01(3) << "\n";
//             postprocess(elcoords,uel,j,NODES,CONN,ut,us,dpsig,gradu,detJ,gpV7,S01,E01,CMat);
//             outputFile << std::setprecision(4) << std::fixed << gpcoords(j,21) << "     " << gpcoords(j,22) << "     " << gpcoords(j,23) << "     ";
//             outputFile << std::setprecision(3) << std::scientific << S01(0) << "     " << S01(1) << "     " << S01(2) << "     " << S01(5) << "     " << S01(4) << "     " << S01(3) << "\n";
//         }
//     }
// }
