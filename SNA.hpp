/* -*- c++ -*- -------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov
   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Aidan Thompson, Christian Trott, SNL
------------------------------------------------------------------------- */

//
// SNA.hpp
//
// LGPL Version 2.1 HEADER START
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
//
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA 02110-1301  USA
//
// LGPL Version 2.1 HEADER END
//

//
// Copyright (c) 2019, Regents of the University of Minnesota.
// All rights reserved.
//
// Contributors:
//    Yaser Afshar
//
// Brief: This file is adapted from the LAMMPS software package
//        `lammps/src/SNAP/sna.h` and amended and updated by
//        Yaser Afshar
//


#ifndef SNA_HPP
#define SNA_HPP

#include "helper.hpp"

#include <vector>

/*!
 * \brief Helper indices to be used by the SNA object
 *
 */
struct SNA_ZINDICES
{
  int j1;
  int j2;
  int j;

  int ma1min;
  int ma2max;
  int na;

  int mb1min;
  int mb2max;
  int nb;

  int jju;
};

/*!
 * \brief Helper indioces for the bispectrum components
 *
 */
struct SNA_BINDICES
{
  int j1;
  int j2;
  int j;
};

/*! \class SNA
 * \brief SNAP (Spectral Neighbor Analysis Potential) [1]: a new formulation for interatomic potentials
 *
 * SNAP provides a powerful framework for automated generation of interatomic potentials
 * fit to QM data. SNAP approach uses GAP’s neighbor bispectrum, but replaces Gaussian process with linear regression.
 *
 *
 * \note
 * This implementation is based on the method outlined in Bartok[2], using formulae from VMK[3].
 * For the Clebsch-Gordan coefficients, we convert the VMK half-integral labels
 * a, b, c, alpha, beta, gamma to array offsets j1, j2, j, m1, m2, m
 * using the following relations:
 *
 * j1 = 2*a
 * j2 = 2*b
 * j =  2*c
 *
 * m1 = alpha+a      2*alpha = 2*m1 - j1
 * m2 = beta+b    or 2*beta = 2*m2 - j2
 * m =  gamma+c      2*gamma = 2*m - j
 *
 * in this way:
 *
 * -a <= alpha <= a
 * -b <= beta <= b
 * -c <= gamma <= c
 *
 * becomes:
 *
 * 0 <= m1 <= j1
 * 0 <= m2 <= j2
 * 0 <= m <= j
 *
 * and the requirement that a+b+c be integral implies that j1+j2+j must be even.
 * The requirement that:
 *
 * gamma = alpha+beta
 *
 * becomes:
 *
 * 2*m - j = 2*m1 - j1 + 2*m2 - j2
 *
 * Similarly, for the Wigner U-functions U(J,m,m') we convert the half-integral labels J,m,m' to
 * array offsets j,ma,mb:
 *
 * j = 2*J
 * ma = J+m
 * mb = J+m'
 *
 * so that:
 *
 * 0 <= j <= 2*Jmax
 * 0 <= ma, mb <= j.
 *
 * For the bispectrum components B(J1,J2,J) we convert to:
 *
 * j1 = 2*J1
 * j2 = 2*J2
 * j = 2*J
 *
 * and the requirement:
 *
 * |J1-J2| <= J <= J1+J2, for j1+j2+j integral
 *
 * becomes:
 *
 * |j1-j2| <= j <= j1+j2, for j1+j2+j even integer
 *
 * or
 *
 * j = |j1-j2|, |j1-j2|+2,...,j1+j2-2,j1+j2
 *
 * [1] Thompson, Swiler, Trott, Foiles, Tucker, J Comp Phys, 285, 316 (2015)
 *
 * [2] Albert Bartok-Partay, Doctoral Thesis, Cambrindge University, (2009)
 *
 * [3] D. A. Varshalovich, A. N. Moskalev, and V. K. Khersonskii, World Scientific (1988)
 */
class SNA
{
public:
  /*!
   * \brief Construct a new SNA object
   *
   * \param rfac0_in Parameter in distance to angle conversion
   * \param twojmax_in Band limit for bispectrum components (non-negative integer)
   * \param rmin0_in Parameter in distance (distance units)
   * \param switchflag_in Flag for switching function
   * \param bzeroflag_in Flag indicating whether bzero subtracted from barray or not
   */
  SNA(double const rfac0_in,
      int const twojmax_in,
      double const rmin0_in,
      int const switchflag_in,
      int const bzeroflag_in);

  /*!
   * \brief Destroy the SNA object
   *
   */
  ~SNA();

  /*!
   * \brief Create index list
   *
   */
  void build_indexlist();

  /*!
   * \brief Initialize the object
   *
   */
  void init();

  // functions for bispectrum coefficients

  /*!
   * \brief Compute Ui by summing over \c jnum neighbors
   *
   * Calculates all expansion coefficients for an atom i
   *
   * \param nneigh Number of neighbors
   */
  void compute_ui(int const nneigh);

  /*!
   * \brief Compute Zi by summing over products of Ui
   *
   *
   */
  void compute_zi();

  /*!
   * \brief Compute Yi from Ui without storing Zi, looping over zlist indices
   *
   * \param beta
   */
  void compute_yi(double const *const beta);

  /*!
   * \brief Compute Bi by summing conj(Ui)*Zi
   *
   */
  void compute_bi();

  // functions for derivatives

  /*!
   * \brief Calculate derivative of Ui w.r.t. atom j
   *
   * Calculate the derivatives of ui with respect to the distance vector between
   * atoms \c i and \c j
   *
   * \param rij_in The distance vector between atoms \c i and \c j
   * \param wj_in Dimensionless weights that are chosen to distinguish atoms of different types
   * \param rcut Cutoff distance
   * \param neigh_j Neighbor index
   */
  void compute_duidrj(double const *const rij_in, double const wj_in, double const rcut, int const neigh_j);

  /*!
   * \brief Calculate derivative of Bi w.r.t. atom j
   *
   * variant using indexlist for j1,j2,j
   * variant using symmetry relation
   *
   * \note
   * This function is the most computationally expensive part of the algorithm.
   * For the parameter sets used here, this function is responsible for approximately
   * 90% of all floating point and memory operations.
   */
  void compute_dbidrj();

  /*!
   * \brief Compute dEidRj
   *
   * \param dedr
   */
  void compute_deidrj(double *const dedr);

  /*!
   * \brief The switching function
   *
   * The switching function ensures that the contribution of each neighbor
   * atom goes smoothly to zero at Rcut.
   *
   * \param r Neighbor atom distnace to the central atom
   * \param rcut Cutoff distance
   *
   * \return double Computed weight from the switching function
   */
  double compute_sfac(double const r, double const rcut);

  /*!
   * \brief The switching function derivative
   *
   * \param r Neighbor atom distnace to the central atom
   * \param rcut Cutoff distance
   *
   * \return double Computed weight from the switching function derivative
   */
  double compute_dsfac(double const r, double const rcut);

  /*!
   * \brief Make sure arrays are big enough, in case the new size is bigger than
   * the allocated sign, and update the maximum size
   *
   * \param newnmax The new size
   */
  void grow_rij(int const newnmax);

private:
  /*!
   * \brief For any positive number n, returns it's factorial
   *
   * \param n A positive integer number
   *
   * \return double Factorial of a positive integer number \c n
   *
   * \note
   * Overflow happens if n > 170
   */
  inline double factorial(int const n);

  /*!
   * \brief Create and allocate data arrays
   *
   */
  void create_twojmax_arrays();

  /*!
   * \brief Assign Clebsch-Gordan coefficients using the quasi-binomial formula VMK 8.2.1(3)
   *
   */
  void init_clebsch_gordan();

  /*!
   * \brief Pre-compute table of sqrt[p/m2], p, q = 1, twojmax
   * the p = 0, q = 0 entries are allocated and skipped for convenience.
   *
   */
  void init_rootpqarray();

  /*!
   * \brief Set the ulisttot arrays to zero
   *
   */
  void zero_uarraytot();

  /*!
   * \brief Add the weight for the central atom to the start of each block
   *
   * \param wself_in Dimensionless weight for the central atom
   */
  void addself_uarraytot(double const wself_in);

  /*!
   * \brief Add Wigner U-functions for one neighbor to the total
   *
   * \param r Neighbor atom distnace to the central atom
   * \param wj_in Dimensionless weights that are chosen to distinguish atoms of different types
   * \param rcut Cutoff distance
   * \param neigh_j Neighbor index
   */
  void add_uarraytot(double const r, double const wj_in, double const rcut, int const neigh_j);

  /*!
   * \brief Compute Wigner U-functions for one neighbor
   *
   * \param dx
   * \param dy
   * \param dz
   * \param z0
   * \param r Neighbor atom distnace to the central atom
   * \param neigh_j Neighbor index
   *
   * Elements of the three-dimensional finite rotation matrix, which are written in terms of
   * Euler angles are called Wigner’s D-functions
   *
   */
  void compute_uarray(double const dx, double const dy, double const dz, double const z0, double const r, int const neigh_j);

  /*!
   * \brief The function delta given by VMK Eq. 8.2(1)
   *
   * \param j1 Non-negative integer value
   * \param j2 Non-negative integer value
   * \param j Half-integer value
   *
   * \return double deltacg
   */
  double deltacg(int const j1, int const j2, int const j);

  /*!
   * \brief Compute the number of cpefficients
   *
   * \return int Number of cpefficients
   */
  int compute_ncoeff();

  /*!
   * \brief Compute derivatives of Wigner U-functions for one neighbor \sa compute_uarray()
   *
   * \param dx
   * \param dy
   * \param dz
   * \param z0
   * \param r Neighbor atom distnace to the central atom
   * \param dz0dr
   * \param wj_in Dimensionless weights that are chosen to distinguish atoms of different types
   * \param rcut Cutoff distance
   * \param neigh_j Neighbor index
   */
  void compute_duarray(double const dx, double const dy, double const dz,
                       double const z0, double const r, double const dz0dr,
                       double const wj_in, double const rcut, int const neigh_j);

public:
  /*! Band limit for bispectrum components (non-negative integer) */
  int twojmax;

  /*! Number of coefficients*/
  int ncoeff;

  /*! Vector of displacements between atom I and it's neighbors */
  Array2D<double> rij;

  /*! Indices of neighbors of atom I within cutoff */
  std::vector<int> inside;

  /*!
   * \brief Weights for neighbors of atom I within cutoff.
   *
   * Dimensionless weights that are chosen to distinguish atoms of different types
   */
  std::vector<double> wj;

  /*! Cutoffs for neighbors of atom I within cutoff */
  std::vector<double> rcutij;

  /*! */
  std::vector<double> blist;

  /*! */
  Array2D<double> dblist;

private:
  /*! Parameter in distance to angle conversion (distance units) */
  double rmin0;

  /*! Parameter in distance to angle conversion (0 < rcutfac < 1) */
  double rfac0;

  /*!
   * \brief Flag for switching function
   *
   * switchflag value = 0 or 1
   * 0 = do not use switching function
   * 1 = use switching function
   *
   * \note
   * The switching function ensures that the contribution of
   * each neighbor atom goes smoothly to zero at Rcut.
   */
  int switchflag;

  /*!
   * \brief This flag is set to 1 if bzero subtracted from barray
   *
   * bzeroflag value = 0 or 1
   * 0 = do not subtract B0
   * 1 = subtract B0
   */
  int bzeroflag;

  /*!
   * \brief Dimensionless weight for the central atom.
   *
   * \note
   * The central atom is arbitrarily assigned a unit weight.
   */
  double wself;

  /*! Index counter */
  int idxcg_max;

  /*! Index counter */
  int idxu_max;

  /*! Index counter */
  int idxz_max;

  /*! Index counter */
  int idxb_max;

  // data for bispectrum coefficients
  std::vector<SNA_ZINDICES> idxz;

  /*! */
  std::vector<SNA_BINDICES> idxb;

  /*! Array of B values for isolated atoms */
  std::vector<double> bzero;

  /*! The Clebsch–Gordan coupling coefficients */
  std::vector<double> cglist;

  /*! U total list arrays */
  std::vector<double> ulisttot_r;
  std::vector<double> ulisttot_i;

  /*! U total list maximum index */
  std::vector<int> idxu_block;

  /*! Coefficient list for the angular part of the density function */
  std::vector<double> ylist_r;
  std::vector<double> ylist_i;

  /*! Z list arrays */
  std::vector<double> zlist_r;
  std::vector<double> zlist_i;

  /*! U list arrays */
  Array2D<double> ulist_r_ij;
  Array2D<double> ulist_i_ij;

  /*! Square root of p, and q indices */
  Array2D<double> rootpqarray;

  /*! Derivatives of u list arrays */
  Array2D<double> dulist_r;
  Array2D<double> dulist_i;

  /*! Index container */
  Array3D<int> idxcg_block;

  /*! Index container */
  Array3D<int> idxz_block;

  /*! Index container */
  Array3D<int> idxb_block;
};

#endif // SNA_HPP
