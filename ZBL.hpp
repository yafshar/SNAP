/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

//
// ZBL.hpp
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
// Copyright (c) 2020, Regents of the University of Minnesota.
// All rights reserved.
//
// Contributors:
//    Yaser Afshar
//
// Brief: This file is adapted from the LAMMPS software package
//        `lammps/src/pair_zbl.h` and amended and updated by
//        Yaser Afshar
//


#ifndef ZBL_HPP
#define ZBL_HPP

#include "KIM_ModelDriverHeaders.hpp"

#include "helper.hpp"

/*!
 * \brief ZBL pair interaction style
 *
 * It computes the Ziegler-Biersack-Littmark (ZBL) screened nuclear repulsion
 * for describing high-energy collisions between atoms.
 *
 * The ZBL interaction is already smoothed to 0.0 at the cutoff.
 */
class ZBL
{
public:
  /*!
   * \brief Construct a new ZBL object
   *
   */
  ZBL();

  /*!
   * \brief Construct a new ZBL object
   *
   * \param inner Distance where switching function begins
   * \param outer Global cutoff for ZBL interaction
   */
  ZBL(double const inner, double const outer);

  /*!
   * \brief Destroy the ZBL object
   *
   */
  ~ZBL();

  /*!
   * \brief Move constructor, Construct a new ZBL object
   *
   * \param other ZBL object
   */
  ZBL(ZBL &&other);

  /*!
   * \brief Move assignment operator
   *
   * \param other ZBL object
   * \return ZBL&
   */
  ZBL &operator=(ZBL &&other);

private:
  ZBL(ZBL const &) = delete;
  ZBL &operator=(ZBL const &) = delete;

public:
  /*!
   * \brief Allocate members to the correct size
   *
   * \param newSize New size
   */
  void allocate(int const newSize);

  /*!
   * \brief Set the coeff
   *
   * \param iSpecies i atom type
   * \param jSpecies j atom type
   * \param z_iSpecies The nuclear charge of the iSpecies atoms
   * \param z_jSpecies The nuclear charge of the jSpecies atoms
   * \param angstrom 1 angstrom in native units
   * \param qqr2e Conversion of q^2/r to energy
   * \param qelectron 1 electron charge abs() in native units
   */
  void set_coeff(int const iSpecies,
                 int const jSpecies,
                 double const z_iSpecies,
                 double const z_jSpecies,
                 double const angstrom,
                 double const qqr2e,
                 double const qelectron);

  /*!
   * \brief Compute ZBL pair energy
   *
   * \param rij Pair atoms distance
   * \param iSpecies i atom type
   * \param jSpecies j atom type
   * \return double
   */
  double e_zbl(double const rij, int const iSpecies, int const jSpecies);

  /*!
   * \brief Compute ZBL first derivative
   *
   * \param rij Pair atoms distance
   * \param iSpecies i atom type
   * \param jSpecies j atom type
   * \return double
   */
  double dzbldr(double const rij, int const iSpecies, int const jSpecies);

  /*!
   * \brief Compute ZBL second derivative
   *
   * \param rij Pair atoms distance
   * \param iSpecies i atom type
   * \param jSpecies j atom type
   * \return double
   */
  double d2zbldr2(double const rij, int const iSpecies, int const jSpecies);

public:
  /*! Distance where switching function begins */
  double cut_inner;

  /*! Distance square where switching function begins */
  double cut_innersq;

  /*! Global cutoff for ZBL interaction */
  double cut_global;

  /*! Global cutoff square for ZBL interaction */
  double cut_globalsq;

private:
  Array2D<double> d1a;
  Array2D<double> d2a;
  Array2D<double> d3a;
  Array2D<double> d4a;

  Array2D<double> zze;

public:
  Array2D<double> sw1;
  Array2D<double> sw2;
  Array2D<double> sw3;
  Array2D<double> sw4;
  Array2D<double> sw5;
};

namespace ZBLConstants
{
// ZBL constants
static double const pzbl = 0.23;
static double const a0 = 0.46850;
static double const c1 = 0.02817;
static double const c2 = 0.28022;
static double const c3 = 0.50986;
static double const c4 = 0.18175;
static double const d1 = 0.20162;
static double const d2 = 0.40290;
static double const d3 = 0.94229;
static double const d4 = 3.19980;
} // namespace ZBLConstants

#endif // ZBL_HPP
