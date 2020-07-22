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
// TABLE.hpp
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
//        `lammps/src/pair_table.h` and amended and updated by
//        Yaser Afshar
//


#ifndef TABLE_HPP
#define TABLE_HPP

#include "helper.hpp"

#include <cstdio>

#include <vector>

/*!
 * \brief Method of interpolation
 *
 */
enum class TABLESTYLE
{
  /*! Find the nearest value to the distance R */
  LOOKUP,
  /*! Find the 2 surrounding values to the distance R */
  LINEAR,
  /*! Find the appropriate set of spline coefficients to evaluate a cubic polynomial at the distance R */
  SPLINE,
  /*! Use a fast bit-mapping technique and a linear interpolation to find the nearest value to the distance R */
  BITMAP
};

/*!
 * \brief Tabulated file distance style
 *
 *
 * \note
 * To match the n number of points in the tabulated file with the requested
 * number of points in the input parameter file, we do a preliminary
 * interpolation by creating splines using the n number of tabulated values
 * in the tabulated file as nodal points.
 * Then we interpolate the energy and force values to the requested number of
 * points in the input parameter file. The resulting table would have the
 * requested resolution (number of points).
 *
 * \note
 * If you use \c R or \c RSQ, the tabulated distance values in the file are
 * effectively ignored.
 * If parameters \c R and \c RSQ are followed by 2 values rlo and rhi, the
 * distance associated with each energy and force value is computed from
 * these 2 values (at high accuracy).
 * For \c R, distances uniformly spaced between rlo and rhi are computed;
 * for \c RSQ, squared distances uniformly spaced between rlo*rlo and rhi*rhi
 * are computed.
 *
 * In the \c BITMAP format, the number of parameter in the file must be equal
 * to \f$ 2^n \f$ where \c n is the number of points requested in the input
 * parameter file. The entire table extent as specified in the file must be
 * used.
 */
enum class TABLEDISTANCESTYLE
{
  /*! The distances in each line of the table are used as-is to perform spline interpolation. */
  NONE,
  /*! R separation distance */
  RLINEAR,
  /*! R separation distance squared (default value)*/
  RSQ,
  /*!  BITMAP format */
  BITMAP
};

/*!
 * \brief Table style information class
 *
 */
struct TABLE_INFO
{
  /*! Style of the table */
  TABLESTYLE tableStyle;

  /*! Length of the table requested for interpolation */
  int tableLength;
};

/*!
 * \brief New union declaration
 *
 */
union union_int_float_t {
  int i;
  float f;
};


/*!
 * \brief TABLE class
 *
 * TABLE class creates interpolation tables from potential energy and force
 * values listed in a file(s) as a function of distance.
 *
 */
class TABLE
{
public:
  /*!
   * \brief Construct a new TABLE object
   *
   */
  TABLE();

  /*!
   * \brief Construct a new TABLE object
   *
   * \param table_style Style of the table
   * \param table_length Length of the table requested for interpolation
   */
  TABLE(TABLESTYLE const &table_style, int const table_length);

  /*!
   * \brief Destroy the TABLE object
   *
   */
  ~TABLE();

  /*!
   * \brief Move constructor, Construct a new TABLE object
   *
   * \param other TABLE object
   */
  TABLE(TABLE &&other);

  /*!
   * \brief Move assignment operator
   *
   * \param other TABLE object
   * \return TABLE&
   */
  TABLE &operator=(TABLE &&other);

private:
  TABLE(TABLE const &) = delete;
  TABLE &operator=(TABLE const &) = delete;

public:
  /*!
   * \brief Read a table section from a tabulated potential file and
   * sets the member values which are:
   * ninput,rfile,efile,ffile,rflag,rlo,rhi,fpflag,fplo,fphi,ntablebits
   *
   * \param filePointers Pointer to the current open file
   * \param keyword The keyword specifies a section of the file.
   *
   * \return int 0|false if everything goes well and 1|true if it fails
   */
  int read_table(std::FILE *const filePointer, char const *keyword);

  /*!
   * \brief compute r,e,f vectors from splined values
   *
   */
  int compute_table();

  /*!
   * \brief Build spline representation of e,f over entire range of read-in
   * table this function sets these values in TABLE: e2file,f2file
   *
   */
  void spline_table();

protected:
  /*!
   * \brief Get the Next Data Line
   *
   * \param filePtr Pointer to the current open file
   * \param nextLinePtr Pointer to the line
   * \param maxSize Maximum size of the line
   * \param endOfFileFlag Flag to indicate if we reached the end of the file
   */
  void GetNextDataLine(std::FILE *const filePtr,
                       char *nextLinePtr,
                       int const maxSize,
                       int *endOfFileFlag);

  /*!
   * \brief Extract attributes from parameter line in table section format
   * of line:
   * N value R/RSQ/BITMAP lo hi FPRIME fplo fphi
   * Where, N is required, other params are optional
   *
   * \return int 0|false if everything goes well and 1|true if it fails
   */
  int param_extract(char *line);

public:
  /*! Style of the table */
  TABLESTYLE tableStyle;

  /*! Length of the table requested for interpolation */
  int tableLength;

  /*! Number of input points in the tabulated potential file */
  int ninput;

  /*! The distance style in the tabulated potential file*/
  TABLEDISTANCESTYLE rflag;

  /*! */
  int fpflag;
  /*! */
  int match;
  /*! */
  int ntablebits;
  /*! */
  int nshiftbits;
  /*! */
  int nmask;

  /*! */
  double rlo;
  /*! */
  double rhi;
  /*! */
  double fplo;
  /*! */
  double fphi;
  /*! */
  double cut;

  /*! */
  std::vector<double> rfile;
  /*! */
  std::vector<double> efile;
  /*! */
  std::vector<double> ffile;
  /*! */
  std::vector<double> e2file;
  /*! */
  std::vector<double> f2file;

  /*! */
  double innersq;
  /*! */
  double delta;
  /*! */
  double invdelta;
  /*! */
  double deltasq6;

  /*! */
  std::vector<double> rsq;
  /*! */
  std::vector<double> drsq;
  /*! */
  std::vector<double> e;
  /*! */
  std::vector<double> de;
  /*! */
  std::vector<double> f;
  /*! */
  std::vector<double> df;
  /*! */
  std::vector<double> e2;
  /*! */
  std::vector<double> f2;
};

#endif // TABLE_HPP
