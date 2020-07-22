/* ----------------------------------------------------------------------
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
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

//
// TABLE.cpp
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
//        `lammps/src/pair_table.cpp` and amended and updated by
//        Yaser Afshar
//


#include "TABLE.hpp"

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <climits>
#include <cfloat>

#include <string>
#include <utility>
#include <algorithm>
#include <vector>

#ifdef MAXLINE
#undef MAXLINE
#endif
#define MAXLINE 1024

#ifdef EPSILONR
#undef EPSILONR
#endif
#define EPSILONR 1.0e-6

/*!
 * \brief Construct a cubic spline given a set of x and y values
 *
 * \param x Input array of length n
 * \param y Input function value array of length n in order
 * \param n Length of the array (last element index)
 * \param yp1 The first derivative of the interpolating function at point 1
 * \param ypn The first derivative of the interpolating function at point n
 * \param y2 Output function value array of length n
 */
static void spline(std::vector<double> const &x,
                   std::vector<double> const &y,
                   int const n,
                   double const yp1,
                   double const ypn,
                   std::vector<double> &y2);

/*!
 * \brief Returns a cubic spline interpolated value y at given x.
 *
 * \param xa Input array of length n
 * \param ya Input function value array of length n in order
 * \param y2a Output function value array of length n from spline \sa spline
 * \param n Length of the array (last element index)
 * \param x A value of x
 *
 * \return double Interpolated value y
 */
static double splint(std::vector<double> const &xa,
                     std::vector<double> const &ya,
                     std::vector<double> const &y2a,
                     int const n,
                     double const x);

/*!
 * \brief Define bitmap parameters based on inner and outer cutoffs
 *
 * \param inner inner cutoff
 * \param outer outer cutoff
 * \param ntablebits number of table bits based on the length in table file
 * \param masklo
 * \param maskhi
 * \param nmask
 * \param nshiftbits
 *
 * \return int 0|false if everything goes well and 1|true if it fails
 */
static int init_bitmap(double const inner,
                       double const outer,
                       int const ntablebits,
                       int &masklo,
                       int &maskhi,
                       int &nmask,
                       int &nshiftbits);

TABLE::TABLE() : tableLength(0),
                 ninput(0),
                 rflag(TABLEDISTANCESTYLE::NONE),
                 fpflag(0),
                 match(0),
                 ntablebits(0),
                 nshiftbits(0),
                 nmask(0),
                 rlo(0.0),
                 rhi(0.0),
                 fplo(0.0),
                 fphi(0.0),
                 cut(0.0),
                 innersq(0.0),
                 delta(0.0),
                 invdelta(0.0),
                 deltasq6(0.0)
{
}

TABLE::TABLE(TABLESTYLE const &table_style, int const table_length) : tableStyle(table_style),
                                                                      tableLength(table_length),
                                                                      ninput(0),
                                                                      rflag(TABLEDISTANCESTYLE::NONE),
                                                                      fpflag(0),
                                                                      match(0),
                                                                      ntablebits(0),
                                                                      nshiftbits(0),
                                                                      nmask(0),
                                                                      rlo(0.0),
                                                                      rhi(0.0),
                                                                      fplo(0.0),
                                                                      fphi(0.0),
                                                                      cut(0.0),
                                                                      innersq(0.0),
                                                                      delta(0.0),
                                                                      invdelta(0.0),
                                                                      deltasq6(0.0)
{
}

TABLE::~TABLE() {}

TABLE::TABLE(TABLE &&other)
{
  tableStyle = other.tableStyle;
  tableLength = other.tableLength;
  ninput = other.ninput;
  rflag = other.rflag;
  fpflag = other.fpflag;
  match = other.match;
  ntablebits = other.ntablebits;
  nshiftbits = other.nshiftbits;
  nmask = other.nmask;
  rlo = other.rlo;
  rhi = other.rhi;
  fplo = other.fplo;
  fphi = other.fphi;
  cut = other.cut;
  rfile = std::move(other.rfile);
  efile = std::move(other.efile);
  ffile = std::move(other.ffile);
  e2file = std::move(other.e2file);
  f2file = std::move(other.f2file);
  innersq = other.innersq;
  delta = other.delta;
  invdelta = other.invdelta;
  deltasq6 = other.deltasq6;
  rsq = std::move(other.rsq);
  drsq = std::move(other.drsq);
  e = std::move(other.e);
  de = std::move(other.de);
  f = std::move(other.f);
  df = std::move(other.df);
  e2 = std::move(other.e2);
  f2 = std::move(other.f2);
}

TABLE &TABLE::operator=(TABLE &&other)
{
  tableStyle = other.tableStyle;
  tableLength = other.tableLength;
  ninput = other.ninput;
  rflag = other.rflag;
  fpflag = other.fpflag;
  match = other.match;
  ntablebits = other.ntablebits;
  nshiftbits = other.nshiftbits;
  nmask = other.nmask;
  rlo = other.rlo;
  rhi = other.rhi;
  fplo = other.fplo;
  fphi = other.fphi;
  cut = other.cut;
  rfile = std::move(other.rfile);
  efile = std::move(other.efile);
  ffile = std::move(other.ffile);
  e2file = std::move(other.e2file);
  f2file = std::move(other.f2file);
  innersq = other.innersq;
  delta = other.delta;
  invdelta = other.invdelta;
  deltasq6 = other.deltasq6;
  rsq = std::move(other.rsq);
  drsq = std::move(other.drsq);
  e = std::move(other.e);
  de = std::move(other.de);
  f = std::move(other.f);
  df = std::move(other.df);
  e2 = std::move(other.e2);
  f2 = std::move(other.f2);

  return *this;
}

void TABLE::GetNextDataLine(std::FILE *const filePtr,
                            char *nextLinePtr,
                            int const maxSize,
                            int *endOfFileFlag)
{
  *endOfFileFlag = 0;
  do
  {
    if (!std::fgets(nextLinePtr, maxSize, filePtr))
    {
      *endOfFileFlag = 1;
      break;
    }

    while (nextLinePtr[0] == ' ' || nextLinePtr[0] == '\t' || nextLinePtr[0] == '\n' || nextLinePtr[0] == '\r')
      nextLinePtr++;

  } while ((std::strncmp("#", nextLinePtr, 1) == 0) || (std::strlen(nextLinePtr) == 0));

  // remove comments starting with `#' in a line
  char *pch = std::strchr(nextLinePtr, '#');
  if (pch)
  {
    *pch = '\0';
  }
}

int TABLE::read_table(std::FILE *const filePointers, char const *keyword)
{
  // Rewind the file in case it is in the middle of the file
  if (std::fseek(filePointers, 0, SEEK_SET))
  {
    HELPER_LOG_ERROR("Failed to rewind the Tabulated file.\n");
    return true;
  }

  char nextLine[MAXLINE];

  int endOfFileFlag;

  // Loop until section found with matching keyword
  while (1)
  {
    GetNextDataLine(filePointers, nextLine, MAXLINE, &endOfFileFlag);
    if (endOfFileFlag)
    {
      HELPER_LOG_ERROR("Reached the end of the Tabulated file.\n"
                       "Unable to find the keyword = " +
                       std::string(keyword) +
                       "\n");
      return true;
    }

    char *word = std::strtok(nextLine, " \t\n\r");
    if (std::strcmp(word, keyword) == 0)
    {
      // matching keyword
      break;
    }
  } // End of loop to find the matching keyword

  {
    // Read args on the 2nd line of section
    // allocate table arrays for file values
    GetNextDataLine(filePointers, nextLine, MAXLINE, &endOfFileFlag);
    if (endOfFileFlag)
    {
      HELPER_LOG_ERROR("Reached the end of Tabulated file.\n");
      return true;
    }

    if (param_extract(nextLine))
    {
      // Failed to extract the parameters
      return true;
    }

    rfile.resize(ninput);
    efile.resize(ninput);
    ffile.resize(ninput);

    // setup bitmap parameters for table to read in
    ntablebits = 0;

    int masklo;
    int maskhi;
    int nmask;
    int nshiftbits;

    if (rflag == TABLEDISTANCESTYLE::BITMAP)
    {
      while (1 << ntablebits < ninput)
      {
        ++ntablebits;
      }

      if (1 << ntablebits != ninput)
      {
        HELPER_LOG_ERROR("Bitmapped table is incorrect length in table file.\n");
        return true;
      }

      if (init_bitmap(rlo, rhi, ntablebits, masklo, maskhi, nmask, nshiftbits))
      {
        return true;
      }
    } // if BITMAP

    // read r,e,f table values from file
    // if rflag set, compute r
    // if rflag not set, use r from file
    int itmp;
    double rold;
    double rnew;

    union_int_float_t rsq_lookup;

    int rerror = 0;
    int cerror = 0;

    for (int i = 0; i < ninput; ++i)
    {
      GetNextDataLine(filePointers, nextLine, MAXLINE, &endOfFileFlag);
      if (endOfFileFlag)
      {
        HELPER_LOG_ERROR("Reached the end of the Tabulated file.\n");
        return true;
      }

      if (4 != std::sscanf(nextLine, "%d %lg %lg %lg", &itmp, &rold, &efile[i], &ffile[i]))
      {
        ++cerror;
      }

      if (rflag == TABLEDISTANCESTYLE::NONE)
      {
        rnew = rold;
      }
      else if (rflag == TABLEDISTANCESTYLE::RLINEAR)
      {
        rnew = rlo + (rhi - rlo) * i / (ninput - 1);
      }
      else if (rflag == TABLEDISTANCESTYLE::RSQ)
      {
        rnew = rlo * rlo + (rhi * rhi - rlo * rlo) * i / (ninput - 1);
        rnew = std::sqrt(rnew);
      }
      else if (rflag == TABLEDISTANCESTYLE::BITMAP)
      {
        rsq_lookup.i = i << nshiftbits;
        rsq_lookup.i |= masklo;
        if (rsq_lookup.f < rlo * rlo)
        {
          rsq_lookup.i = i << nshiftbits;
          rsq_lookup.i |= maskhi;
        }
        float const rsq_lookup_f = std::sqrt(rsq_lookup.f);
        rnew = static_cast<double>(rsq_lookup_f);
      }

      if (std::fabs(rnew - rold) / rold > EPSILONR)
      {
        ++rerror;
      }

      rfile[i] = rnew;
    } // End of loop through the inputs

    // warn if force != dE/dr at any point that is not an inflection point
    // check via secant approximation to dE/dr
    // skip two end points since do not have surrounding secants
    // inflection point is where curvature changes sign

    double r;
    double e;
    double f;
    double rprev;
    double rnext;
    double eprev;
    double enext;
    double fleft;
    double fright;

    int ferror = 0;

    // bitmapped tables do not follow regular ordering, so we cannot check them here
    if (rflag != TABLEDISTANCESTYLE::BITMAP)
    {
      for (int i = 1; i < ninput - 1; ++i)
      {
        r = rfile[i];
        rprev = rfile[i - 1];
        rnext = rfile[i + 1];

        e = efile[i];
        eprev = efile[i - 1];
        enext = efile[i + 1];

        f = ffile[i];
        fleft = -(e - eprev) / (r - rprev);
        fright = -(enext - e) / (rnext - r);

        if (f < fleft && f < fright)
        {
          ++ferror;
        }

        if (f > fleft && f > fright)
        {
          ++ferror;
        }
      }
    }

    if (ferror)
    {
      HELPER_LOG_WARNING(std::to_string(ferror) +
                         " of " +
                         std::to_string(ninput) +
                         " force values in table are inconsistent with -dE/dr.\n" +
                         " Should only be flagged at inflection points.\n");
    }

    // warn if re-computed distance values differ from file values
    if (rerror)
    {
      HELPER_LOG_WARNING(std::to_string(rerror) +
                         " of " +
                         std::to_string(ninput) +
                         " distance values in table with relative error.\n" +
                         " over " +
                         std::to_string(EPSILONR) +
                         " to re-computed values.\n");
    }

    // warn if data was read incompletely, e.g. columns were missing
    if (cerror)
    {
      HELPER_LOG_WARNING(std::to_string(cerror) +
                         " of " +
                         std::to_string(ninput) +
                         " lines in table were incomplete.\n" +
                         " or could not be parsed completely.\n");
    }
  }

  // everything is good
  return false;
}

int TABLE::param_extract(char *line)
{
  ninput = 0;
  rflag = TABLEDISTANCESTYLE::NONE;
  fpflag = 0;

  char *word = std::strtok(line, " \t\n\r\f");
  while (word)
  {
    if (std::strcmp(word, "N") == 0)
    {
      word = std::strtok(NULL, " \t\n\r\f");
      ninput = std::atoi(word);
    }
    else if (std::strcmp(word, "R") == 0 ||
             std::strcmp(word, "RSQ") == 0 ||
             std::strcmp(word, "BITMAP") == 0)
    {
      if (std::strcmp(word, "R") == 0)
      {
        rflag = TABLEDISTANCESTYLE::RLINEAR;
      }
      else if (std::strcmp(word, "RSQ") == 0)
      {
        rflag = TABLEDISTANCESTYLE::RSQ;
      }
      else if (std::strcmp(word, "BITMAP") == 0)
      {
        rflag = TABLEDISTANCESTYLE::BITMAP;
      }

      word = std::strtok(NULL, " \t\n\r\f");
      rlo = std::atof(word);

      word = std::strtok(NULL, " \t\n\r\f");
      rhi = std::atof(word);
    }
    // The parameter “FPRIME” is followed by 2 values fplo and fphi which
    // are the derivative of the force at the innermost and outermost
    // distances listed in the table.
    //  If not specified by the “FPRIME” parameter, they are estimated
    // (less accurately) by the first 2 and last 2 force values in the
    // table.
    else if (std::strcmp(word, "FPRIME") == 0)
    {
      fpflag = 1;
      word = std::strtok(NULL, " \t\n\r\f");
      fplo = std::atof(word);

      word = std::strtok(NULL, " \t\n\r\f");
      fphi = std::atof(word);
    }
    else
    {
      HELPER_LOG_ERROR("Invalid keyword (" +
                       std::string(word) +
                       ") in the tabulated file!\n");
      return true;
    }

    word = std::strtok(NULL, " \t\n\r\f");
  }

  if (ninput)
  {
    return false;
  }

  HELPER_LOG_ERROR("There is no N indicating the number of entries in the tabulated file!\n");
  return true;
}

int TABLE::compute_table()
{
  int const tlm1 = tableLength - 1;

  // inner = inner table bound
  double const inner = (rflag == TABLEDISTANCESTYLE::NONE) ? rfile[0] : rlo;

  // cut = outer table bound
  // delta = table spacing in rsq for N-1 bins
  innersq = inner * inner;
  delta = (cut * cut - innersq) / tlm1;
  invdelta = 1.0 / delta;

  // direct lookup tables
  // N-1 evenly spaced bins in rsq from inner to cut
  // e,f = value at midpt of bin
  // e,f are N-1 in length since store 1 value at bin midpt
  // f is converted to f/r when stored in f[i]
  // e,f are never a match to read-in values, always computed via spline interp

  if (tableStyle == TABLESTYLE::LOOKUP)
  {
    e.resize(tlm1);
    f.resize(tlm1);

    for (int i = 0; i < tlm1; i++)
    {
      double const _rsq = innersq + (i + 0.5) * delta;
      double const r = std::sqrt(_rsq);
      e[i] = splint(rfile, efile, e2file, ninput, r);
      f[i] = splint(rfile, ffile, f2file, ninput, r) / r;
    }
    return false;
  }

  // linear tables
  // N-1 evenly spaced bins in rsq from inner to cut
  // rsq,e,f = value at lower edge of bin
  // de,df values = delta from lower edge to upper edge of bin
  // rsq,e,f are N in length so de,df arrays can compute difference
  // f is converted to f/r when stored in f[i]
  // e,f can match read-in values, else compute via spline interp

  if (tableStyle == TABLESTYLE::LINEAR)
  {
    rsq.resize(tableLength);
    e.resize(tableLength);
    f.resize(tableLength);
    de.resize(tlm1);
    df.resize(tlm1);

    for (int i = 0; i < tableLength; ++i)
    {
      double const _rsq = innersq + i * delta;
      double const r = std::sqrt(_rsq);
      rsq[i] = _rsq;
      if (match)
      {
        e[i] = efile[i];
        f[i] = ffile[i] / r;
      }
      else
      {
        e[i] = splint(rfile, efile, e2file, ninput, r);
        f[i] = splint(rfile, ffile, f2file, ninput, r) / r;
      }
    }

    for (int i = 0; i < tlm1; ++i)
    {
      de[i] = e[i + 1] - e[i];
      df[i] = f[i + 1] - f[i];
    }

    return false;
  }

  // cubic spline tables
  // N-1 evenly spaced bins in rsq from inner to cut
  // rsq,e,f = value at lower edge of bin
  // e2,f2 = spline coefficient for each bin
  // rsq,e,f,e2,f2 are N in length so have N-1 spline bins
  // f is converted to f/r after e is splined
  // e,f can match read-in values, else compute via spline interp

  if (tableStyle == TABLESTYLE::SPLINE)
  {
    rsq.resize(tableLength);
    e.resize(tableLength);
    f.resize(tableLength);
    e2.resize(tableLength);
    f2.resize(tableLength);

    deltasq6 = delta * delta / 6.0;

    for (int i = 0; i < tableLength; ++i)
    {
      double const _rsq = innersq + i * delta;
      double const r = std::sqrt(_rsq);
      rsq[i] = _rsq;
      if (match)
      {
        e[i] = efile[i];
        f[i] = ffile[i] / r;
      }
      else
      {
        e[i] = splint(rfile, efile, e2file, ninput, r);
        f[i] = splint(rfile, ffile, f2file, ninput, r);
      }
    }

    // ep0,epn = dh/dg at inner and at cut
    // h(r) = e(r) and g(r) = r^2
    // dh/dg = (de/dr) / 2r = -f/2r

    double const ep0 = -f[0] / (2.0 * std::sqrt(innersq));
    double const epn = -f[tlm1] / (2.0 * cut);
    spline(rsq, e, tableLength, ep0, epn, e2);

    // fp0,fpn = dh/dg at inner and at cut
    // h(r) = f(r)/r and g(r) = r^2
    // dh/dg = (1/r df/dr - f/r^2) / 2r
    // dh/dg in secant approx = (f(r2)/r2 - f(r1)/r1) / (g(r2) - g(r1))

    double fp0, fpn;
    double secant_factor = 0.1;
    if (fpflag)
    {
      fp0 = (fplo / std::sqrt(innersq) - f[0] / innersq) / (2.0 * std::sqrt(innersq));
    }
    else
    {
      double const rsq1 = innersq;
      double const rsq2 = rsq1 + secant_factor * delta;
      fp0 = (splint(rfile, ffile, f2file, ninput, std::sqrt(rsq2)) / std::sqrt(rsq2) -
             f[0] / std::sqrt(rsq1)) /
            (secant_factor * delta);
    }

    if (fpflag && cut == rfile[ninput - 1])
    {
      fpn = (fphi / cut - f[tlm1] / (cut * cut)) / (2.0 * cut);
    }
    else
    {
      double const rsq2 = cut * cut;
      double const rsq1 = rsq2 - secant_factor * delta;
      fpn = (f[tlm1] / std::sqrt(rsq2) - splint(rfile, ffile, f2file, ninput, std::sqrt(rsq1)) / std::sqrt(rsq1)) /
            (secant_factor * delta);
    }

    for (int i = 0; i < tableLength; ++i)
    {
      f[i] /= std::sqrt(rsq[i]);
    }

    spline(rsq, f, tableLength, fp0, fpn, f2);

    return false;
  }

  // bitmapped linear tables
  // 2^N bins from inner to cut, spaced in bitmapped manner
  // f is converted to f/r when stored in f[i]
  // e,f can match read-in values, else compute via spline interp

  if (tableStyle == TABLESTYLE::BITMAP)
  {
    // linear lookup tables of length ntable = 2^n
    // stored value = value at lower edge of bin
    int masklo;
    int maskhi;

    if (init_bitmap(inner, cut, tableLength, masklo, maskhi, nmask, nshiftbits))
    {
      return true;
    }

    int const ntable = 1 << tableLength;
    int const ntablem1 = ntable - 1;

    rsq.resize(ntable);
    drsq.resize(ntable);
    e.resize(ntable);
    f.resize(ntable);
    de.resize(ntable);
    df.resize(ntable);

    union_int_float_t rsq_lookup;
    union_int_float_t minrsq_lookup;

    minrsq_lookup.i = 0 << nshiftbits;
    minrsq_lookup.i |= maskhi;

    for (int i = 0; i < ntable; ++i)
    {
      rsq_lookup.i = i << nshiftbits;
      rsq_lookup.i |= masklo;
      if (rsq_lookup.f < innersq)
      {
        rsq_lookup.i = i << nshiftbits;
        rsq_lookup.i |= maskhi;
      }

      float const rsq_lookup_f = std::sqrt(rsq_lookup.f);
      double const r = static_cast<double>(rsq_lookup_f);

      rsq[i] = rsq_lookup.f;

      if (match)
      {
        e[i] = efile[i];
        f[i] = ffile[i] / r;
      }
      else
      {
        e[i] = splint(rfile, efile, e2file, ninput, r);
        f[i] = splint(rfile, ffile, f2file, ninput, r) / r;
      }

      minrsq_lookup.f = std::min(minrsq_lookup.f, rsq_lookup.f);
    }

    innersq = minrsq_lookup.f;

    for (int i = 0; i < ntablem1; ++i)
    {
      de[i] = e[i + 1] - e[i];
      df[i] = f[i + 1] - f[i];

      drsq[i] = 1.0 / (rsq[i + 1] - rsq[i]);
    }

    // get the delta values for the last table entries
    // tables are connected periodically between 0 and ntablem1

    de[ntablem1] = e[0] - e[ntablem1];
    df[ntablem1] = f[0] - f[ntablem1];

    drsq[ntablem1] = 1.0 / (rsq[0] - rsq[ntablem1]);

    // get the correct delta values at itablemax
    // smallest r is in bin itablemin
    // largest r is in bin itablemax, which is itablemin-1,
    //   or ntablem1 if itablemin=0

    // deltas at itablemax only needed if corresponding rsq < cut*cut
    // if so, compute deltas between rsq and cut*cut
    //   if match, data at cut*cut is unavailable, so we'll take
    //   deltas at itablemax-1 as a good approximation

    int itablemin = minrsq_lookup.i & nmask;
    itablemin >>= nshiftbits;

    int itablemax = itablemin - 1;

    if (itablemin == 0)
    {
      itablemax = ntablem1;
    }

    int itablemaxm1 = itablemax - 1;
    if (itablemax == 0)
    {
      itablemaxm1 = ntablem1;
    }

    rsq_lookup.i = itablemax << nshiftbits;
    rsq_lookup.i |= maskhi;

    if (rsq_lookup.f < cut * cut)
    {
      if (match)
      {
        de[itablemax] = de[itablemaxm1];
        df[itablemax] = df[itablemaxm1];

        drsq[itablemax] = drsq[itablemaxm1];
      }
      else
      {
        rsq_lookup.f = cut * cut;
        float const rsq_lookup_f = std::sqrt(rsq_lookup.f);
        double const r = static_cast<double>(rsq_lookup_f);

        double const e_tmp = splint(rfile, efile, e2file, ninput, r);
        double const f_tmp = splint(rfile, ffile, f2file, ninput, r) / r;

        de[itablemax] = e_tmp - e[itablemax];
        df[itablemax] = f_tmp - f[itablemax];

        drsq[itablemax] = 1.0 / (rsq_lookup.f - rsq[itablemax]);
      }
    }

    return false;
  }

  // TABLE style did not match
  return true;
}

void TABLE::spline_table()
{
  e2file.resize(ninput);
  f2file.resize(ninput);

  double const ep0 = -ffile[0];
  double const epn = -ffile[ninput - 1];

  spline(rfile, efile, ninput, ep0, epn, e2file);

  if (!fpflag)
  {
    fplo = (ffile[1] - ffile[0]) / (rfile[1] - rfile[0]);
    fphi = (ffile[ninput - 1] - ffile[ninput - 2]) /
           (rfile[ninput - 1] - rfile[ninput - 2]);
  }

  double const fp0 = fplo;
  double const fpn = fphi;

  spline(rfile, ffile, ninput, fp0, fpn, f2file);
}

void spline(std::vector<double> const &x,
            std::vector<double> const &y,
            int const n,
            double const yp1,
            double const ypn,
            std::vector<double> &y2)
{
  std::vector<double> u(n);

  // The largest anticipated value of
  if (yp1 > 0.99e30)
  {
    y2[0] = 0.0;
    u[0] = 0.0;
  }
  else
  {
    y2[0] = -0.5;
    u[0] = (3.0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
  }

  for (int i = 1; i < n - 1; i++)
  {
    double const sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
    double const p = sig * y2[i - 1] + 2.0;
    y2[i] = (sig - 1.0) / p;
    u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
    u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
  }

  double qn, un;

  if (ypn > 0.99e30)
  {
    qn = 0.0;
    un = 0.0;
  }
  else
  {
    qn = 0.5;
    un = (3.0 / (x[n - 1] - x[n - 2])) * (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
  }

  y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);

  for (int k = n - 2; k >= 0; --k)
  {
    y2[k] = y2[k] * y2[k + 1] + u[k];
  }
}

double splint(std::vector<double> const &xa,
              std::vector<double> const &ya,
              std::vector<double> const &y2a,
              int const n,
              double const x)
{
  int klo = 0;
  int khi = n - 1;
  while (khi - klo > 1)
  {
    int const k = (khi + klo) >> 1;
    if (xa[k] > x)
    {
      khi = k;
    }
    else
    {
      klo = k;
    }
  }

  double const h = xa[khi] - xa[klo];
  double const a = (xa[khi] - x) / h;
  double const b = (x - xa[klo]) / h;
  double const y = a * ya[klo] + b * ya[khi] + ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * (h * h) / 6.0;
  return y;
}

int init_bitmap(double const inner,
                double const outer,
                int const ntablebits,
                int &masklo,
                int &maskhi,
                int &nmask,
                int &nshiftbits)
{
  if (sizeof(int) != sizeof(float))
  {
    HELPER_LOG_ERROR("Bitmapped lookup tables require int/float be same size.\n");
    return true;
  }

  if (ntablebits > static_cast<int>(sizeof(float)) * CHAR_BIT)
  {
    HELPER_LOG_ERROR("Too many total bits for bitmapped lookup table.\n");
    return true;
  }

  if (inner >= outer)
  {
    HELPER_LOG_ERROR("TABLE inner cutoff >= outer cutoff.\n");
    return true;
  }

  int nlowermin = 1;
  while (!((std::pow(static_cast<double>(2), static_cast<double>(nlowermin)) <= inner * inner) &&
           (std::pow(static_cast<double>(2), static_cast<double>(nlowermin) + 1.0) > inner * inner)))
  {
    if (std::pow(static_cast<double>(2), static_cast<double>(nlowermin)) <= inner * inner)
    {
      ++nlowermin;
    }
    else
    {
      --nlowermin;
    }
  }

  double const required_range = outer * outer / std::pow(static_cast<double>(2), static_cast<double>(nlowermin));

  int nexpbits = 0;
  double available_range = 2.0;
  while (available_range < required_range)
  {
    ++nexpbits;
    available_range = std::pow(static_cast<double>(2), std::pow(static_cast<double>(2), static_cast<double>(nexpbits)));
  }

  int const nmantbits = ntablebits - nexpbits;

  if (nexpbits > static_cast<int>(sizeof(float)) * CHAR_BIT - FLT_MANT_DIG)
  {
    HELPER_LOG_ERROR("Too many exponent bits for lookup table.\n");
    return true;
  }
  if (nmantbits + 1 > FLT_MANT_DIG)
  {
    HELPER_LOG_ERROR("Too many mantissa bits for lookup table.\n");
    return true;
  }
  if (nmantbits < 3)
  {
    HELPER_LOG_ERROR("Too few bits for lookup table.\n");
    return true;
  }

  nshiftbits = FLT_MANT_DIG - (nmantbits + 1);

  nmask = 1;
  for (int j = 0; j < ntablebits + nshiftbits; j++)
  {
    nmask *= 2;
  }
  --nmask;

  union_int_float_t rsq_lookup;

  rsq_lookup.f = outer * outer;
  maskhi = rsq_lookup.i & ~(nmask);

  rsq_lookup.f = inner * inner;
  masklo = rsq_lookup.i & ~(nmask);

  // everythhing is OK
  return false;
}

#undef EPSILONR
#undef MAXLINE
