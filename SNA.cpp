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
   Contributing authors: Aidan Thompson, Christian Trott, SNL
------------------------------------------------------------------------- */

#include "SNA.hpp"

#include <cmath>

#ifdef MY_PI
#undef MY_PI
#endif

#define MY_PI 3.14159265358979323846 // pi

SNA::SNA(double const rfac0_in,
         int const twojmax_in,
         double const rmin0_in,
         int const switchflag_in,
         int const bzeroflag_in) : twojmax(twojmax_in),
                                   rmin0(rmin0_in),
                                   rfac0(rfac0_in),
                                   switchflag(switchflag_in),
                                   wself(1.0),
                                   bzeroflag(bzeroflag_in)
{
  ncoeff = compute_ncoeff();

  build_indexlist();

  create_arrays();

  if (bzeroflag)
  {
    double const www = wself * wself * wself;
    for (int j = 0; j <= twojmax; ++j)
    {
      bzero[j] = www * (j + 1);
    }
  }
}

SNA::~SNA() {}

void SNA::build_indexlist()
{
  // index list for cglist
  int const jdim = twojmax + 1;

  idxcg_block.resize(jdim, jdim, jdim);

  int counter = 0;
  for (int j1 = 0; j1 <= twojmax; ++j1)
  {
    for (int j2 = 0; j2 <= j1; ++j2)
    {
      for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2)
      {
        idxcg_block(j1, j2, j) = counter;
        counter += (j1 + 1) * (j2 + 1);
      }
    }
  }
  idxcg_max = counter;

  // index list for uarray need to include both halves
  idxu_block.resize(jdim, 0);

  counter = 0;
  for (int j = 0; j <= twojmax; ++j)
  {
    idxu_block[j] = counter;
    counter += (j + 1) * (j + 1);
  }
  idxu_max = counter;

  // index list for beta and B
  counter = 0;
  for (int j1 = 0; j1 <= twojmax; ++j1)
  {
    for (int j2 = 0; j2 <= j1; ++j2)
    {
      for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2)
      {
        if (j >= j1)
        {
          ++counter;
        }
      }
    }
  }
  idxb_max = counter;

  idxb.resize(idxb_max);

  // reverse index list for beta and b
  idxb_block.resize(jdim, jdim, jdim);

  counter = 0;
  for (int j1 = 0; j1 <= twojmax; ++j1)
  {
    for (int j2 = 0; j2 <= j1; ++j2)
    {
      for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2)
      {
        if (j >= j1)
        {
          idxb[counter].j1 = j1;
          idxb[counter].j2 = j2;
          idxb[counter].j = j;

          idxb_block(j1, j2, j) = counter++;
        }
      }
    }
  }

  // index list for zlist
  counter = 0;
  for (int j1 = 0; j1 <= twojmax; ++j1)
  {
    for (int j2 = 0; j2 <= j1; ++j2)
    {
      for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2)
      {
        for (int mb = 0; 2 * mb <= j; ++mb)
        {
          counter += j + 1;
        }
      }
    }
  }
  idxz_max = counter;

  idxz.resize(idxz_max);

  idxz_block.resize(jdim, jdim, jdim);

  counter = 0;
  for (int j1 = 0; j1 <= twojmax; ++j1)
  {
    for (int j2 = 0; j2 <= j1; ++j2)
    {
      for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2)
      {
        idxz_block(j1, j2, j) = counter;

        // find right beta[jjb] entry multiply and divide by j+1 factors
        // account for multiplicity of 1, 2, or 3
        for (int mb = 0; 2 * mb <= j; ++mb)
        {
          for (int ma = 0; ma <= j; ++ma, ++counter)
          {
            idxz[counter].j1 = j1;
            idxz[counter].j2 = j2;
            idxz[counter].j = j;

            idxz[counter].ma_min = std::max(0, (2 * ma - j - j2 + j1) / 2);
            idxz[counter].ma_max = (2 * ma - j - (2 * idxz[counter].ma_min - j1) + j2) / 2;

            idxz[counter].na = std::min(j1, (2 * ma - j + j2 + j1) / 2) - idxz[counter].ma_min + 1;

            idxz[counter].mb_min = std::max(0, (2 * mb - j - j2 + j1) / 2);
            idxz[counter].mb_max = (2 * mb - j - (2 * idxz[counter].mb_min - j1) + j2) / 2;

            idxz[counter].nb = std::min(j1, (2 * mb - j + j2 + j1) / 2) - idxz[counter].mb_min + 1;

            // apply to z(j1,j2,j,ma,mb) to unique element of y(j)

            idxz[counter].jju = idxu_block[j] + (j + 1) * mb + ma;
          }
        }
      }
    }
  }
}

void SNA::init()
{
  init_clebsch_gordan();
  init_rootpqarray();
}

void SNA::grow_rij(int const newnmax)
{
  if (newnmax <= (int)rcutij.size())
  {
    return;
  }

  rij.resize(newnmax, 3);
  inside.resize(newnmax);
  wj.resize(newnmax);
  rcutij.resize(newnmax);
  ulist_r_ij.resize(newnmax, idxu_max);
  ulist_i_ij.resize(newnmax, idxu_max);
}

void SNA::compute_ui(int const jnum)
{
  zero_uarraytot();

  addself_uarraytot(wself);

  // Loop over jnum neighbors
  for (int j = 0; j < jnum; ++j)
  {
    double const dx = rij(j, 0);
    double const dy = rij(j, 1);
    double const dz = rij(j, 2);
    double const rsq = dx * dx + dy * dy + dz * dz;
    double const r = std::sqrt(rsq);

    double const theta0 = (r - rmin0) * rfac0 * MY_PI / (rcutij[j] - rmin0);
    double const z0 = r / std::tan(theta0);

    compute_uarray(dx, dy, dz, z0, r, j);

    add_uarraytot(r, wj[j], rcutij[j], j);
  }
}

void SNA::compute_zi()
{
  for (int jjz = 0; jjz < idxz_max; ++jjz)
  {
    int const j1 = idxz[jjz].j1;
    int const j2 = idxz[jjz].j2;
    int const j = idxz[jjz].j;

    int const ma_min = idxz[jjz].ma_min;
    int const ma_max = idxz[jjz].ma_max;
    int const na = idxz[jjz].na;

    int const mb_min = idxz[jjz].mb_min;
    int const mb_max = idxz[jjz].mb_max;
    int const nb = idxz[jjz].nb;

    double const *cgblock = cglist.data() + idxcg_block(j1, j2, j);

    zlist_r[jjz] = 0.0;
    zlist_i[jjz] = 0.0;

    int jju1 = idxu_block[j1] + (j1 + 1) * mb_min;
    int jju2 = idxu_block[j2] + (j2 + 1) * mb_max;
    int icgb = mb_min * (j2 + 1) + mb_max;

    for (int ib = 0; ib < nb; ++ib)
    {
      double suma1_r = 0.0;
      double suma1_i = 0.0;

      double const *u1_r = &ulisttot_r[jju1];
      double const *u1_i = &ulisttot_i[jju1];
      double const *u2_r = &ulisttot_r[jju2];
      double const *u2_i = &ulisttot_i[jju2];

      int icga = ma_min * (j2 + 1) + ma_max;

      for (int ia = 0, ma1 = ma_min, ma2 = ma_max; ia < na; ++ia, ++ma1, --ma2)
      {
        suma1_r += cgblock[icga] * (u1_r[ma1] * u2_r[ma2] - u1_i[ma1] * u2_i[ma2]);
        suma1_i += cgblock[icga] * (u1_r[ma1] * u2_i[ma2] + u1_i[ma1] * u2_r[ma2]);
        icga += j2;
      }

      zlist_r[jjz] += cgblock[icgb] * suma1_r;
      zlist_i[jjz] += cgblock[icgb] * suma1_i;

      jju1 += j1 + 1;
      jju2 -= j2 + 1;
      icgb += j2;
    }
  }
}

void SNA::compute_yi(double const *const beta)
{
  // Zero yarray
  for (int j = 0; j <= twojmax; ++j)
  {
    int jju = idxu_block[j];

    for (int mb = 0; 2 * mb <= j; ++mb)
    {
      for (int ma = 0; ma <= j; ++ma, ++jju)
      {
        ylist_r[jju] = 0.0;
        ylist_i[jju] = 0.0;
      }
    }
  }

  double betaj;
  for (int jjz = 0; jjz < idxz_max; ++jjz)
  {
    int const j1 = idxz[jjz].j1;
    int const j2 = idxz[jjz].j2;
    int const j = idxz[jjz].j;

    int const ma_min = idxz[jjz].ma_min;
    int const ma_max = idxz[jjz].ma_max;
    int const na = idxz[jjz].na;

    int const mb_min = idxz[jjz].mb_min;
    int const mb_max = idxz[jjz].mb_max;
    int const nb = idxz[jjz].nb;

    double const *cgblock = cglist.data() + idxcg_block(j1, j2, j);

    // int mb = (2 * (mb_min + mb_max) - j1 - j2 + j) / 2;
    // int ma = (2 * (ma_min + ma_max) - j1 - j2 + j) / 2;

    double ztmp_r = 0.0;
    double ztmp_i = 0.0;

    int jju1 = idxu_block[j1] + (j1 + 1) * mb_min;
    int jju2 = idxu_block[j2] + (j2 + 1) * mb_max;
    int icgb = mb_min * (j2 + 1) + mb_max;

    for (int ib = 0; ib < nb; ++ib)
    {
      double suma1_r = 0.0;
      double suma1_i = 0.0;

      double const *u1_r = &ulisttot_r[jju1];
      double const *u1_i = &ulisttot_i[jju1];
      double const *u2_r = &ulisttot_r[jju2];
      double const *u2_i = &ulisttot_i[jju2];

      int icga = ma_min * (j2 + 1) + ma_max;

      for (int ia = 0, ma1 = ma_min, ma2 = ma_max; ia < na; ++ia, ++ma1, --ma2)
      {
        suma1_r += cgblock[icga] * (u1_r[ma1] * u2_r[ma2] - u1_i[ma1] * u2_i[ma2]);
        suma1_i += cgblock[icga] * (u1_r[ma1] * u2_i[ma2] + u1_i[ma1] * u2_r[ma2]);
        icga += j2;
      }

      ztmp_r += cgblock[icgb] * suma1_r;
      ztmp_i += cgblock[icgb] * suma1_i;

      jju1 += j1 + 1;
      jju2 -= j2 + 1;
      icgb += j2;
    }

    // apply to z(j1,j2,j,ma,mb) to unique element of y(j) find the right y_list[jju] and beta[jjb] entries
    // multiply and divide by j+1 factors account for multiplicity of 1, 2, or 3

    int const jju = idxz[jjz].jju;

    // pick out right beta value

    if (j >= j1)
    {
      int const jjb = idxb_block(j1, j2, j);
      if (j1 == j)
      {
        betaj = ((j2 == j) ? 3 : 2) * beta[jjb];
      }
      else
      {
        betaj = beta[jjb];
      }
    }
    else if (j >= j2)
    {
      int const jjb = idxb_block(j, j2, j1);
      if (j2 == j)
      {
        betaj = 2 * beta[jjb] * (j1 + 1) / (j + 1.0);
      }
      else
      {
        betaj = beta[jjb] * (j1 + 1) / (j + 1.0);
      }
    }
    else
    {
      int const jjb = idxb_block(j2, j, j1);
      betaj = beta[jjb] * (j1 + 1) / (j + 1.0);
    }

    ylist_r[jju] += betaj * ztmp_r;
    ylist_i[jju] += betaj * ztmp_i;
  }
}

void SNA::compute_deidrj(double *const dedr)
{
  dedr[0] = 0.0;
  dedr[1] = 0.0;
  dedr[2] = 0.0;

  for (int j = 0; j <= twojmax; ++j)
  {
    int jju = idxu_block[j];

    for (int mb = 0; 2 * mb < j; ++mb)
    {
      for (int ma = 0; ma <= j; ++ma, ++jju)
      {
        double const jjjmambyarray_r = ylist_r[jju];
        double const jjjmambyarray_i = ylist_i[jju];

        auto dudr_r = dulist_r.data_1D(jju);
        auto dudr_i = dulist_i.data_1D(jju);

        dedr[0] += dudr_r[0] * jjjmambyarray_r + dudr_i[0] * jjjmambyarray_i;
        dedr[1] += dudr_r[1] * jjjmambyarray_r + dudr_i[1] * jjjmambyarray_i;
        dedr[2] += dudr_r[2] * jjjmambyarray_r + dudr_i[2] * jjjmambyarray_i;
      }
    }

    // For j even, handle middle column

    if (j % 2 == 0)
    {
      int const mb = j / 2;
      for (int ma = 0; ma < mb; ++ma, ++jju)
      {
        double const jjjmambyarray_r = ylist_r[jju];
        double const jjjmambyarray_i = ylist_i[jju];

        auto dudr_r = dulist_r.data_1D(jju);
        auto dudr_i = dulist_i.data_1D(jju);

        dedr[0] += dudr_r[0] * jjjmambyarray_r + dudr_i[0] * jjjmambyarray_i;
        dedr[1] += dudr_r[1] * jjjmambyarray_r + dudr_i[1] * jjjmambyarray_i;
        dedr[2] += dudr_r[2] * jjjmambyarray_r + dudr_i[2] * jjjmambyarray_i;
      }

      {
        double const jjjmambyarray_r = ylist_r[jju];
        double const jjjmambyarray_i = ylist_i[jju];

        auto dudr_r = dulist_r.data_1D(jju);
        auto dudr_i = dulist_i.data_1D(jju);

        dedr[0] += (dudr_r[0] * jjjmambyarray_r + dudr_i[0] * jjjmambyarray_i) * 0.5;
        dedr[1] += (dudr_r[1] * jjjmambyarray_r + dudr_i[1] * jjjmambyarray_i) * 0.5;
        dedr[2] += (dudr_r[2] * jjjmambyarray_r + dudr_i[2] * jjjmambyarray_i) * 0.5;
      }
    }
  }

  dedr[0] *= 2.0;
  dedr[1] *= 2.0;
  dedr[2] *= 2.0;
}

void SNA::compute_bi()
{
  for (int jjb = 0; jjb < idxb_max; ++jjb)
  {
    int const j1 = idxb[jjb].j1;
    int const j2 = idxb[jjb].j2;
    int const j = idxb[jjb].j;

    int jjz = idxz_block(j1, j2, j);
    int jju = idxu_block[j];

    double sumzu = 0.0;
    for (int mb = 0; 2 * mb < j; ++mb)
    {
      for (int ma = 0; ma <= j; ++ma, ++jjz, ++jju)
      {
        sumzu += ulisttot_r[jju] * zlist_r[jjz] + ulisttot_i[jju] * zlist_i[jjz];
      }
    }

    // For j even, handle middle column

    if (j % 2 == 0)
    {
      int const mb = j / 2;
      for (int ma = 0; ma < mb; ++ma, ++jjz, ++jju)
      {
        sumzu += ulisttot_r[jju] * zlist_r[jjz] + ulisttot_i[jju] * zlist_i[jjz];
      }

      sumzu += 0.5 * (ulisttot_r[jju] * zlist_r[jjz] + ulisttot_i[jju] * zlist_i[jjz]);
    }

    blist[jjb] = 2.0 * sumzu;

    // apply bzero shift
    if (bzeroflag)
    {
      blist[jjb] -= bzero[j];
    }
  }
}

void SNA::compute_dbidrj()
{
  double sumzdu_r[3];

  for (int jjb = 0; jjb < idxb_max; ++jjb)
  {
    int const j1 = idxb[jjb].j1;
    int const j2 = idxb[jjb].j2;
    int const j = idxb[jjb].j;

    auto dbdr = dblist.data_1D(jjb);

    dbdr[0] = 0.0;
    dbdr[1] = 0.0;
    dbdr[2] = 0.0;

    // Sum terms Conj(dudr(j,ma,mb))*z(j1,j2,j,ma,mb)

    int jjz = idxz_block(j1, j2, j);
    int jju = idxu_block[j];

    sumzdu_r[0] = 0.0;
    sumzdu_r[1] = 0.0;
    sumzdu_r[2] = 0.0;

    for (int mb = 0; 2 * mb < j; ++mb)
    {
      for (int ma = 0; ma <= j; ++ma, ++jjz, ++jju)
      {
        auto dudr_r = dulist_r.data_1D(jju);
        auto dudr_i = dulist_i.data_1D(jju);

        sumzdu_r[0] += dudr_r[0] * zlist_r[jjz] + dudr_i[0] * zlist_i[jjz];
        sumzdu_r[1] += dudr_r[1] * zlist_r[jjz] + dudr_i[1] * zlist_i[jjz];
        sumzdu_r[2] += dudr_r[2] * zlist_r[jjz] + dudr_i[2] * zlist_i[jjz];
      }
    }

    // For j even, handle middle column
    if (j % 2 == 0)
    {
      int const mb = j / 2;
      for (int ma = 0; ma < mb; ++ma, ++jjz, ++jju)
      {
        auto dudr_r = dulist_r.data_1D(jju);
        auto dudr_i = dulist_i.data_1D(jju);

        sumzdu_r[0] += dudr_r[0] * zlist_r[jjz] + dudr_i[0] * zlist_i[jjz];
        sumzdu_r[1] += dudr_r[1] * zlist_r[jjz] + dudr_i[1] * zlist_i[jjz];
        sumzdu_r[2] += dudr_r[2] * zlist_r[jjz] + dudr_i[2] * zlist_i[jjz];
      }

      {
        auto dudr_r = dulist_r.data_1D(jju);
        auto dudr_i = dulist_i.data_1D(jju);

        sumzdu_r[0] += (dudr_r[0] * zlist_r[jjz] + dudr_i[0] * zlist_i[jjz]) * 0.5;
        sumzdu_r[1] += (dudr_r[1] * zlist_r[jjz] + dudr_i[1] * zlist_i[jjz]) * 0.5;
        sumzdu_r[2] += (dudr_r[2] * zlist_r[jjz] + dudr_i[2] * zlist_i[jjz]) * 0.5;
      }
    } // end if jeven

    dbdr[0] += 2.0 * sumzdu_r[0];
    dbdr[1] += 2.0 * sumzdu_r[1];
    dbdr[2] += 2.0 * sumzdu_r[2];

    // Sum over Conj(dudr(j1,ma1,mb1))*z(j,j2,j1,ma1,mb1)

    jjz = idxz_block(j, j2, j1);
    jju = idxu_block[j1];

    sumzdu_r[0] = 0.0;
    sumzdu_r[1] = 0.0;
    sumzdu_r[2] = 0.0;

    for (int mb = 0; 2 * mb < j1; ++mb)
    {
      for (int ma = 0; ma <= j1; ++ma, ++jjz, ++jju)
      {
        auto dudr_r = dulist_r.data_1D(jju);
        auto dudr_i = dulist_i.data_1D(jju);

        sumzdu_r[0] += dudr_r[0] * zlist_r[jjz] + dudr_i[0] * zlist_i[jjz];
        sumzdu_r[1] += dudr_r[1] * zlist_r[jjz] + dudr_i[1] * zlist_i[jjz];
        sumzdu_r[2] += dudr_r[2] * zlist_r[jjz] + dudr_i[2] * zlist_i[jjz];
      }
    }

    // For j1 even, handle middle column
    if (j1 % 2 == 0)
    {
      int mb = j1 / 2;
      for (int ma = 0; ma < mb; ++ma, ++jjz, ++jju)
      {
        auto dudr_r = dulist_r.data_1D(jju);
        auto dudr_i = dulist_i.data_1D(jju);

        sumzdu_r[0] += dudr_r[0] * zlist_r[jjz] + dudr_i[0] * zlist_i[jjz];
        sumzdu_r[1] += dudr_r[1] * zlist_r[jjz] + dudr_i[1] * zlist_i[jjz];
        sumzdu_r[2] += dudr_r[2] * zlist_r[jjz] + dudr_i[2] * zlist_i[jjz];
      }

      {
        auto dudr_r = dulist_r.data_1D(jju);
        auto dudr_i = dulist_i.data_1D(jju);

        sumzdu_r[0] += (dudr_r[0] * zlist_r[jjz] + dudr_i[0] * zlist_i[jjz]) * 0.5;
        sumzdu_r[1] += (dudr_r[1] * zlist_r[jjz] + dudr_i[1] * zlist_i[jjz]) * 0.5;
        sumzdu_r[2] += (dudr_r[2] * zlist_r[jjz] + dudr_i[2] * zlist_i[jjz]) * 0.5;
      }
    } // end if j1even

    double const j1fac = (j + 1) / (j1 + 1.0);

    dbdr[0] += 2.0 * sumzdu_r[0] * j1fac;
    dbdr[1] += 2.0 * sumzdu_r[1] * j1fac;
    dbdr[2] += 2.0 * sumzdu_r[2] * j1fac;

    //

    jjz = idxz_block(j, j1, j2);
    jju = idxu_block[j2];

    sumzdu_r[0] = 0.0;
    sumzdu_r[1] = 0.0;
    sumzdu_r[2] = 0.0;

    for (int mb = 0; 2 * mb < j2; ++mb)
    {
      for (int ma = 0; ma <= j2; ++ma, ++jjz, ++jju)
      {
        auto dudr_r = dulist_r.data_1D(jju);
        auto dudr_i = dulist_i.data_1D(jju);

        sumzdu_r[0] += dudr_r[0] * zlist_r[jjz] + dudr_i[0] * zlist_i[jjz];
        sumzdu_r[1] += dudr_r[1] * zlist_r[jjz] + dudr_i[1] * zlist_i[jjz];
        sumzdu_r[2] += dudr_r[2] * zlist_r[jjz] + dudr_i[2] * zlist_i[jjz];
      }
    }

    // For j2 even, handle middle column
    if (j2 % 2 == 0)
    {
      int const mb = j2 / 2;
      for (int ma = 0; ma < mb; ++ma, ++jjz, ++jju)
      {
        auto dudr_r = dulist_r.data_1D(jju);
        auto dudr_i = dulist_i.data_1D(jju);

        sumzdu_r[0] += dudr_r[0] * zlist_r[jjz] + dudr_i[0] * zlist_i[jjz];
        sumzdu_r[1] += dudr_r[1] * zlist_r[jjz] + dudr_i[1] * zlist_i[jjz];
        sumzdu_r[2] += dudr_r[2] * zlist_r[jjz] + dudr_i[2] * zlist_i[jjz];
      }

      {
        auto dudr_r = dulist_r.data_1D(jju);
        auto dudr_i = dulist_i.data_1D(jju);

        sumzdu_r[0] += (dudr_r[0] * zlist_r[jjz] + dudr_i[0] * zlist_i[jjz]) * 0.5;
        sumzdu_r[1] += (dudr_r[1] * zlist_r[jjz] + dudr_i[1] * zlist_i[jjz]) * 0.5;
        sumzdu_r[2] += (dudr_r[2] * zlist_r[jjz] + dudr_i[2] * zlist_i[jjz]) * 0.5;
      }
    } // end if j2even

    double const j2fac = (j + 1) / (j2 + 1.0);

    dbdr[0] += 2.0 * sumzdu_r[0] * j2fac;
    dbdr[1] += 2.0 * sumzdu_r[1] * j2fac;
    dbdr[2] += 2.0 * sumzdu_r[2] * j2fac;
  } //end loop over j1 j2 j
}

void SNA::compute_duidrj(double const *const rij_in, double const wj_in, double const rcut, int const atom_j)
{
  double const dx = rij_in[0];
  double const dy = rij_in[1];
  double const dz = rij_in[2];

  double const rsq = dx * dx + dy * dy + dz * dz;
  double const r = std::sqrt(rsq);

  double const rscale0 = rfac0 * MY_PI / (rcut - rmin0);
  double const theta0 = (r - rmin0) * rscale0;

  double const cs = std::cos(theta0);
  double const sn = std::sin(theta0);

  double const z0 = r * cs / sn;

  double const dz0dr = z0 / r - (r * rscale0) * (rsq + z0 * z0) / rsq;

  compute_duarray(dx, dy, dz, z0, r, dz0dr, wj_in, rcut, atom_j);
}

void SNA::zero_uarraytot()
{
  for (int j = 0; j <= twojmax; ++j)
  {
    int jju = idxu_block[j];
    for (int ma_mb = 0; ma_mb <= j * j; ++ma_mb, ++jju)
    {
      ulisttot_r[jju] = 0.0;
      ulisttot_i[jju] = 0.0;
    }
  }
}

void SNA::addself_uarraytot(double const wself_in)
{
  for (int j = 0; j <= twojmax; ++j)
  {
    int jju = idxu_block[j];
    for (int ma = 0; ma <= j; ++ma)
    {
      ulisttot_r[jju] = wself_in;
      ulisttot_i[jju] = 0.0;
      jju += j + 2;
    }
  }
}

void SNA::add_uarraytot(double const r, double const wj_in, double const rcut, int const atom_j)
{
  double sfac = compute_sfac(r, rcut);
  sfac *= wj_in;

  auto ulist_r = ulist_r_ij.data_1D(atom_j);
  auto ulist_i = ulist_i_ij.data_1D(atom_j);

  for (int j = 0; j <= twojmax; ++j)
  {
    int jju = idxu_block[j];
    for (int mb = 0; mb <= j; ++mb)
    {
      for (int ma = 0; ma <= j; ++ma, ++jju)
      {
        ulisttot_r[jju] += sfac * ulist_r[jju];
        ulisttot_i[jju] += sfac * ulist_i[jju];
      }
    }
  }
}

void SNA::compute_uarray(double const dx, double const dy, double const dz, double const z0, double const r, int const atom_j)
{
  // Compute Cayley-Klein parameters for unit quaternion

  double const r0inv = 1.0 / std::sqrt(r * r + z0 * z0);

  double const a_r = r0inv * z0;
  double const a_i = -r0inv * dz;
  double const b_r = r0inv * dy;
  double const b_i = -r0inv * dx;

  // VMK Section 4.8.2

  auto ulist_r = ulist_r_ij.data_1D(atom_j);
  auto ulist_i = ulist_i_ij.data_1D(atom_j);

  ulist_r[0] = 1.0;
  ulist_i[0] = 0.0;

  double rootpq;

  for (int j = 1; j <= twojmax; ++j)
  {
    int jju = idxu_block[j];
    int jjup = idxu_block[j - 1];

    // Fill the left side of matrix layer from previous layer

    for (int mb = 0; 2 * mb <= j; ++mb, ++jju)
    {
      ulist_r[jju] = 0.0;
      ulist_i[jju] = 0.0;

      for (int ma = 0; ma < j; ++ma, ++jju, ++jjup)
      {
        rootpq = rootpqarray(j - ma, j - mb);

        ulist_r[jju] += rootpq * (a_r * ulist_r[jjup] + a_i * ulist_i[jjup]);
        ulist_i[jju] += rootpq * (a_r * ulist_i[jjup] - a_i * ulist_r[jjup]);

        rootpq = rootpqarray(ma + 1, j - mb);

        ulist_r[jju + 1] = -rootpq * (b_r * ulist_r[jjup] + b_i * ulist_i[jjup]);
        ulist_i[jju + 1] = -rootpq * (b_r * ulist_i[jjup] - b_i * ulist_r[jjup]);
      }
    }

    // Copy the left side to the right side with inversion symmetry VMK 4.4(2)
    // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])

    jju = idxu_block[j];
    jjup = jju + (j + 1) * (j + 1) - 1;

    for (int mb = 0, mbpar = 1; 2 * mb <= j; ++mb)
    {
      for (int ma = 0, mapar = mbpar; ma <= j; ++ma, ++jju, --jjup)
      {
        if (mapar)
        {
          ulist_r[jjup] = ulist_r[jju];
          ulist_i[jjup] = -ulist_i[jju];
        }
        else
        {
          ulist_r[jjup] = -ulist_r[jju];
          ulist_i[jjup] = ulist_i[jju];
        }
        mapar = !mapar;
      }
      mbpar = !mbpar;
    }
  }
}

void SNA::compute_duarray(double const dx, double const dy, double const dz,
                          double const z0, double const r, double const dz0dr,
                          double const wj, double const rcut, int const atom_j)
{
  double da_r[3];
  double da_i[3];
  double db_r[3];
  double db_i[3];
  double dz0[3];
  double dr0inv[3];

  double rootpq;

  double const rinv = 1.0 / r;

  double const ux = dx * rinv;
  double const uy = dy * rinv;
  double const uz = dz * rinv;

  double const r0inv = 1.0 / std::sqrt(r * r + z0 * z0);

  double const a_r = z0 * r0inv;
  double const a_i = -dz * r0inv;
  double const b_r = dy * r0inv;
  double const b_i = -dx * r0inv;

  double const dr0invdr = -r0inv * r0inv * r0inv * (r + z0 * dz0dr);

  dr0inv[0] = dr0invdr * ux;
  dr0inv[1] = dr0invdr * uy;
  dr0inv[2] = dr0invdr * uz;

  dz0[0] = dz0dr * ux;
  dz0[1] = dz0dr * uy;
  dz0[2] = dz0dr * uz;

  da_r[0] = dz0[0] * r0inv + z0 * dr0inv[0];
  da_r[1] = dz0[1] * r0inv + z0 * dr0inv[1];
  da_r[2] = dz0[2] * r0inv + z0 * dr0inv[2];

  da_i[0] = -dz * dr0inv[0];
  da_i[1] = -dz * dr0inv[1];
  da_i[2] = -dz * dr0inv[2] - r0inv;

  db_r[0] = dy * dr0inv[0];
  db_r[1] = dy * dr0inv[1] + r0inv;
  db_r[2] = dy * dr0inv[2];

  db_i[0] = -dx * dr0inv[0] - r0inv;
  db_i[1] = -dx * dr0inv[1];
  db_i[2] = -dx * dr0inv[2];

  auto ulist_r = ulist_r_ij.data_1D(atom_j);
  auto ulist_i = ulist_i_ij.data_1D(atom_j);

  dulist_r(0, 0) = 0.0;
  dulist_r(0, 1) = 0.0;
  dulist_r(0, 2) = 0.0;

  dulist_i(0, 0) = 0.0;
  dulist_i(0, 1) = 0.0;
  dulist_i(0, 2) = 0.0;

  for (int j = 1; j <= twojmax; ++j)
  {
    int jju = idxu_block[j];
    int jjup = idxu_block[j - 1];

    for (int mb = 0; 2 * mb <= j; ++mb, ++jju)
    {
      dulist_r(jju, 0) = 0.0;
      dulist_r(jju, 1) = 0.0;
      dulist_r(jju, 2) = 0.0;

      dulist_i(jju, 0) = 0.0;
      dulist_i(jju, 1) = 0.0;
      dulist_i(jju, 2) = 0.0;

      for (int ma = 0; ma < j; ++ma, ++jju, ++jjup)
      {
        rootpq = rootpqarray(j - ma, j - mb);
        for (int k = 0; k < 3; ++k)
        {
          dulist_r(jju, k) += rootpq *
                              (da_r[k] * ulist_r[jjup] + da_i[k] * ulist_i[jjup] +
                               a_r * dulist_r(jjup, k) + a_i * dulist_i(jjup, k));
          dulist_i(jju, k) += rootpq *
                              (da_r[k] * ulist_i[jjup] - da_i[k] * ulist_r[jjup] +
                               a_r * dulist_i(jjup, k) - a_i * dulist_r(jjup, k));
        }

        rootpq = rootpqarray(ma + 1, j - mb);
        for (int k = 0; k < 3; ++k)
        {
          dulist_r(jju + 1, k) = -rootpq *
                                 (db_r[k] * ulist_r[jjup] + db_i[k] * ulist_i[jjup] +
                                  b_r * dulist_r(jjup, k) + b_i * dulist_i(jjup, k));
          dulist_i(jju + 1, k) = -rootpq *
                                 (db_r[k] * ulist_i[jjup] - db_i[k] * ulist_r[jjup] +
                                  b_r * dulist_i(jjup, k) - b_i * dulist_r(jjup, k));
        }
      }
    }

    // copy left side to right side with inversion symmetry VMK 4.4(2)
    // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])

    jju = idxu_block[j];
    jjup = jju + (j + 1) * (j + 1) - 1;

    for (int mb = 0, mbpar = 1; 2 * mb <= j; ++mb)
    {
      for (int ma = 0, mapar = mbpar; ma <= j; ++ma, ++jju, --jjup)
      {
        if (mapar)
        {
          for (int k = 0; k < 3; ++k)
          {
            dulist_r(jjup, k) = dulist_r(jju, k);
            dulist_i(jjup, k) = -dulist_i(jju, k);
          }
        }
        else
        {
          for (int k = 0; k < 3; ++k)
          {
            dulist_r(jjup, k) = -dulist_r(jju, k);
            dulist_i(jjup, k) = dulist_i(jju, k);
          }
        }
        mapar = !mapar;
      }
      mbpar = !mbpar;
    }
  }

  double const sfac = compute_sfac(r, rcut) * wj;

  double const dsfac = compute_dsfac(r, rcut) * wj;

  for (int j = 0; j <= twojmax; ++j)
  {
    int jju = idxu_block[j];

    for (int mb = 0; 2 * mb <= j; ++mb)
    {
      for (int ma = 0; ma <= j; ++ma, ++jju)
      {
        dulist_r(jju, 0) = dsfac * ulist_r[jju] * ux + sfac * dulist_r(jju, 0);
        dulist_r(jju, 1) = dsfac * ulist_r[jju] * uy + sfac * dulist_r(jju, 1);
        dulist_r(jju, 2) = dsfac * ulist_r[jju] * uz + sfac * dulist_r(jju, 2);

        dulist_i(jju, 0) = dsfac * ulist_i[jju] * ux + sfac * dulist_i(jju, 0);
        dulist_i(jju, 1) = dsfac * ulist_i[jju] * uy + sfac * dulist_i(jju, 1);
        dulist_i(jju, 2) = dsfac * ulist_i[jju] * uz + sfac * dulist_i(jju, 2);
      }
    }
  }
}

void SNA::create_arrays()
{
  rootpqarray.resize(twojmax + 2, twojmax + 2);

  cglist.resize(idxcg_max);

  ulisttot_r.resize(idxu_max);
  ulisttot_i.resize(idxu_max);

  dulist_r.resize(idxu_max, 3);
  dulist_i.resize(idxu_max, 3);

  zlist_r.resize(idxz_max);
  zlist_i.resize(idxz_max);

  blist.resize(idxb_max);

  dblist.resize(idxb_max, 3);

  ylist_r.resize(idxu_max);
  ylist_i.resize(idxu_max);

  if (bzeroflag)
  {
    bzero.resize(twojmax + 1);
  }
}

inline double SNA::factorial(int const n)
{
  return std::tgamma(n + 1);
}

double SNA::deltacg(int const j1, int const j2, int const j)
{
  double const factorial_1 = factorial((j1 + j2 - j) / 2);
  double const factorial_2 = factorial((j1 - j2 + j) / 2);
  double const factorial_3 = factorial((-j1 + j2 + j) / 2);
  double const sfaccg = factorial((j1 + j2 + j) / 2 + 1);

  return std::sqrt(factorial_1 * factorial_2 * factorial_3 / sfaccg);
}

void SNA::init_clebsch_gordan()
{
  for (int j1 = 0, counter = 0; j1 <= twojmax; ++j1)
  {
    for (int j2 = 0; j2 <= j1; ++j2)
    {
      for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2)
      {
        for (int m1 = 0; m1 <= j1; ++m1)
        {
          int const aa2 = 2 * m1 - j1;

          for (int m2 = 0; m2 <= j2; ++m2, ++counter)
          {
            // -c <= cc <= c

            int const bb2 = 2 * m2 - j2;
            int const m = (aa2 + bb2 + j) / 2;

            if (m < 0 || m > j)
            {
              cglist[counter] = 0.0;
              continue;
            }

            double sum = 0.0;

            for (int z = std::max(0, std::max(-(j - j2 + aa2) / 2, -(j - j1 - bb2) / 2));
                 z <= std::min((j1 + j2 - j) / 2, std::min((j1 - aa2) / 2, (j2 + bb2) / 2));
                 ++z)
            {
              double const factorial_1 = factorial(z);
              double const factorial_2 = factorial((j1 + j2 - j) / 2 - z);
              double const factorial_3 = factorial((j1 - aa2) / 2 - z);
              double const factorial_4 = factorial((j2 + bb2) / 2 - z);
              double const factorial_5 = factorial((j - j2 + aa2) / 2 + z);
              double const factorial_6 = factorial((j - j1 - bb2) / 2 + z);

              double const ifac = (z % 2) ? -1.0 : 1.0;

              sum += ifac / (factorial_1 * factorial_2 * factorial_3 * factorial_4 * factorial_5 * factorial_6);
            }

            double const dcg = deltacg(j1, j2, j);

            double sfaccg;
            {
              double const factorial_1 = factorial((j1 + aa2) / 2);
              double const factorial_2 = factorial((j1 - aa2) / 2);
              double const factorial_3 = factorial((j2 + bb2) / 2);
              double const factorial_4 = factorial((j2 - bb2) / 2);

              int const cc2 = 2 * m - j;

              double const factorial_5 = factorial((j + cc2) / 2);
              double const factorial_6 = factorial((j - cc2) / 2);

              sfaccg = factorial_1 * factorial_2 * factorial_3 * factorial_4 * factorial_5 * factorial_6 * (j + 1);
            }

            cglist[counter] = sum * dcg * std::sqrt(sfaccg);
          }
        }
      }
    }
  }
}

void SNA::init_rootpqarray()
{
  for (int p = 1; p <= twojmax; ++p)
  {
    for (int q = 1; q <= twojmax; ++q)
    {
      rootpqarray(p, q) = std::sqrt(static_cast<double>(p) / q);
    }
  }
}

int SNA::compute_ncoeff()
{
  int counter = 0;
  for (int j1 = 0; j1 <= twojmax; ++j1)
  {
    for (int j2 = 0; j2 <= j1; ++j2)
    {
      for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2)
      {
        if (j >= j1)
        {
          ++counter;
        }
      }
    }
  }
  return counter;
}

double SNA::compute_sfac(double const r, double const rcut)
{
  if (!switchflag)
  {
    return 1.0;
  }

  if (r <= rmin0)
  {
    return 1.0;
  }
  else if (r > rcut)
  {
    return 0.0;
  }
  else
  {
    double const rcutfac = MY_PI / (rcut - rmin0);
    return 0.5 * (std::cos((r - rmin0) * rcutfac) + 1.0);
  }
}

double SNA::compute_dsfac(double const r, double const rcut)
{
  if (!switchflag)
  {
    return 0.0;
  }

  if (r <= rmin0 || r > rcut)
  {
    return 0.0;
  }
  else
  {
    double const rcutfac = MY_PI / (rcut - rmin0);
    return -0.5 * std::sin((r - rmin0) * rcutfac) * rcutfac;
  }
}

#undef MY_PI
