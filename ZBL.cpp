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
   Contributing authors: Stephen Foiles, Aidan Thompson (SNL)
------------------------------------------------------------------------- */

//
// ZBL.cpp
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
//        `lammps/src/pair_zbl.cpp` and amended and updated by
//        Yaser Afshar
//


#include "ZBL.hpp"

#include <cmath>

#include <utility>

ZBL::ZBL() : cut_inner(0.0),
             cut_innersq(0.0),
             cut_global(0.0),
             cut_globalsq(0.0)
{
}

ZBL::ZBL(double const inner, double const outer) : cut_inner(inner),
                                                   cut_innersq(inner * inner),
                                                   cut_global(outer),
                                                   cut_globalsq(outer * outer)
{
}

ZBL::~ZBL() {}

ZBL::ZBL(ZBL &&other)
{
    cut_inner = other.cut_inner;
    cut_innersq = other.cut_innersq;
    cut_global = other.cut_global;
    cut_globalsq = other.cut_globalsq;
    d1a = std::move(other.d1a);
    d2a = std::move(other.d2a);
    d3a = std::move(other.d3a);
    d4a = std::move(other.d4a);
    zze = std::move(other.zze);
    sw1 = std::move(other.sw1);
    sw2 = std::move(other.sw2);
    sw3 = std::move(other.sw3);
    sw4 = std::move(other.sw4);
    sw5 = std::move(other.sw5);
}

ZBL &ZBL::operator=(ZBL &&other)
{
    cut_inner = other.cut_inner;
    cut_innersq = other.cut_innersq;
    cut_global = other.cut_global;
    cut_globalsq = other.cut_globalsq;
    d1a = std::move(other.d1a);
    d2a = std::move(other.d2a);
    d3a = std::move(other.d3a);
    d4a = std::move(other.d4a);
    zze = std::move(other.zze);
    sw1 = std::move(other.sw1);
    sw2 = std::move(other.sw2);
    sw3 = std::move(other.sw3);
    sw4 = std::move(other.sw4);
    sw5 = std::move(other.sw5);
    return *this;
}

void ZBL::allocate(int const newSize)
{
    d1a.resize(newSize, newSize, 0.0);
    d2a.resize(newSize, newSize, 0.0);
    d3a.resize(newSize, newSize, 0.0);
    d4a.resize(newSize, newSize, 0.0);
    zze.resize(newSize, newSize, 0.0);
    sw1.resize(newSize, newSize, 0.0);
    sw2.resize(newSize, newSize, 0.0);
    sw3.resize(newSize, newSize, 0.0);
    sw4.resize(newSize, newSize, 0.0);
    sw5.resize(newSize, newSize, 0.0);
}

void ZBL::set_coeff(int const iSpecies,
                    int const jSpecies,
                    double const z_iSpecies,
                    double const z_jSpecies,
                    double const angstrom,
                    double const qqr2e,
                    double const qelectron)
{
    double const ainv = (std::pow(z_iSpecies, ZBLConstants::pzbl) + std::pow(z_jSpecies, ZBLConstants::pzbl)) / (ZBLConstants::a0 * angstrom);

    d1a(iSpecies, jSpecies) = ZBLConstants::d1 * ainv;
    d2a(iSpecies, jSpecies) = ZBLConstants::d2 * ainv;
    d3a(iSpecies, jSpecies) = ZBLConstants::d3 * ainv;
    d4a(iSpecies, jSpecies) = ZBLConstants::d4 * ainv;

    zze(iSpecies, jSpecies) = z_iSpecies * z_jSpecies * qqr2e * qelectron * qelectron;

    if (iSpecies != jSpecies)
    {
        d1a(jSpecies, iSpecies) = d1a(iSpecies, jSpecies);
        d2a(jSpecies, iSpecies) = d2a(iSpecies, jSpecies);
        d3a(jSpecies, iSpecies) = d3a(iSpecies, jSpecies);
        d4a(jSpecies, iSpecies) = d4a(iSpecies, jSpecies);
        zze(jSpecies, iSpecies) = zze(iSpecies, jSpecies);
    }

    double const tc = cut_global - cut_inner;

    double const fc = e_zbl(cut_global, iSpecies, jSpecies);
    double const fcp = dzbldr(cut_global, iSpecies, jSpecies);
    double const fcpp = d2zbldr2(cut_global, iSpecies, jSpecies);

    double const swa = (-3.0 * fcp + tc * fcpp) / (tc * tc);
    double const swb = (2.0 * fcp - tc * fcpp) / (tc * tc * tc);
    double const swc = -fc + (tc / 2.0) * fcp - (tc * tc / 12.0) * fcpp;

    sw1(iSpecies, jSpecies) = swa;
    sw2(iSpecies, jSpecies) = swb;
    sw3(iSpecies, jSpecies) = swa / 3.0;
    sw4(iSpecies, jSpecies) = swb / 4.0;
    sw5(iSpecies, jSpecies) = swc;

    if (iSpecies != jSpecies)
    {
        sw1(jSpecies, iSpecies) = sw1(iSpecies, jSpecies);
        sw2(jSpecies, iSpecies) = sw2(iSpecies, jSpecies);
        sw3(jSpecies, iSpecies) = sw3(iSpecies, jSpecies);
        sw4(jSpecies, iSpecies) = sw4(iSpecies, jSpecies);
        sw5(jSpecies, iSpecies) = sw5(iSpecies, jSpecies);
    }
}

double ZBL::e_zbl(double const rij, int const iSpecies, int const jSpecies)
{
    double const d1aij = d1a(iSpecies, jSpecies);
    double const d2aij = d2a(iSpecies, jSpecies);
    double const d3aij = d3a(iSpecies, jSpecies);
    double const d4aij = d4a(iSpecies, jSpecies);
    double const zzeij = zze(iSpecies, jSpecies);

    double const sum1 = ZBLConstants::c1 * std::exp(-d1aij * rij);
    double const sum2 = ZBLConstants::c2 * std::exp(-d2aij * rij);
    double const sum3 = ZBLConstants::c3 * std::exp(-d3aij * rij);
    double const sum4 = ZBLConstants::c4 * std::exp(-d4aij * rij);
    double const sum = sum1 + sum2 + sum3 + sum4;

    double const rinv = 1.0 / rij;

    double const result = zzeij * sum * rinv;
    return result;
}

double ZBL::dzbldr(double const rij, int const iSpecies, int const jSpecies)
{
    double const d1aij = d1a(iSpecies, jSpecies);
    double const d2aij = d2a(iSpecies, jSpecies);
    double const d3aij = d3a(iSpecies, jSpecies);
    double const d4aij = d4a(iSpecies, jSpecies);

    double const zzeij = zze(iSpecies, jSpecies);

    double const e1 = std::exp(-d1aij * rij);
    double const e2 = std::exp(-d2aij * rij);
    double const e3 = std::exp(-d3aij * rij);
    double const e4 = std::exp(-d4aij * rij);

    double const sum1 = ZBLConstants::c1 * e1;
    double const sum2 = ZBLConstants::c2 * e2;
    double const sum3 = ZBLConstants::c3 * e3;
    double const sum4 = ZBLConstants::c4 * e4;
    double const sum = sum1 + sum2 + sum3 + sum4;

    double const sum_p1 = ZBLConstants::c1 * d1aij * e1;
    double const sum_p2 = ZBLConstants::c2 * d2aij * e2;
    double const sum_p3 = ZBLConstants::c3 * d3aij * e3;
    double const sum_p4 = ZBLConstants::c4 * d4aij * e4;
    double const sum_p = sum_p1 + sum_p2 + sum_p3 + sum_p4;

    double const rinv = 1.0 / rij;

    double const result = zzeij * (-sum_p - sum * rinv) * rinv;

    return result;
}

double ZBL::d2zbldr2(double const rij, int const iSpecies, int const jSpecies)
{
    double const d1aij = d1a(iSpecies, jSpecies);
    double const d2aij = d2a(iSpecies, jSpecies);
    double const d3aij = d3a(iSpecies, jSpecies);
    double const d4aij = d4a(iSpecies, jSpecies);
    double const zzeij = zze(iSpecies, jSpecies);

    double const e1 = std::exp(-d1aij * rij);
    double const e2 = std::exp(-d2aij * rij);
    double const e3 = std::exp(-d3aij * rij);
    double const e4 = std::exp(-d4aij * rij);

    double const sum1 = ZBLConstants::c1 * e1;
    double const sum2 = ZBLConstants::c2 * e2;
    double const sum3 = ZBLConstants::c3 * e3;
    double const sum4 = ZBLConstants::c4 * e4;
    double const sum = sum1 + sum2 + sum3 + sum4;

    double const sum_p1 = ZBLConstants::c1 * e1 * d1aij;
    double const sum_p2 = ZBLConstants::c2 * e2 * d2aij;
    double const sum_p3 = ZBLConstants::c3 * e3 * d3aij;
    double const sum_p4 = ZBLConstants::c4 * e4 * d4aij;
    double const sum_p = sum_p1 + sum_p2 + sum_p3 + sum_p4;

    double const sum_pp1 = ZBLConstants::c1 * e1 * d1aij * d1aij;
    double const sum_pp2 = ZBLConstants::c2 * e2 * d2aij * d2aij;
    double const sum_pp3 = ZBLConstants::c3 * e3 * d3aij * d3aij;
    double const sum_pp4 = ZBLConstants::c4 * e4 * d4aij * d4aij;
    double const sum_pp = sum_pp1 + sum_pp2 + sum_pp3 + sum_pp4;

    double const rinv = 1.0 / rij;

    double const result = zzeij * (sum_pp + 2.0 * sum_p * rinv + 2.0 * sum * rinv * rinv) * rinv;

    return result;
}
