//
// SNAP.hpp
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
//    Ryan S. Elliott
//    Yaser Afshar
//


#ifndef SNAP_HPP
#define SNAP_HPP

#include "KIM_ModelDriverHeaders.hpp"

#include <memory>

extern "C"
{
  int model_driver_create(KIM::ModelDriverCreate *const modelDriverCreate,
                          KIM::LengthUnit const requestedLengthUnit,
                          KIM::EnergyUnit const requestedEnergyUnit,
                          KIM::ChargeUnit const requestedChargeUnit,
                          KIM::TemperatureUnit const requestedTemperatureUnit,
                          KIM::TimeUnit const requestedTimeUnit);
}

// Forward declaration
class SNAPImplementation;

/*!
 * \brief SNAP model driver class for %KIM API
 *
 *
 * \note
 * There is no need to make these "extern" since KIM will only access them
 * via function pointers.  "static" is required so that there is not
 * n implicit this pointer added to the prototype by the C++ compiler
 */
class SNAP
{
public:
  /*!
   * \brief Construct a new SNAP object
   *
   * \param modelDriverCreate
   * \param requestedLengthUnit
   * \param requestedEnergyUnit
   * \param requestedChargeUnit
   * \param requestedTemperatureUnit
   * \param requestedTimeUnit
   * \param ierr
   */
  SNAP(KIM::ModelDriverCreate *const modelDriverCreate,
       KIM::LengthUnit const requestedLengthUnit,
       KIM::EnergyUnit const requestedEnergyUnit,
       KIM::ChargeUnit const requestedChargeUnit,
       KIM::TemperatureUnit const requestedTemperatureUnit,
       KIM::TimeUnit const requestedTimeUnit,
       int *const ierr);

  /*!
   * \brief Destroy the SNAP object
   *
   */
  ~SNAP();

  static int Destroy(KIM::ModelDestroy *const modelDestroy);

  static int Refresh(KIM::ModelRefresh *const modelRefresh);

  static int WriteParameterizedModel(KIM::ModelWriteParameterizedModel const *const modelWriteParameterizedModel);

  static int Compute(KIM::ModelCompute const *const modelCompute, KIM::ModelComputeArguments const *const modelComputeArguments);

  static int ComputeArgumentsCreate(KIM::ModelCompute const *const modelCompute, KIM::ModelComputeArgumentsCreate *const modelComputeArgumentsCreate);

  static int ComputeArgumentsDestroy(KIM::ModelCompute const *const modelCompute, KIM::ModelComputeArgumentsDestroy *const modelComputeArgumentsDestroy);

private:
  std::unique_ptr<SNAPImplementation> snap_implementation;
};

#endif // SNAP_HPP
