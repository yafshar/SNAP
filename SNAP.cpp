//
// SNAP.cpp
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


#include "SNAP.hpp"
#include "SNAPImplementation.hpp"

// The standard interface to KIM Model Drivers

extern "C"
{
  int model_driver_create(KIM::ModelDriverCreate *const modelDriverCreate,
                          KIM::LengthUnit const requestedLengthUnit,
                          KIM::EnergyUnit const requestedEnergyUnit,
                          KIM::ChargeUnit const requestedChargeUnit,
                          KIM::TemperatureUnit const requestedTemperatureUnit,
                          KIM::TimeUnit const requestedTimeUnit)
  {
    if (!modelDriverCreate)
    {
      HELPER_LOG_ERROR("The ModelDriverCreate pointer is not assigned");
      return true;
    }

    int ier;

    // read input files, convert units if needed, compute
    // interpolation coefficients, set cutoff, and publish parameters
    SNAP *const modelObject = new SNAP(modelDriverCreate,
                                       requestedLengthUnit,
                                       requestedEnergyUnit,
                                       requestedChargeUnit,
                                       requestedTemperatureUnit,
                                       requestedTimeUnit,
                                       &ier);
    if (ier)
    {
      // constructor has already reported the error
      delete modelObject;
      return true;
    }

    // register pointer to SNAP object in KIM object
    modelDriverCreate->SetModelBufferPointer(static_cast<void *>(modelObject));

    // everything is good
    return false;
  }
}

// Implementation of SNAP public wrapper functions

SNAP::SNAP(KIM::ModelDriverCreate *const modelDriverCreate,
           KIM::LengthUnit const requestedLengthUnit,
           KIM::EnergyUnit const requestedEnergyUnit,
           KIM::ChargeUnit const requestedChargeUnit,
           KIM::TemperatureUnit const requestedTemperatureUnit,
           KIM::TimeUnit const requestedTimeUnit,
           int *const ier) : snap_implementation(new SNAPImplementation(modelDriverCreate,
                                                                        requestedLengthUnit,
                                                                        requestedEnergyUnit,
                                                                        requestedChargeUnit,
                                                                        requestedTemperatureUnit,
                                                                        requestedTimeUnit,
                                                                        ier)) {}

SNAP::~SNAP() {}

// static member function
int SNAP::Destroy(KIM::ModelDestroy *const modelDestroy)
{
  if (!modelDestroy)
  {
    HELPER_LOG_ERROR("The ModelDestroy pointer is not assigned");
    return true;
  }

  SNAP *modelObject = NULL;

  modelDestroy->GetModelBufferPointer(reinterpret_cast<void **>(&modelObject));

  if (modelObject)
  {
    // delete the object itself
    delete modelObject;
  }

  // everything is good
  return false;
}

// static member function
int SNAP::Refresh(KIM::ModelRefresh *const modelRefresh)
{
  if (!modelRefresh)
  {
    HELPER_LOG_ERROR("The ModelRefresh pointer is not assigned");
    return true;
  }

  SNAP *modelObject = NULL;

  modelRefresh->GetModelBufferPointer(reinterpret_cast<void **>(&modelObject));

  if (modelObject)
  {
    return modelObject->snap_implementation->Refresh(modelRefresh);
  }
  else
  {
    HELPER_LOG_ERROR("The Model pointer returned from GetModelBufferPointer is not assigned");
    return true;
  }
}

int SNAP::WriteParameterizedModel(KIM::ModelWriteParameterizedModel const *const modelWriteParameterizedModel)
{
  if (!modelWriteParameterizedModel)
  {
    HELPER_LOG_ERROR("The ModelWriteParameterizedModel pointer is not assigned");
    return true;
  }

  SNAP *modelObject = NULL;

  modelWriteParameterizedModel->GetModelBufferPointer(reinterpret_cast<void **>(&modelObject));

  if (modelObject)
  {
    return modelObject->snap_implementation->WriteParameterizedModel(modelWriteParameterizedModel);
  }
  else
  {
    HELPER_LOG_ERROR("The Model pointer returned from GetModelBufferPointer is not assigned");
    return true;
  }
}

// static member function
int SNAP::Compute(KIM::ModelCompute const *const modelCompute,
                  KIM::ModelComputeArguments const *const modelComputeArguments)
{
  if (!modelCompute || !modelComputeArguments)
  {
    if (!modelCompute)
    {
      HELPER_LOG_ERROR("The ModelCompute object pointer is not assigned");
    }
    if (!modelComputeArguments)
    {
      HELPER_LOG_ERROR("The ModelComputeArguments object pointer is not assigned");
    }
    return true;
  }

  SNAP *modelObject = NULL;

  modelCompute->GetModelBufferPointer(reinterpret_cast<void **>(&modelObject));

  if (modelObject)
  {
    return modelObject->snap_implementation->Compute(modelCompute, modelComputeArguments);
  }
  else
  {
    HELPER_LOG_ERROR("The Model pointer returned from GetModelBufferPointer is not assigned");
    return true;
  }
}

// static member function
int SNAP::ComputeArgumentsCreate(KIM::ModelCompute const *const modelCompute,
                                 KIM::ModelComputeArgumentsCreate *const modelComputeArgumentsCreate)
{
  if (!modelCompute || !modelComputeArgumentsCreate)
  {
    if (!modelCompute)
    {
      HELPER_LOG_ERROR("The ModelCompute object pointer is not assigned");
    }
    if (!modelComputeArgumentsCreate)
    {
      HELPER_LOG_ERROR("The ModelComputeArgumentsCreate object pointer is not assigned");
    }
    return true;
  }

  SNAP *modelObject = NULL;

  modelCompute->GetModelBufferPointer(reinterpret_cast<void **>(&modelObject));

  if (modelObject)
  {
    return modelObject->snap_implementation->ComputeArgumentsCreate(modelComputeArgumentsCreate);
  }
  else
  {
    HELPER_LOG_ERROR("The Model pointer returned from GetModelBufferPointer is not assigned");
    return true;
  }
}

// static member function
int SNAP::ComputeArgumentsDestroy(KIM::ModelCompute const *const modelCompute,
                                  KIM::ModelComputeArgumentsDestroy *const modelComputeArgumentsDestroy)
{
  if (!modelCompute || !modelComputeArgumentsDestroy)
  {
    if (!modelCompute)
    {
      HELPER_LOG_ERROR("The ModelCompute object pointer is not assigned");
    }
    if (!modelComputeArgumentsDestroy)
    {
      HELPER_LOG_ERROR("The ModelComputeArgumentsDestroy object pointer is not assigned");
    }
    return true;
  }

  SNAP *modelObject = NULL;

  modelCompute->GetModelBufferPointer(reinterpret_cast<void **>(&modelObject));

  if (modelObject)
  {
    return modelObject->snap_implementation->ComputeArgumentsDestroy(modelComputeArgumentsDestroy);
  }
  else
  {
    HELPER_LOG_ERROR("The Model pointer returned from GetModelBufferPointer is not assigned");
    return true;
  }
}

#undef HELPER_LOG_ERROR
