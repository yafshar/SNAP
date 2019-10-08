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
