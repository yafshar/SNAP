//
// CDDL HEADER START
//
// The contents of this file are subject to the terms of the Common Development
// and Distribution License Version 1.0 (the "License").
//
// You can obtain a copy of the license at
// http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
// specific language governing permissions and limitations under the License.
//
// When distributing Covered Code, include this CDDL HEADER in each file and
// include the License file in a prominent location with the name LICENSE.CDDL.
// If applicable, add the following below this CDDL HEADER, with the fields
// enclosed by brackets "[]" replaced with your own identifying information:
//
// Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
//
// CDDL HEADER END
//

//
// Copyright (c) 2019, Regents of the University of Minnesota.
// All rights reserved.
//
// Contributors:
//    Yaser Afshar
//    Ryan S. Elliott
//

#ifndef SNAP_IMPLEMENTATION_HPP
#define SNAP_IMPLEMENTATION_HPP

#include "KIM_ModelDriverHeaders.hpp"

#include "helper.hpp"
#include "SNA.hpp"

#include <vector>
#include <memory>

#ifdef NUM_PARAMETER_FILES
#undef NUM_PARAMETER_FILES
#endif

/*!
 * \brief SNAP has two input files, which are:
 *  1) a SNAP coefficient file, usually ends in the '.snapcoeff' extension.
 *  2) a SNAP parameter file, usually ends in the '.snapparam' extension.
 */
#define NUM_PARAMETER_FILES 2

/*! \class SNAPImplementation
 * \brief SNAP model driver Implementation class
 */
class SNAPImplementation
{
public:
  /*!
   * \brief Construct a new SNAPImplementation object
   *
   * \param modelDriverCreate A %KIM API Model object
   * \param requestedLengthUnit User LengthUnit
   * \param requestedEnergyUnit User EnergyUnit
   * \param requestedChargeUnit User ChargeUnit
   * \param requestedTemperatureUnit User TemperatureUnit
   * \param requestedTimeUnit User TimeUnit
   * \param ier 0|false if everything goes well and 1|true if it fails
   */
  SNAPImplementation(KIM::ModelDriverCreate *const modelDriverCreate,
                     KIM::LengthUnit const requestedLengthUnit,
                     KIM::EnergyUnit const requestedEnergyUnit,
                     KIM::ChargeUnit const requestedChargeUnit,
                     KIM::TemperatureUnit const requestedTemperatureUnit,
                     KIM::TimeUnit const requestedTimeUnit,
                     int *const ier);

  /*!
   * \brief Destroy the SNAPImplementation object
   *
   */
  ~SNAPImplementation();

  // no explicit Destroy() needed here

  /*!
   * \brief A refresh routine
   *
   * \param modelRefresh A refresh routine
   *
   * \return int 0|false if everything goes well and 1|true if it fails
   */
  int Refresh(KIM::ModelRefresh *const modelRefresh);

  /*!
   * \brief A routine to write the parameterized model files
   *
   * \param modelWriteParameterizedModel A %KIM API Model object
   *
   * \return int 0|false if everything goes well and 1|true if it fails
   */
  int WriteParameterizedModel(KIM::ModelWriteParameterizedModel const *const modelWriteParameterizedModel) const;

  /*!
   * \brief Compute routine
   *
   * \param modelCompute A %KIM API Model object
   * \param modelComputeArguments A %KIM API ComputeArguments object
   *
   * \return int 0|false if everything goes well and 1|true if it fails
   */
  int Compute(KIM::ModelCompute const *const modelCompute, KIM::ModelComputeArguments const *const modelComputeArguments);

  /*!
   * \brief Create a compute arguments routine
   *
   * \param modelComputeArgumentsCreate
   *
   * \return int 0|false if everything goes well and 1|true if it fails
   */
  int ComputeArgumentsCreate(KIM::ModelComputeArgumentsCreate *const modelComputeArgumentsCreate) const;

  /*!
   * \brief Destroy a compute arguments routine
   *
   * \param modelComputeArgumentsDestroy
   *
   * \return int 0|false if everything goes well and 1|true if it fails
   */
  int ComputeArgumentsDestroy(KIM::ModelComputeArgumentsDestroy *const modelComputeArgumentsDestroy) const;

private:
  /*!
   * \brief Open the SNAP parameter files
   *
   * \param modelDriverCreate A %KIM API Model object
   * \param numberParameterFiles Number of parameter files to open
   * \param parameterFilePointers FILE pointer to the opened files
   *
   * \return int 0|false if everything goes well and 1|true if it fails
   */
  int OpenParameterFiles(KIM::ModelDriverCreate *const modelDriverCreate,
                         int const numberParameterFiles,
                         std::FILE **parameterFilePointers);

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
   * \brief Parse the input files and process the read parameters
   *
   * \param modelDriverCreate A %KIM API Model object
   * \param numberParameterFiles Number of parameter files to open
   * \param parameterFilePointers FILE pointer to the opened files
   *
   * \return int 0|false if everything goes well and 1|true if it fails
   */
  int ProcessParameterFiles(KIM::ModelDriverCreate *const modelDriverCreate,
                            int const numberParameterFiles,
                            std::FILE *const *parameterFilePointers);

  /*!
   * \brief Close all the opened files
   *
   * \param numberParameterFiles Number of parameter files to close
   * \param parameterFilePointers FILE pointer to the opened files
   */
  void CloseParameterFiles(int const numberParameterFiles,
                           std::FILE *const *parameterFilePointers);

  /*!
   * \brief Convert units of the parameters from the input files
   *
   * \param modelDriverCreate A %KIM API Model object
   * \param requestedLengthUnit
   * \param requestedEnergyUnit
   * \param requestedChargeUnit
   * \param requestedTemperatureUnit
   * \param requestedTimeUnit
   *
   * \return int 0|false if everything goes well and 1|true if it fails
   */
  int ConvertUnits(KIM::ModelDriverCreate *const modelDriverCreate,
                   KIM::LengthUnit const &requestedLengthUnit,
                   KIM::EnergyUnit const &requestedEnergyUnit,
                   KIM::ChargeUnit const &requestedChargeUnit,
                   KIM::TemperatureUnit const &requestedTemperatureUnit,
                   KIM::TimeUnit const &requestedTimeUnit);

  /*!
   * \brief Set the Compute Mutable Values object
   *
   * \param modelComputeArguments A %KIM API ComputeArguments object
   * \param isComputeProcess_dEdr ComputeProcess_dEdr flag
   * \param isComputeProcess_d2Edr2 ComputeProcess_d2Edr2 flag
   * \param isComputeEnergy ComputeEnergy flag
   * \param isComputeForces ComputeForces flag
   * \param isComputeParticleEnergy ComputeParticleEnergy flag
   * \param isComputeVirial ComputeVirial flag
   * \param isComputeParticleVirial ComputeParticleVirial flag
   * \param particleSpeciesCodes Particle species code
   * \param particleContributing Particle contirubuting flag list
   * \param coordinates Particles' coordinates
   * \param energy System energy
   * \param forces Particles' forces
   * \param particleEnergy Particles' energy
   * \param virial System virial
   * \param particleVirial Particles' virial
   *
   * \return int 0|false if everything goes well and 1|true if it fails
   */
  int setComputeMutableValues(KIM::ModelComputeArguments const *const modelComputeArguments,
                              bool &isComputeProcess_dEdr,
                              bool &isComputeProcess_d2Edr2,
                              bool &isComputeEnergy,
                              bool &isComputeForces,
                              bool &isComputeParticleEnergy,
                              bool &isComputeVirial,
                              bool &isComputeParticleVirial,
                              int const *&particleSpeciesCodes,
                              int const *&particleContributing,
                              VectorOfSizeDIM const *&coordinates,
                              double *&energy,
                              VectorOfSizeDIM *&forces,
                              double *&particleEnergy,
                              VectorOfSizeSix *&virial,
                              VectorOfSizeSix *&particleVirial);

  /*!
   * \brief Get the Compute Index
   *
   * \param isComputeProcess_dEdr ComputeProcess_dEdr flag
   * \param isComputeProcess_d2Edr2 ComputeProcess_d2Edr2 flag
   * \param isComputeEnergy ComputeEnergy flag
   * \param isComputeForces ComputeForces flag
   * \param isComputeParticleEnergy ComputeParticleEnergy flag
   * \param isComputeVirial ComputeVirial flag
   * \param isComputeParticleVirial ComputeParticleVirial flag
   *
   * \return int The computed Index
   */
  int getComputeIndex(bool const isComputeProcess_dEdr,
                      bool const isComputeProcess_d2Edr2,
                      bool const isComputeEnergy,
                      bool const isComputeForces,
                      bool const isComputeParticleEnergy,
                      bool const isComputeVirial,
                      bool const isComputeParticleVirial) const;

  /*!
   * \brief Use (possibly) new values of parameters to compute other quantities
   *
   * \tparam ModelObj A %KIM API Model object
   *
   * \param modelObj It is a %KIM API ModelDriverCreate object during initialization or
   *                 a ModelRefresh object when the Model's parameters have been altered
   *
   * \return int 0|false if everything goes well and 1|true if it fails
   */
  template <class ModelObj>
  int setRefreshMutableValues(ModelObj *const modelObj);

  /*!
   * \brief Set the Model's particle Numbering
   *
   * \param modelDriverCreate A %KIM API Model object
   *
   * \return int 0|false if everything goes well and 1|true if it fails
   */
  int RegisterKIMModelSettings(KIM::ModelDriverCreate *const modelDriverCreate) const;

  /*!
   * \brief Register %KIM API parameters
   *
   * \param modelDriverCreate A %KIM API Model object
   *
   * \return int 0|false if everything goes well and 1|true if it fails
   */
  int RegisterKIMParameters(KIM::ModelDriverCreate *const modelDriverCreate);

  /*!
   * \brief Register %KIM API Functions
   *
   * \param modelDriverCreate A %KIM API Model object
   *
   * \return int 0|false if everything goes well and 1|true if it fails
   */
  int RegisterKIMFunctions(KIM::ModelDriverCreate *const modelDriverCreate) const;

  /*!
   * \brief Register %KIM API arguments settings
   *
   * \param modelComputeArgumentsCreate
   *
   * \return int 0|false if everything goes well and 1|true if it fails
   */
  int RegisterKIMComputeArgumentsSettings(KIM::ModelComputeArgumentsCreate *const modelComputeArgumentsCreate) const;

  /*!
   * \brief Compute the bispectrum component of each atom
   *
   * \param modelComputeArguments The %KIM API compute argument object
   * \param particleSpeciesCodes The array of particle species code
   * \param particleContributing The array of particle contributing
   * \param coordinates The KIM \c coordinates argument (atoms' coordinates)
   */
  void computeBispectrum(KIM::ModelComputeArguments const *const modelComputeArguments,
                         int const *particleSpeciesCodes,
                         int const *particleContributing,
                         VectorOfSizeDIM const *coordinates);

  /*!
   * \brief Compute the linear coefficient corresponding to
   * the bispectrum component of each atom
   *
   * \param particleSpeciesCodes The array of particle species code
   * \param particleContributing The array of particle contributing
   */
  void computeBeta(int const *particleSpeciesCodes, int const *particleContributing);

  /*!
   * \brief Main compute routine
   *
   * \tparam isComputeProcess_dEdr ComputeProcess_dEdr flag
   * \tparam isComputeProcess_d2Edr2ComputeProcess_d2Edr2 flag
   * \tparam isComputeEnergy ComputeEnergy flag
   * \tparam isComputeForces ComputeForces flag
   * \tparam isComputeParticleEnergy ComputeParticleEnergy flag
   * \tparam isComputeVirial ComputeVirial flag
   * \tparam isComputeParticleVirial ComputeParticleVirial flag
   *
   * \param modelCompute A %KIM API Model object
   * \param modelComputeArguments A %KIM API ComputeArguments object
   * \param particleSpeciesCodes Species code
   * \param particleContributing Contribution array
   * \param coordinates Particles coordinates
   * \param energy Global energy value
   * \param forces Particles forces
   * \param particleEnergy Particles energy
   * \param virial Global virial term
   * \param particleVirial Particles virial
   *
   * \return int 0|false if everything goes well and 1|true if it fails
   */
  template <bool isComputeProcess_dEdr,
            bool isComputeProcess_d2Edr2,
            bool isComputeEnergy,
            bool isComputeForces,
            bool isComputeParticleEnergy,
            bool isComputeVirial,
            bool isComputeParticleVirial>
  int Compute(KIM::ModelCompute const *const modelCompute,
              KIM::ModelComputeArguments const *const modelComputeArguments,
              int const *const particleSpeciesCodes,
              int const *const particleContributing,
              const VectorOfSizeDIM *const coordinates,
              double *const energy,
              VectorOfSizeDIM *const forces,
              double *const particleEnergy,
              VectorOfSizeSix virial,
              VectorOfSizeSix *const particleVirial) const;

private:
  /*!
   * \brief Number of particles
   *
   * \note
   * This is a Mutable value that can change with each call to \b Refresh() and \b Compute().
   */
  int cachedNumberOfParticles_;

  /*!
   * \brief Number of contributing particles
   *
   * \note
   * This is a Mutable value that can change with each call to \b Refresh() and \b Compute().
   */
  int numberOfContributingParticles_;

  /*!
   * \brief Flag indicating if we need to request neighbors from non-contributing particles
   *
   * \note
   * This is a Mutable value that can change with each call to \b Refresh() and \b Compute().
   */
  int modelWillNotRequestNeighborsOfNoncontributingParticles_;

  /*!
   * \brief Cutoff value in %KIM API object
   *
   * \note
   * This is a Mutable value that can change with each call to \b Refresh() and \b Compute().
   */
  double influenceDistance_;

  /*! Number of unique elements */
  int nelements;

  /*! Number of all coefficients */
  int ncoeffall;

  /*! Band limit for bispectrum components (non-negative integer) */
  int twojmax;

  /*! Number of coefficients */
  int ncoeff;

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
   * \brief quadraticflag flag
   *
   * quadraticflag value = 0 or 1
   * 0 = do not generate quadratic terms
   * 1 = generate quadratic terms
   */
  int quadraticflag;

  /*! length of beta */
  int beta_max;

  /*! Parameter in distance to angle conversion [optional keyword] */
  double rfac0;

  /*! Parameter in distance (distance units) [optional keyword] */
  double rmin0;

  /*! Scale factor applied to all cutoff radii (0 < rcutfac < 1) (positive real) [optional keyword] */
  double rcutfac;

  /*! Names of unique elements */
  std::vector<std::string> elements;

  /*! Element radii */
  std::vector<double> radelem;

  /*! Element weights */
  std::vector<double> wjelem;

  /*! Element bispectrum coefficients */
  Array2D<double> coeffelem;

  /*! betas for all atoms in list */
  Array2D<double> beta;

  /*!
   * \brief Bispectrum components for all atoms in list
   *
   * Bispectrum components of an atom are order parameters characterizing
   * the radial and angular distribution of neighbor atoms. Geometry described
   * by bispectrum components.
   */
  Array2D<double> bispectrum;

  /*! Cutoff square for pair interactions */
  Array2D<double> cutsq;

  /*! SNAP object */
  std::unique_ptr<SNA> snap;
};

#endif // SNAP_IMPLEMENTATION_HPP
