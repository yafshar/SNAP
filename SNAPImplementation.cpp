//
// SNAPImplementation.cpp
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
//    Ryan S. Elliott
//


#include "KIM_ModelDriverHeaders.hpp"
#include "KIM_LogMacros.hpp"

#include "SNAP.hpp"
#include "SNAPImplementation.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <string>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <map>
#include <algorithm>
#include <numeric>
#include <limits>

#ifdef MAXLINE
#undef MAXLINE
#endif
#define MAXLINE 2048

#ifdef ONE
#undef ONE
#endif
#define ONE 1.0

#ifdef KIM_LOGGER_OBJECT_NAME
#undef KIM_LOGGER_OBJECT_NAME
#endif

// public member functions

#define KIM_LOGGER_OBJECT_NAME modelDriverCreate
SNAPImplementation::SNAPImplementation(
    KIM::ModelDriverCreate *const modelDriverCreate,
    KIM::LengthUnit const requestedLengthUnit,
    KIM::EnergyUnit const requestedEnergyUnit,
    KIM::ChargeUnit const requestedChargeUnit,
    KIM::TemperatureUnit const requestedTemperatureUnit,
    KIM::TimeUnit const requestedTimeUnit,
    int *const ier) : cachedNumberOfParticles_(0),
                      numberOfContributingParticles_(0),
                      modelWillNotRequestNeighborsOfNoncontributingParticles_(1),
                      influenceDistance_(0),
                      nelements(0),
                      ncoeffall(0),
                      twojmax(0),
                      ncoeff(0),
                      switchflag(1),    // Set defaults for optional keywords
                      bzeroflag(1),     // Set defaults for optional keywords
                      quadraticflag(0), // Set defaults for optional keywords
                      beta_max(0),
                      rfac0(0.99363), // Set defaults for optional keywords
                      rmin0(0.0),
                      rcutfac(0.0)
{
  // Everything is good
  *ier = false;

  if (!modelDriverCreate)
  {
    HELPER_LOG_ERROR("The ModelDriverCreate object pointer is not assigned \n");
    *ier = true;
    return;
  }

  {
    std::FILE *parameterFilePointers[NUM_PARAMETER_FILES];

    int numberParameterFiles(0);

    // Getting the number of parameter files, open the files and process them
    modelDriverCreate->GetNumberOfParameterFiles(&numberParameterFiles);
    if (numberParameterFiles != NUM_PARAMETER_FILES)
    {
      LOG_ERROR("Wrong number of parameter files! \n");
      return;
    }

    *ier = OpenParameterFiles(modelDriverCreate, numberParameterFiles, parameterFilePointers);
    if (*ier)
      return;

    *ier = ProcessParameterFiles(modelDriverCreate, numberParameterFiles, parameterFilePointers);

    CloseParameterFiles(numberParameterFiles, parameterFilePointers);

    if (*ier)
      return;
  }

  *ier = ConvertUnits(modelDriverCreate,
                      requestedLengthUnit,
                      requestedEnergyUnit,
                      requestedChargeUnit,
                      requestedTemperatureUnit,
                      requestedTimeUnit);
  if (*ier)
    return;

  *ier = setRefreshMutableValues(modelDriverCreate);
  if (*ier)
    return;

  *ier = RegisterKIMModelSettings(modelDriverCreate);
  if (*ier)
    return;

  *ier = RegisterKIMParameters(modelDriverCreate);
  if (*ier)
    return;

  *ier = RegisterKIMFunctions(modelDriverCreate);
  if (*ier)
    return;
}
#undef KIM_LOGGER_OBJECT_NAME

SNAPImplementation::~SNAPImplementation()
{
}

int SNAPImplementation::Refresh(KIM::ModelRefresh *const modelRefresh)
{
  int ier = setRefreshMutableValues(modelRefresh);
  return ier;
}

int SNAPImplementation::WriteParameterizedModel(KIM::ModelWriteParameterizedModel const *const modelWriteParameterizedModel) const
{
  if (!modelWriteParameterizedModel)
  {
    HELPER_LOG_ERROR("The modelWriteParameterizedModel object pointer is not assigned \n");
    return true;
  }

  std::string const *path = NULL;
  std::string const *modelName = NULL;

  modelWriteParameterizedModel->GetPath(&path);
  modelWriteParameterizedModel->GetModelName(&modelName);

  std::string coefficient_file(*modelName);
  std::string parameter_file(*modelName);

  coefficient_file += ".snapcoeff";
  parameter_file += ".snapparam";

  // Set the file name for the next parameter file.
  // It must be called once for each parameter file. The order of these calls is important
  // and determines the order in which the parameter files will be listed in the automatically
  // generated CMakeLists.txt file.
  modelWriteParameterizedModel->SetParameterFileName(coefficient_file);
  modelWriteParameterizedModel->SetParameterFileName(parameter_file);

  // SNAP coefficient file

  std::string buffer = *path + "/" + coefficient_file;

  std::fstream fs;

  fs.open(buffer.c_str(), std::fstream::out);

  if (!fs.is_open())
  {
    HELPER_LOG_ERROR("Unable to open the coefficient file for writing.");
    return true;
  }
  else
  {
    if (fs.fail())
    {
      HELPER_LOG_ERROR("An error has occurred from opening the file on the associated stream!");
      return true;
    }

    fs << "# SNAP coefficients for ";
    for (int ielem = 0; ielem < nelements; ++ielem)
    {
      fs << elements[ielem] << " ";
    }
    fs << "\n"
       << std::endl;

    // nelements & ncoeffall
    fs << nelements << " " << ncoeffall << std::endl;

    fs << std::setprecision(std::numeric_limits<double>::digits10 + 1) << std::fixed;

    // Loop over nelements blocks in the SNAP coefficient file
    for (int ielem = 0; ielem < nelements; ++ielem)
    {
      fs << elements[ielem] << " " << radelem[ielem] << " " << wjelem[ielem] << std::endl;

      for (int icoeff = 0; icoeff < ncoeffall; ++icoeff)
      {
        fs << coeffelem(ielem, icoeff) << std::endl;
      }
    }

    fs << "\n\n"
       << "#\n"
       << "# SNAP coefficient file"
       << "#\n"
       << "# First line: nelements & ncoeffall\n"
       << "#\n"
       << "# This SNAP coefficient file has " << nelements << ((nelements > 1) ? " blocks, where \n" : " block, where\n")
       << "#    First line of each block contains: 'species name' & 'species radius' & 'species weight' \n"
       << "#    followed by " << ncoeffall << " SNAP coefficients. \n"
       << "#\n";

    fs.close();
  }

  // SNAP parameter file

  buffer = *path + "/" + parameter_file;

  fs.open(buffer.c_str(), std::fstream::out);

  if (!fs.is_open())
  {
    HELPER_LOG_ERROR("Unable to open the parameter file for writing.");
    return true;
  }
  else
  {
    if (fs.fail())
    {
      HELPER_LOG_ERROR("An error has occurred from opening the file on the associated stream!");
      return true;
    }

    fs << "# SNAP parameters for ";
    for (int ielem = 0; ielem < nelements; ++ielem)
    {
      fs << elements[ielem] << " ";
    }
    fs << "\n";
    fs << std::endl;

    fs << std::setprecision(std::numeric_limits<double>::digits10 + 1) << std::fixed;

    // rcutfac
    fs << "rcutfac       " << rcutfac << std::endl;

    // twojmax
    fs << "twojmax       " << twojmax << std::endl;

    // rfac0
    fs << "rfac0         " << rfac0 << std::endl;

    // rmin0
    fs << "rmin0         " << rmin0 << std::endl;

    // switchflag
    fs << "switchflag    " << switchflag << std::endl;

    // bzeroflag
    fs << "bzeroflag     " << bzeroflag << std::endl;

    // quadraticflag
    fs << "quadraticflag " << quadraticflag << std::endl;

    fs << "\n\n"
       << "#\n"
       << "# SNAP parameter file"
       << "#\n";

    fs.close();
  }

  return false;
}

int SNAPImplementation::Compute(KIM::ModelCompute const *const modelCompute,
                                KIM::ModelComputeArguments const *const modelComputeArguments)
{
  // KIM API Model Input compute flags
  bool isComputeProcess_dEdr = false;
  bool isComputeProcess_d2Edr2 = false;

  // KIM API Model Output compute flags
  bool isComputeEnergy = false;
  bool isComputeForces = false;
  bool isComputeParticleEnergy = false;
  bool isComputeVirial = false;
  bool isComputeParticleVirial = false;

  // KIM API Model Input
  int const *particleSpeciesCodes = NULL;
  int const *particleContributing = NULL;

  VectorOfSizeDIM const *coordinates = NULL;

  // KIM API Model Output
  double *energy = NULL;
  VectorOfSizeDIM *forces = NULL;
  double *particleEnergy = NULL;
  VectorOfSizeSix *virial = NULL;
  VectorOfSizeSix *particleVirial = NULL;

  int ier = setComputeMutableValues(modelComputeArguments,
                                    isComputeProcess_dEdr,
                                    isComputeProcess_d2Edr2,
                                    isComputeEnergy,
                                    isComputeForces,
                                    isComputeParticleEnergy,
                                    isComputeVirial,
                                    isComputeParticleVirial,
                                    particleSpeciesCodes,
                                    particleContributing,
                                    coordinates,
                                    energy,
                                    forces,
                                    particleEnergy,
                                    virial,
                                    particleVirial);
  if (ier)
  {
    HELPER_LOG_ERROR("setComputeMutableValues fails.");
    return ier;
  }

  if (beta_max < numberOfContributingParticles_)
  {
    beta_max = numberOfContributingParticles_;

    beta.resize(beta_max, ncoeff);
    bispectrum.resize(beta_max, ncoeff);
  }

  // Compute the bispectrum component of atom
  if (quadraticflag || isComputeEnergy || isComputeParticleEnergy)
  {
    computeBispectrum(modelComputeArguments, particleSpeciesCodes, particleContributing, coordinates);
  }

  // Compute dE_i/dB_i = beta_i for all i in list
  computeBeta(particleSpeciesCodes, particleContributing);

#include "SNAPImplementationComputeDispatch.cpp"

  return ier;
}

int SNAPImplementation::ComputeArgumentsCreate(KIM::ModelComputeArgumentsCreate *const modelComputeArgumentsCreate) const
{
  int ier = RegisterKIMComputeArgumentsSettings(modelComputeArgumentsCreate);
  if (ier)
  {
    return ier;
  }
  // Nothing else to do for this case
  // Everything is good
  return false;
}

int SNAPImplementation::ComputeArgumentsDestroy(KIM::ModelComputeArgumentsDestroy *const modelComputeArgumentsDestroy) const
{
  // Avoid not used warning
  (void)modelComputeArgumentsDestroy;
  // Nothing else to do for this case
  // Everything is good
  return false;
}

// Private member functions

#define KIM_LOGGER_OBJECT_NAME modelDriverCreate
int SNAPImplementation::OpenParameterFiles(KIM::ModelDriverCreate *const modelDriverCreate,
                                           int const numberParameterFiles,
                                           std::FILE **parameterFilePointers)
{
  for (int i = 0; i < numberParameterFiles; ++i)
  {
    std::string const *parameterFileName;

    int ier = modelDriverCreate->GetParameterFileName(i, &parameterFileName);
    if (ier)
    {
      LOG_ERROR("Unable to get parameter file name \n");
      return ier;
    }

    parameterFilePointers[i] = std::fopen(parameterFileName->c_str(), "r");
    if (!parameterFilePointers[i])
    {
      HELPER_LOG_ERROR("SNAP parameter file number " + std::to_string(i) + " can not be opened \n");
      for (int j = i - 1; i <= 0; --i)
      {
        std::fclose(parameterFilePointers[j]);
      }
      return true;
    }
  }

  // everything is good
  return false;
}
#undef KIM_LOGGER_OBJECT_NAME

void SNAPImplementation::GetNextDataLine(std::FILE *const filePtr,
                                         char *nextLinePtr,
                                         int const maxSize,
                                         int *endOfFileFlag)
{
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
    *pch = '\0';
}

#define KIM_LOGGER_OBJECT_NAME modelDriverCreate
int SNAPImplementation::ProcessParameterFiles(KIM::ModelDriverCreate *const modelDriverCreate,
                                              int const /* numberParameterFiles */,
                                              std::FILE *const *parameterFilePointers)
{
  char nextLine[MAXLINE];

  int endOfFileFlag(0);

  // Read SNAP coefficient file
  GetNextDataLine(parameterFilePointers[0], nextLine, MAXLINE, &endOfFileFlag);
  if (endOfFileFlag)
  {
    HELPER_LOG_ERROR("End of file in SNAP coefficient file \n");
    return true;
  }

  int ier = std::sscanf(nextLine, "%d %d", &nelements, &ncoeffall);
  if (ier != 2)
  {
    std::string msg = "unable to read nelements & ncoeffall from the line: \n";
    msg += nextLine;
    msg += "\n";
    HELPER_LOG_ERROR(msg);
    return true;
  }

  // Allocate memory based on the number of elements

  // Set up memory for the element lists
  elements.resize(nelements);
  radelem.resize(nelements);
  wjelem.resize(nelements);
  coeffelem.resize(nelements, ncoeffall);

  // Loop over nelements blocks in the SNAP coefficient file
  for (int ielem = 0; ielem < nelements; ++ielem)
  {
    GetNextDataLine(parameterFilePointers[0], nextLine, MAXLINE, &endOfFileFlag);
    if (endOfFileFlag)
    {
      HELPER_LOG_ERROR("End of file in SNAP coefficient file \n");
      return true;
    }

    char spec[MAXLINE];
    ier = std::sscanf(nextLine, "%s %lf %lf", spec, &radelem[ielem], &wjelem[ielem]);
    if (ier != 3)
    {
      std::string msg = "Incorrect format in SNAP coefficient file \n";
      msg += "Unable to read spec, radius & weight from the line: \n";
      msg += nextLine;
      msg += "\n";
      HELPER_LOG_ERROR(msg);
      return true;
    }

    KIM::SpeciesName const speciesName(spec);

    if (speciesName.ToString() == "unknown")
    {
      std::string msg = "Incorrect format in SNAP coefficient file \n";
      msg += "Input species of '";
      msg += spec;
      msg += "' is unknown. \n";
      HELPER_LOG_ERROR(msg);
      return true;
    }

    for (int i = 0; i < ielem; ++i)
    {
      KIM::SpeciesName const tmpName(elements[i]);
      if (speciesName == tmpName)
      {
        std::string msg = "Incorrect format in SNAP coefficient file \n";
        msg += "Species '";
        msg += spec;
        msg += "' is already defined and exist. \n";
        HELPER_LOG_ERROR(msg);
        return true;
      }
    }

    elements[ielem] = speciesName.ToString();

    ier = modelDriverCreate->SetSpeciesCode(speciesName, ielem);
    if (ier)
    {
      LOG_ERROR("SetSpeciesCode failed to set the new species. \n");
      return ier;
    }

    for (int icoeff = 0; icoeff < ncoeffall; ++icoeff)
    {
      GetNextDataLine(parameterFilePointers[0], nextLine, MAXLINE, &endOfFileFlag);
      if (endOfFileFlag)
      {
        HELPER_LOG_ERROR("End of file in SNAP coefficient file. \n");
        return true;
      }

      double coeff;
      ier = std::sscanf(nextLine, "%lf", &coeff);
      if (ier != 1)
      {
        HELPER_LOG_ERROR("Incorrect format in SNAP coefficient file. \n");
        return true;
      }
      coeffelem(ielem, icoeff) = coeff;
    }
  }

  // set flags for required keywords
  bool rcutfacflag(false);
  bool twojmaxflag(false);

  // Read SNAP parameter file
  while (1)
  {
    GetNextDataLine(parameterFilePointers[1], nextLine, MAXLINE, &endOfFileFlag);
    if (endOfFileFlag)
      break;

    // words = ptrs to all words in line
    // strip single and double quotes from words
    char *keywd = std::strtok(nextLine, "' \t\n\r\f");
    char *keyval = std::strtok(NULL, "' \t\n\r\f");

    if (!std::strcmp(keywd, "rcutfac"))
    {
      rcutfac = std::atof(keyval);
      rcutfacflag = true;
    }
    else if (!std::strcmp(keywd, "twojmax"))
    {
      twojmax = std::atoi(keyval);
      twojmaxflag = true;
    }
    else if (!std::strcmp(keywd, "rfac0"))
    {
      rfac0 = std::atof(keyval);
    }
    else if (!std::strcmp(keywd, "rmin0"))
    {
      rmin0 = std::atof(keyval);
    }
    else if (!std::strcmp(keywd, "switchflag"))
    {
      switchflag = std::atoi(keyval);
      if (switchflag != 0 && switchflag != 1)
        switchflag = 1;
    }
    else if (!std::strcmp(keywd, "bzeroflag"))
    {
      bzeroflag = std::atoi(keyval);
      if (bzeroflag != 0 && bzeroflag != 1)
        bzeroflag = 1;
    }
    else if (!std::strcmp(keywd, "quadraticflag"))
    {
      quadraticflag = std::atoi(keyval);
      if (quadraticflag != 0 && quadraticflag != 1)
        quadraticflag = 1;
    }
    else
    {
      std::string msg = "Incorrect SNAP parameter file \n";
      msg += "Wrong SNAP keywords of ";
      msg += keywd;
      msg += ", ";
      msg += keyval;
      msg += "\n";
      HELPER_LOG_ERROR(msg);
      return true;
    }
  }

  if (!rcutfacflag || !twojmaxflag)
  {
    std::string msg = "Incorrect SNAP parameter file \n";
    msg += rcutfacflag ? "'twojmaxflag' is not given \n" : "'rcutfacflag' is not given \n";
    HELPER_LOG_ERROR(msg);
    return true;
  }

  // everything is good
  return false;
}
#undef KIM_LOGGER_OBJECT_NAME

void SNAPImplementation::CloseParameterFiles(int const numberParameterFiles,
                                             std::FILE *const *parameterFilePointers)
{
  for (int i = 0; i < numberParameterFiles; ++i)
    std::fclose(parameterFilePointers[i]);
}

#define KIM_LOGGER_OBJECT_NAME modelDriverCreate
int SNAPImplementation::ConvertUnits(KIM::ModelDriverCreate *const modelDriverCreate,
                                     KIM::LengthUnit const &requestedLengthUnit,
                                     KIM::EnergyUnit const &requestedEnergyUnit,
                                     KIM::ChargeUnit const &requestedChargeUnit,
                                     KIM::TemperatureUnit const &requestedTemperatureUnit,
                                     KIM::TimeUnit const &requestedTimeUnit)
{
  // Define metal unit as a default base units
  KIM::LengthUnit const fromLength = KIM::LENGTH_UNIT::A;
  KIM::EnergyUnit const fromEnergy = KIM::ENERGY_UNIT::eV;
  KIM::ChargeUnit const fromCharge = KIM::CHARGE_UNIT::e;
  KIM::TemperatureUnit const fromTemperature = KIM::TEMPERATURE_UNIT::K;
  KIM::TimeUnit const fromTime = KIM::TIME_UNIT::ps;

  // changing units of sigma, gamma, and cutoff
  double convertLength = ONE;
  int ier = modelDriverCreate->ConvertUnit(fromLength,
                                           fromEnergy,
                                           fromCharge,
                                           fromTemperature,
                                           fromTime,
                                           requestedLengthUnit,
                                           requestedEnergyUnit,
                                           requestedChargeUnit,
                                           requestedTemperatureUnit,
                                           requestedTimeUnit,
                                           1.0,
                                           0.0,
                                           0.0,
                                           0.0,
                                           0.0,
                                           &convertLength);
  if (ier)
  {
    LOG_ERROR("Unable to convert unit");
    return ier;
  }

  // convert to active units
  if (convertLength != ONE)
  {
    std::for_each(radelem.begin(), radelem.end(), [&](double &rc) { rc *= convertLength; });
    rmin0 *= convertLength;
  }

  // changing units of A and lambda
  double convertEnergy = ONE;
  ier = modelDriverCreate->ConvertUnit(fromLength,
                                       fromEnergy,
                                       fromCharge,
                                       fromTemperature,
                                       fromTime,
                                       requestedLengthUnit,
                                       requestedEnergyUnit,
                                       requestedChargeUnit,
                                       requestedTemperatureUnit,
                                       requestedTimeUnit,
                                       0.0,
                                       1.0,
                                       0.0,
                                       0.0,
                                       0.0,
                                       &convertEnergy);
  if (ier)
  {
    LOG_ERROR("Unable to convert energy unit");
    return ier;
  }

  // convert to active units
  if (convertEnergy != ONE)
  {
    // Loop over nelements blocks in the SNAP coefficient file
    for (int ielem = 0; ielem < nelements; ++ielem)
    {
      for (int icoeff = 0; icoeff < ncoeffall; ++icoeff)
      {
        coeffelem(ielem, icoeff) *= convertEnergy;
      }
    }
  }

  // Register units
  ier = modelDriverCreate->SetUnits(requestedLengthUnit,
                                    requestedEnergyUnit,
                                    KIM::CHARGE_UNIT::unused,
                                    KIM::TEMPERATURE_UNIT::unused,
                                    KIM::TIME_UNIT::unused);
  if (ier)
  {
    LOG_ERROR("Unable to set units to requested values");
    return ier;
  }

  // Everything is good
  return false;
}
#undef KIM_LOGGER_OBJECT_NAME

template <class ModelObj>
int SNAPImplementation::setRefreshMutableValues(ModelObj *const modelObj)
{
  // Correct the total number of coefficients
  if (quadraticflag)
  {
    ncoeff = std::sqrt(2 * ncoeffall) - 1;

    int const ntmp = 1 + ncoeff * (ncoeff + 3) / 2;

    if (ntmp != ncoeffall)
    {
      std::string msg = "Incorrect args. \n";
      msg += "ncoeffall = " + std::to_string(ncoeffall);
      msg += " ntmp = " + std::to_string(ntmp);
      msg += " ncoeff = " + std::to_string(ncoeff);
      msg += "\n";
      HELPER_LOG_ERROR(msg);
      return true;
    }
  }
  else
  {
    ncoeff = ncoeffall - 1;
  }

  // We have to creat the SNAP object and do the extra check on
  // the number of coefficients

  // Construct the SNAP object
  snap.reset(new SNA(rfac0, twojmax, rmin0, switchflag, bzeroflag));

  // Extra check
  if (ncoeff != snap->ncoeff)
  {
    std::string msg = "Wrong number of coefficients \n";
    msg += "Number of coefficients to the model and the one created in SNAP object do not match \n";
    msg += "ncoeff = " + std::to_string(ncoeff);
    msg += " & SNAP ncoeff = " + std::to_string(snap->ncoeff);
    msg += "\n";
    HELPER_LOG_ERROR(msg);
    return true;
  }

  // Set up memory for square of cutoff
  cutsq.resize(nelements, nelements);

  // Calculate maximum cutoff for all elements & update parameters
  double rcutmax(0.0);

  for (int ielem = 0; ielem < nelements; ++ielem)
  {
    double const radi = radelem[ielem];

    rcutmax = std::max(2.0 * radi * rcutfac, rcutmax);

    for (int jelem = 0; jelem <= ielem; ++jelem)
    {
      double const cut = (radi + radelem[jelem]) * rcutfac;

      // Set the square of cutoff for all pairs
      cutsq(ielem, jelem) = cutsq(jelem, ielem) = cut * cut;
    }
  }

  influenceDistance_ = rcutmax;

  // Update the maximum cutoff or influence distance value in KIM API object
  modelObj->SetInfluenceDistancePointer(&influenceDistance_);

  // In this model we do not need to request neighbors from non-contributing particles
  // modelWillNotRequestNeighborsOfNoncontributingParticles_ = 1
  // Update the cutoff value in KIM API object
  modelObj->SetNeighborListPointers(1, &influenceDistance_, &modelWillNotRequestNeighborsOfNoncontributingParticles_);

  // everything is good
  return false;
}

int SNAPImplementation::RegisterKIMFunctions(KIM::ModelDriverCreate *const modelDriverCreate) const
{
  // Register functions
  int ier =
      modelDriverCreate->SetRoutinePointer(KIM::MODEL_ROUTINE_NAME::Destroy, KIM::LANGUAGE_NAME::cpp, true, reinterpret_cast<KIM::Function *>(SNAP::Destroy)) ||
      modelDriverCreate->SetRoutinePointer(KIM::MODEL_ROUTINE_NAME::Refresh, KIM::LANGUAGE_NAME::cpp, true, reinterpret_cast<KIM::Function *>(SNAP::Refresh)) ||
      modelDriverCreate->SetRoutinePointer(KIM::MODEL_ROUTINE_NAME::WriteParameterizedModel, KIM::LANGUAGE_NAME::cpp, false, reinterpret_cast<KIM::Function *>(SNAP::WriteParameterizedModel)) ||
      modelDriverCreate->SetRoutinePointer(KIM::MODEL_ROUTINE_NAME::Compute, KIM::LANGUAGE_NAME::cpp, true, reinterpret_cast<KIM::Function *>(SNAP::Compute)) ||
      modelDriverCreate->SetRoutinePointer(KIM::MODEL_ROUTINE_NAME::ComputeArgumentsCreate, KIM::LANGUAGE_NAME::cpp, true, reinterpret_cast<KIM::Function *>(SNAP::ComputeArgumentsCreate)) ||
      modelDriverCreate->SetRoutinePointer(KIM::MODEL_ROUTINE_NAME::ComputeArgumentsDestroy, KIM::LANGUAGE_NAME::cpp, true, reinterpret_cast<KIM::Function *>(SNAP::ComputeArgumentsDestroy));
  return ier;
}

#define KIM_LOGGER_OBJECT_NAME modelComputeArguments
int SNAPImplementation::setComputeMutableValues(KIM::ModelComputeArguments const *const modelComputeArguments,
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
                                                VectorOfSizeSix *&particleVirial)
{
  // get compute flags
  int compProcess_dEdr;
  int compProcess_d2Edr2;

  modelComputeArguments->IsCallbackPresent(KIM::COMPUTE_CALLBACK_NAME::ProcessDEDrTerm, &compProcess_dEdr);
  modelComputeArguments->IsCallbackPresent(KIM::COMPUTE_CALLBACK_NAME::ProcessD2EDr2Term, &compProcess_d2Edr2);

  isComputeProcess_dEdr = compProcess_dEdr;
  isComputeProcess_d2Edr2 = compProcess_d2Edr2;

  int const *numberOfParticles = NULL;
  int ier = modelComputeArguments->GetArgumentPointer(KIM::COMPUTE_ARGUMENT_NAME::numberOfParticles, &numberOfParticles) ||
            modelComputeArguments->GetArgumentPointer(KIM::COMPUTE_ARGUMENT_NAME::particleSpeciesCodes, &particleSpeciesCodes) ||
            modelComputeArguments->GetArgumentPointer(KIM::COMPUTE_ARGUMENT_NAME::particleContributing, &particleContributing) ||
            modelComputeArguments->GetArgumentPointer(KIM::COMPUTE_ARGUMENT_NAME::coordinates, (double const **)&coordinates) ||
            modelComputeArguments->GetArgumentPointer(KIM::COMPUTE_ARGUMENT_NAME::partialEnergy, &energy) ||
            modelComputeArguments->GetArgumentPointer(KIM::COMPUTE_ARGUMENT_NAME::partialForces, (double const **)&forces) ||
            modelComputeArguments->GetArgumentPointer(KIM::COMPUTE_ARGUMENT_NAME::partialParticleEnergy, &particleEnergy) ||
            modelComputeArguments->GetArgumentPointer(KIM::COMPUTE_ARGUMENT_NAME::partialVirial, (double const **)&virial) ||
            modelComputeArguments->GetArgumentPointer(KIM::COMPUTE_ARGUMENT_NAME::partialParticleVirial, (double const **)&particleVirial);
  if (ier)
  {
    LOG_ERROR("GetArgumentPointer");
    return true;
  }

  isComputeEnergy = (energy);
  isComputeForces = (forces);
  isComputeParticleEnergy = (particleEnergy);
  isComputeVirial = (virial);
  isComputeParticleVirial = (particleVirial);

  // Update values
  cachedNumberOfParticles_ = *numberOfParticles;

  // Get the number of contributing particles
  numberOfContributingParticles_ = std::accumulate(particleContributing, particleContributing + cachedNumberOfParticles_, 0);

  // everything is good
  return false;
}
#undef KIM_LOGGER_OBJECT_NAME

int SNAPImplementation::getComputeIndex(bool const isComputeProcess_dEdr,
                                        bool const isComputeProcess_d2Edr2,
                                        bool const isComputeEnergy,
                                        bool const isComputeForces,
                                        bool const isComputeParticleEnergy,
                                        bool const isComputeVirial,
                                        bool const isComputeParticleVirial) const
{
  int const processd2E = 2;
  int const energy = 2;
  int const force = 2;
  int const particleEnergy = 2;
  int const virial = 2;
  int const particleVirial = 2;

  int index = 0;
  // processdE
  index += (int(isComputeProcess_dEdr)) * processd2E * energy * force * particleEnergy * virial * particleVirial;
  // processd2E
  index += (int(isComputeProcess_d2Edr2)) * energy * force * particleEnergy * virial * particleVirial;
  // energy
  index += (int(isComputeEnergy)) * force * particleEnergy * virial * particleVirial;
  // force
  index += (int(isComputeForces)) * particleEnergy * virial * particleVirial;
  // particleEnergy
  index += (int(isComputeParticleEnergy)) * virial * particleVirial;
  // virial
  index += (int(isComputeVirial)) * particleVirial;
  // particleVirial
  index += (int(isComputeParticleVirial));
  return index;
}

int SNAPImplementation::RegisterKIMModelSettings(KIM::ModelDriverCreate *const modelDriverCreate) const
{
  // Register numbering
  return modelDriverCreate->SetModelNumbering(KIM::NUMBERING::zeroBased);
}

#define KIM_LOGGER_OBJECT_NAME modelDriverCreate
int SNAPImplementation::RegisterKIMParameters(KIM::ModelDriverCreate *const modelDriverCreate)
{
  // Publish parameters
  int ier;

  ier = modelDriverCreate->SetParameterPointer(nelements, radelem.data(), "R",
                                               "List of cutoff radii, one for each type (distance units) 'radelem'.");
  if (ier)
  {
    LOG_ERROR("SetParameterPointer radelem");
    return ier;
  }

  ier = modelDriverCreate->SetParameterPointer(nelements, wjelem.data(), "W",
                                               "List of neighbor weights, one for each type 'wjelem'.");
  if (ier)
  {
    LOG_ERROR("SetParameterPointer wjelem");
    return ier;
  }

  ier = modelDriverCreate->SetParameterPointer(nelements * ncoeffall, coeffelem.data(),
                                               "coeffs",
                                               "Element bispectrum coefficients : "
                                               "(matrix of size M x N = " +
                                                   std::to_string(nelements) +
                                                   " x " +
                                                   std::to_string(ncoeffall) +
                                                   ") " +
                                                   "in row-major storage. Ordering is according to each element. "
                                                   "For example, to find the coefficient 'j' related to element 'i', "
                                                   "use index = ((i - 1) * nelements + j). 'coeffelem'.");
  if (ier)
  {
    LOG_ERROR("SetParameterPointer coeffelem");
    return ier;
  }

  ier = modelDriverCreate->SetParameterPointer(1, &rcutfac, "rcutfac",
                                               "Scale factor applied to all cutoff radii (positive real) 'rcutfac'.");
  if (ier)
  {
    LOG_ERROR("SetParameterPointer rcutfac");
    return ier;
  }

  ier = modelDriverCreate->SetParameterPointer(1, &twojmax, "twojmax",
                                               "Band limit for bispectrum components (non-negative integer) 'twojmax'.");
  if (ier)
  {
    LOG_ERROR("SetParameterPointer twojmax");
    return ier;
  }

  ier = modelDriverCreate->SetParameterPointer(1, &rfac0, "rfac0",
                                               "Parameter in distance to angle conversion (0 < rcutfac < 1) 'rfac0'.");
  if (ier)
  {
    LOG_ERROR("SetParameterPointer rfac0");
    return ier;
  }

  ier = modelDriverCreate->SetParameterPointer(1, &rmin0, "rmin0",
                                               "Parameter in distance to angle conversion (distance units) 'rmin0'.");
  if (ier)
  {
    LOG_ERROR("SetParameterPointer rmin0");
    return ier;
  }

  ier = modelDriverCreate->SetParameterPointer(1, &switchflag, "switchflag",
                                               "0 or 1 \n "
                                               "0 = do not use switching function, \n "
                                               "1 = use switching function, 'switchflag'.");
  if (ier)
  {
    LOG_ERROR("SetParameterPointer switchflag");
    return ier;
  }

  ier = modelDriverCreate->SetParameterPointer(1, &bzeroflag, "bzeroflag",
                                               "0 or 1 \n "
                                               "0 = do not subtract B0, \n "
                                               "1 = subtract B0, 'bzeroflag'.");
  if (ier)
  {
    LOG_ERROR("SetParameterPointer bzeroflag");
    return ier;
  }

  ier = modelDriverCreate->SetParameterPointer(1, &quadraticflag, "quadraticflag",
                                               "0 or 1 \n "
                                               "0 = do not generate quadratic terms, \n "
                                               "1 = generate quadratic terms, 'quadraticflag'.");
  if (ier)
  {
    LOG_ERROR("SetParameterPointer quadraticflag");
    return ier;
  }

  // Everything is good
  return false;
}
#undef KIM_LOGGER_OBJECT_NAME

#define KIM_LOGGER_OBJECT_NAME modelComputeArgumentsCreate
int SNAPImplementation::RegisterKIMComputeArgumentsSettings(KIM::ModelComputeArgumentsCreate *const modelComputeArgumentsCreate) const
{
  // register arguments
  LOG_INFORMATION("Register argument supportStatus");

  int ier =
      modelComputeArgumentsCreate->SetArgumentSupportStatus(KIM::COMPUTE_ARGUMENT_NAME::partialEnergy, KIM::SUPPORT_STATUS::optional) ||
      modelComputeArgumentsCreate->SetArgumentSupportStatus(KIM::COMPUTE_ARGUMENT_NAME::partialForces, KIM::SUPPORT_STATUS::optional) ||
      modelComputeArgumentsCreate->SetArgumentSupportStatus(KIM::COMPUTE_ARGUMENT_NAME::partialParticleEnergy, KIM::SUPPORT_STATUS::optional) ||
      modelComputeArgumentsCreate->SetArgumentSupportStatus(KIM::COMPUTE_ARGUMENT_NAME::partialVirial, KIM::SUPPORT_STATUS::optional) ||
      modelComputeArgumentsCreate->SetArgumentSupportStatus(KIM::COMPUTE_ARGUMENT_NAME::partialParticleVirial, KIM::SUPPORT_STATUS::optional);

  // register callbacks
  LOG_INFORMATION("Register callback supportStatus");

  ier = ier ||
        modelComputeArgumentsCreate->SetCallbackSupportStatus(KIM::COMPUTE_CALLBACK_NAME::ProcessDEDrTerm, KIM::SUPPORT_STATUS::optional) ||
        modelComputeArgumentsCreate->SetCallbackSupportStatus(KIM::COMPUTE_CALLBACK_NAME::ProcessD2EDr2Term, KIM::SUPPORT_STATUS::optional);

  return ier;
}
#undef KIM_LOGGER_OBJECT_NAME

void SNAPImplementation::computeBispectrum(KIM::ModelComputeArguments const *const modelComputeArguments,
                                           int const *particleSpeciesCodes,
                                           int const *particleContributing,
                                           VectorOfSizeDIM const *coordinates)
{
  // calculate contribution from pair function
  int numnei = 0;
  int const *n1atom = NULL;

  // Setup loop over contributing particles
  for (int i = 0, contributing_index = 0; i < cachedNumberOfParticles_; ++i)
  {
    if (!particleContributing[i])
      continue;

    // Get neighbors of i
    modelComputeArguments->GetNeighborList(0, i, &numnei, &n1atom);

    // Get the species index for atom i
    int const iSpecies = particleSpeciesCodes[i];
    double const radi = radelem[iSpecies];

    double const x = coordinates[i][0];
    double const y = coordinates[i][1];
    double const z = coordinates[i][2];

    // Make sure {rij, inside, wj, and rcutij} arrays are big enough for numnei atoms
    snap->grow_rij(numnei);

    // number of neighbors of I within cutoff
    int ninside = 0;

    // note Rij sign convention => dU/dRij = dU/dRj = -dU/dRi

    // Setup loop over neighbors of current particle
    for (int n = 0; n < numnei; ++n)
    {
      int const j = n1atom[n];
      int const jSpecies = particleSpeciesCodes[j];

      double const dx = coordinates[j][0] - x;
      double const dy = coordinates[j][1] - y;
      double const dz = coordinates[j][2] - z;

      double const rsq = dx * dx + dy * dy + dz * dz;

      if (rsq < cutsq(iSpecies, jSpecies) && rsq > 1e-20)
      {
        snap->rij(ninside, 0) = dx;
        snap->rij(ninside, 1) = dy;
        snap->rij(ninside, 2) = dz;
        snap->inside[ninside] = j;
        snap->wj[ninside] = wjelem[jSpecies];
        snap->rcutij[ninside] = (radi + radelem[jSpecies]) * rcutfac;
        ++ninside;
      }
    }

    snap->compute_ui(ninside);

    snap->compute_zi();

    snap->compute_bi();

    auto Bi = bispectrum.data_1D(contributing_index++);
    for (int icoeff = 0; icoeff < ncoeff; ++icoeff)
      Bi[icoeff] = snap->blist[icoeff];
  }
}

void SNAPImplementation::computeBeta(int const *particleSpeciesCodes, int const *particleContributing)
{
  // Setup loop over contributing particles
  if (quadraticflag)
  {
    for (int i = 0, contributing_index = 0; i < cachedNumberOfParticles_; ++i)
    {
      if (!particleContributing[i])
        continue;

      // Get the species index for atom i
      int const iSpecies = particleSpeciesCodes[i];

      // Get the 1D view to the 2D coeffelem array at row iSpecies
      auto coeffi = coeffelem.data_1D(iSpecies);

      // Get the pointer to the raw data + 1 to avoid extra sum
      double *Ci = coeffi.data() + 1;

      // Get the 1D view to the 2D beta array at row i
      auto bi = beta.data_1D(contributing_index);

      for (int icoeff = 0; icoeff < ncoeff; ++icoeff)
        bi[icoeff] = Ci[icoeff];

      // Get the pointer to the start of coeffi array of data
      --Ci;

      // Get the 1D view to the 2D bispectrum array at row i
      auto Bi = bispectrum.data_1D(contributing_index++);

      int k = ncoeff + 1;

      for (int icoeff = 0; icoeff < ncoeff; ++icoeff)
      {
        double const bveci = Bi[icoeff];

        bi[icoeff] += Ci[k++] * bveci;

        for (int jcoeff = icoeff + 1; jcoeff < ncoeff; ++jcoeff, ++k)
        {
          bi[icoeff] += Ci[k] * Bi[jcoeff];
          bi[jcoeff] += Ci[k] * bveci;
        }
      }
    }
  }
  else
  {
    for (int i = 0, contributing_index = 0; i < cachedNumberOfParticles_; ++i)
    {
      if (!particleContributing[i])
        continue;

      // Get the species index for atom i
      int const iSpecies = particleSpeciesCodes[i];

      // Get the 1D view to the 2D coeffelem array at row iSpecies
      auto coeffi = coeffelem.data_1D(iSpecies);

      // Get the pointer to the raw data + 1 to avoid extra sum
      double *Ci = coeffi.data() + 1;

      // Get the 1D view to the 2D beta array at row i
      auto bi = beta.data_1D(contributing_index++);

      for (int icoeff = 0; icoeff < ncoeff; ++icoeff)
        bi[icoeff] = Ci[icoeff];
    }
  }
}

#define KIM_LOGGER_OBJECT_NAME modelComputeArguments
template <bool isComputeProcess_dEdr,
          bool isComputeProcess_d2Edr2,
          bool isComputeEnergy,
          bool isComputeForces,
          bool isComputeParticleEnergy,
          bool isComputeVirial,
          bool isComputeParticleVirial>
int SNAPImplementation::Compute(KIM::ModelCompute const *const /* modelCompute */,
                                KIM::ModelComputeArguments const *const modelComputeArguments,
                                int const *const particleSpeciesCodes,
                                int const *const particleContributing,
                                const VectorOfSizeDIM *const coordinates,
                                double *const energy,
                                VectorOfSizeDIM *const forces,
                                double *const particleEnergy,
                                VectorOfSizeSix virial,
                                VectorOfSizeSix *const particleVirial) const
{
  // Everything is good
  int ier = false;

  // Nothing to compute
  if (!isComputeEnergy &&
      !isComputeForces &&
      !isComputeParticleEnergy &&
      !isComputeVirial &&
      !isComputeParticleVirial &&
      !isComputeProcess_dEdr &&
      !isComputeProcess_d2Edr2)
    return ier;

  // Initialize energy
  if (isComputeEnergy)
  {
    *energy = 0.0;
  }

  // Initialize forces
  if (isComputeForces)
  {
    for (int i = 0; i < cachedNumberOfParticles_; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        forces[i][j] = 0.0;
      }
    }
  }

  // Initialize particle energy
  if (isComputeParticleEnergy)
  {
    for (int i = 0; i < cachedNumberOfParticles_; ++i)
    {
      particleEnergy[i] = 0.0;
    }
  }

  // Initialize virial
  if (isComputeVirial)
  {
    for (int i = 0; i < 6; ++i)
    {
      virial[i] = 0.0;
    }
  }

  // Initialize particle virial
  if (isComputeParticleVirial)
  {
    for (int i = 0; i < cachedNumberOfParticles_; ++i)
    {
      for (int j = 0; j < 6; ++j)
      {
        particleVirial[i][j] = 0.0;
      }
    }
  }

  int numnei = 0;
  int const *n1atom = NULL;

  // Loop over all the contributing particles
  for (int i = 0, contributing_index = 0; i < cachedNumberOfParticles_; ++i)
  {
    if (!particleContributing[i])
      continue;

    // Get the species index for atom i
    int const iSpecies = particleSpeciesCodes[i];

    // Get the iSpecies cutoff
    double const radi = radelem[iSpecies];

    double const xi = coordinates[i][0];
    double const yi = coordinates[i][1];
    double const zi = coordinates[i][2];

    // Calculate contribution

    // Get the neighbors
    modelComputeArguments->GetNeighborList(0, i, &numnei, &n1atom);

    // Make sure {rij, inside, wj, and rcutij} arrays are big enough for numnei atoms
    snap->grow_rij(numnei);

    // number of neighbors of I within cutoff
    int ninside = 0;

    // note Rij sign convention => dU/dRij = dU/dRj = -dU/dRi

    // Setup loop over neighbors of current particle
    for (int n = 0; n < numnei; ++n)
    {
      // Index of the neighbor atom
      int const j = n1atom[n];

      // Get the species index for atom j
      int const jSpecies = particleSpeciesCodes[j];

      double const dx = coordinates[j][0] - xi;
      double const dy = coordinates[j][1] - yi;
      double const dz = coordinates[j][2] - zi;

      double const rsq = dx * dx + dy * dy + dz * dz;

      if (rsq < cutsq(iSpecies, jSpecies) && rsq > 1e-20)
      {
        snap->rij(ninside, 0) = dx;
        snap->rij(ninside, 1) = dy;
        snap->rij(ninside, 2) = dz;
        snap->inside[ninside] = j;
        snap->wj[ninside] = wjelem[jSpecies];
        // Get the pair of iSpecies & jSpecies cutoff
        snap->rcutij[ninside] = (radi + radelem[jSpecies]) * rcutfac;
        ++ninside;
      }
    }

    // compute Ui, Yi for atom i
    {
      snap->compute_ui(ninside);

      // Get the 1D view to the 2D beta array at row i
      auto betai = beta.data_1D(contributing_index);

      // Get the pointer to the beta array of data for atom i
      double const *const bi = const_cast<double *>(betai.data());

      snap->compute_yi(bi);
    }

    // Compute contribution to force, etc.

    // For the neighbors of particle i within the cutoff:
    // Compute deidrj = dEi/dRj = -dEi/dRi add to Fi, subtract from Fj

    VectorOfSizeDIM deidrj;

    // Setup loop over neighbors of particle i
    for (int n = 0; n < ninside; ++n)
    {
      // Get the 1D view to the 2D snap->rij array at row n
      auto rij = snap->rij.data_1D(n);

      // Get the pointer to the rij_const array of data
      double const *const rij_const = const_cast<double *>(rij.data());

      snap->compute_duidrj(rij_const, snap->wj[n], snap->rcutij[n], n);

      snap->compute_deidrj(deidrj);

      // Index of the neighbor atom
      int const j = snap->inside[n];

      // Contribution to forces
      if (isComputeForces)
      {
        forces[i][0] += deidrj[0];
        forces[i][1] += deidrj[1];
        forces[i][2] += deidrj[2];

        forces[j][0] -= deidrj[0];
        forces[j][1] -= deidrj[1];
        forces[j][2] -= deidrj[2];
      }

      if (isComputeProcess_dEdr)
      {
        double const rrsq = std::sqrt(rij_const[0] * rij_const[0] + rij_const[1] * rij_const[1] + rij_const[2] * rij_const[2]);
        double dedr = std::sqrt(deidrj[0] * deidrj[0] + deidrj[1] * deidrj[1] + deidrj[2] * deidrj[2]);
        if (!particleContributing[j])
          dedr *= 0.5;
        ier = modelComputeArguments->ProcessDEDrTerm(dedr, rrsq, rij_const, i, j);
        if (ier)
        {
          LOG_ERROR("ProcessDEDrTerm");
          return ier;
        }
      }

      if (isComputeVirial || isComputeParticleVirial)
      {
        // Virial has 6 components and is stored as a 6-element
        // vector in the following order: xx, yy, zz, yz, xz, xy.

        VectorOfSizeSix v;

        v[0] = rij_const[0] * deidrj[0];
        v[1] = rij_const[1] * deidrj[1];
        v[2] = rij_const[2] * deidrj[2];
        v[3] = rij_const[1] * deidrj[2];
        v[4] = rij_const[0] * deidrj[2];
        v[5] = rij_const[0] * deidrj[1];

        if (isComputeVirial)
        {
          virial[0] += v[0];
          virial[1] += v[1];
          virial[2] += v[2];
          virial[3] += v[3];
          virial[4] += v[4];
          virial[5] += v[5];
        }

        if (isComputeParticleVirial)
        {
          v[0] *= 0.5;
          v[1] *= 0.5;
          v[2] *= 0.5;
          v[3] *= 0.5;
          v[4] *= 0.5;
          v[5] *= 0.5;

          particleVirial[i][0] += v[0];
          particleVirial[i][1] += v[1];
          particleVirial[i][2] += v[2];
          particleVirial[i][3] += v[3];
          particleVirial[i][4] += v[4];
          particleVirial[i][5] += v[5];

          particleVirial[j][0] += v[0];
          particleVirial[j][1] += v[1];
          particleVirial[j][2] += v[2];
          particleVirial[j][3] += v[3];
          particleVirial[j][4] += v[4];
          particleVirial[j][5] += v[5];
        } // isComputeParticleVirial
      }   // isComputeVirial || isComputeParticleVirial
    }     // End of loop over neighbors of particle i

    // Energy contribution
    if (isComputeEnergy || isComputeParticleEnergy)
    {
      // Compute contribution to energy.

      // Energy of particle i, sum over coeffs_k * Bi_k

      // Get the 1D view to the 2D coeffelem array at row iSpecies
      auto coeffi = coeffelem.data_1D(iSpecies);

      // Compute phi
      double phi = coeffi[0];

      // Get the pointer to the raw data + 1 to avoid extra summation
      double *Ci = coeffi.data() + 1;

      // Get the bispectrum of particle i
      // Get the 1D view to the 2D bispectrum array at row i
      auto Bi = bispectrum.data_1D(contributing_index);

      // E = beta.B + 0.5*B^t.alpha.B

      // Linear contributions
      for (int icoeff = 0; icoeff < ncoeff; ++icoeff)
        phi += Ci[icoeff] * Bi[icoeff];

      // Quadratic contributions
      if (quadraticflag)
      {
        // Get the pointer to the start of coeffi array of data
        --Ci;

        int k = ncoeff + 1;

        for (int icoeff = 0; icoeff < ncoeff; ++icoeff)
        {
          double const bveci = Bi[icoeff];

          phi += 0.5 * Ci[k++] * bveci * bveci;

          for (int jcoeff = icoeff + 1; jcoeff < ncoeff; ++jcoeff)
          {
            double const bvecj = Bi[jcoeff];

            phi += Ci[k++] * bveci * bvecj;
          }
        }
      } // quadraticflag

      // Contribution to energy
      if (isComputeEnergy)
      {
        *energy += phi;
      }

      // Contribution to particleEnergy
      if (isComputeParticleEnergy)
      {
        particleEnergy[i] += phi;
      }
    } // isComputeEnergy || isComputeParticleEnergy

    ++contributing_index;

  } // End of loop over contributing particles

  // everything is good
  return false;
}
#undef KIM_LOGGER_OBJECT_NAME

#undef ONE
#undef MAXLINE
#undef HELPER_LOG_ERROR
