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
// Copyright (c) 2019--2020, Regents of the University of Minnesota.
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
#define MAXLINE 1024

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
                      chemflag(0),
                      bnormflag(0),
                      wselfallflag(0),
                      beta_max(0),
                      rfac0(0.99363), // Set defaults for optional keywords
                      rmin0(0.0),
                      rcutfac(0.0),
                      snap(nullptr),
                      nAllSpecies(0),
                      nHybridStyleSpecies(0),
                      nzbls(0),
                      // Metal unit as a default base units
                      angstrom(1.0),
                      qqr2e(14.399645354084361),
                      qelectron(1.0),
                      inner(0.0),
                      outer(0.0),
                      zbl(nullptr),
                      ntables(0)
{
  // Everything is good
  *ier = false;

  if (!modelDriverCreate)
  {
    HELPER_LOG_ERROR("The ModelDriverCreate object pointer is not assigned\n");
    *ier = true;
    return;
  }

  {
    std::FILE *parameterFilePointers[MAX_NUM_PARAMETER_FILES];

    int numberParameterFiles(0);

    // Getting the number of parameter files, open the files and process them
    modelDriverCreate->GetNumberOfParameterFiles(&numberParameterFiles);
    if (numberParameterFiles > MAX_NUM_PARAMETER_FILES)
    {
      LOG_ERROR("Too many parameter files!\n");
      *ier = true;
      return;
    }

    if (!numberParameterFiles)
    {
      LOG_ERROR("There is no parameter file!\n");
      *ier = true;
      return;
    }

    *ier = OpenParameterFiles(modelDriverCreate, numberParameterFiles, parameterFilePointers);
    if (*ier)
    {
      return;
    }

    *ier = ProcessParameterFiles(modelDriverCreate, numberParameterFiles, parameterFilePointers);

    CloseParameterFiles(numberParameterFiles, parameterFilePointers);

    if (*ier)
    {
      return;
    }
  }

  *ier = ConvertUnits(modelDriverCreate,
                      requestedLengthUnit,
                      requestedEnergyUnit,
                      requestedChargeUnit,
                      requestedTemperatureUnit,
                      requestedTimeUnit);
  if (*ier)
  {
    return;
  }

  *ier = setRefreshMutableValues(modelDriverCreate);
  if (*ier)
  {
    return;
  }

  *ier = RegisterKIMModelSettings(modelDriverCreate);
  if (*ier)
  {
    return;
  }

  *ier = RegisterKIMParameters(modelDriverCreate);
  if (*ier)
  {
    return;
  }

  *ier = RegisterKIMFunctions(modelDriverCreate);
  if (*ier)
  {
    return;
  }
}
#undef KIM_LOGGER_OBJECT_NAME

SNAPImplementation::~SNAPImplementation()
{
}

int SNAPImplementation::Refresh(KIM::ModelRefresh *const modelRefresh)
{
  return setRefreshMutableValues(modelRefresh);
}

int SNAPImplementation::WriteParameterizedModel(KIM::ModelWriteParameterizedModel const *const modelWriteParameterizedModel) const
{
  if (!modelWriteParameterizedModel)
  {
    HELPER_LOG_ERROR("The modelWriteParameterizedModel object pointer is not assigned\n");
    return true;
  }

  std::string const *path = NULL;
  std::string const *modelName = NULL;

  modelWriteParameterizedModel->GetPath(&path);
  modelWriteParameterizedModel->GetModelName(&modelName);

  std::string coefficient_file(*modelName);
  std::string parameter_file(*modelName);
  std::string hybrid_file(*modelName);

  coefficient_file += ".snapcoeff";
  parameter_file += ".snapparam";
  hybrid_file += ".hybridparam";

  // Set the file name for the next parameter file.
  // It must be called once for each parameter file. The order of these calls is important
  // and determines the order in which the parameter files will be listed in the automatically
  // generated CMakeLists.txt file.
  modelWriteParameterizedModel->SetParameterFileName(coefficient_file);
  modelWriteParameterizedModel->SetParameterFileName(parameter_file);
  if (nzbls || ntables)
  {
    modelWriteParameterizedModel->SetParameterFileName(hybrid_file);
  }

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
       << "# This SNAP coefficient file has " << nelements << ((nelements > 1) ? " blocks, where\n" : " block, where\n")
       << "#    First line of each block contains: 'species name' & 'species radius' & 'species weight'\n"
       << "#    followed by " << ncoeffall << " SNAP coefficients.\n"
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

    // bnormflag
    fs << "bnormflag     " << bnormflag << std::endl;

    // chemflag
    fs << "chemflag      " << chemflag << std::endl;

    // wselfallflag
    fs << "wselfallflag  " << wselfallflag << std::endl;

    fs << "\n\n"
       << "#\n"
       << "# SNAP parameter file"
       << "#\n";

    fs.close();
  }

  // HYBRID parameter file
  if (nzbls || ntables)
  {
    buffer = *path + "/" + hybrid_file;

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

      fs << "# HYBRID style parameters (extra 'zbl' and two-body pair styles\n"
         << "# as a table style can be assigned)\n"
         << "#\n"
         << "# NOTE:\n"
         << "#   No mixing rule will be used here.\n"
         << "#\n"
         << "#   Each pair (I,J) or (J,I) can be assigned to one style\n"
         << "#   If you specify the same pair for the second time, it wipes\n"
         << "#   out all the previous assignments of that pair and the\n"
         << "#   second one will be calculated for the two interacting atoms\n"
         << "#   of those types.\n"
         << "\n\n"
         << "# Number of elements for the hybrid style\n";

      fs << nHybridStyleSpecies;

      fs << "# Number of elements names, (atom names)\n";
      for (int iSpec = 0; iSpec < nHybridStyleSpecies; ++iSpec)
      {
        fs << hybridStyleSpeciesNames[iSpec] << "  ";
      }
      fs << "\n"
         << std::endl;

      fs << std::setprecision(std::numeric_limits<double>::digits10 + 1) << std::fixed;

      if (nzbls)
      {
        fs << "# In the ZBL style, the inner and outer cutoff are the\n"
           << "# same for all pairs of atom types.\n"
           << "# zbl  inner  outer\n";
        fs << "zbl  " << inner << "  " << outer << std::endl;
        fs << "# Element_1  Element_2  zbl  Z_1  Z_2\n";
        for (int iSpecies = 0; iSpecies < nAllSpecies; ++iSpecies)
        {
          for (int jSpecies = iSpecies; jSpecies < nAllSpecies; ++jSpecies)
          {
            if (setflag(iSpecies, jSpecies) == HYBRIDSTYLE::ZBL)
            {
              fs << allSpeciesNames[iSpecies] << "  " << allSpeciesNames[jSpecies] << "  zbl  " << atomicNumber[iSpecies] << "  " << atomicNumber[jSpecies] << std::endl;
            }
          }
        }
      }

      if (ntables)
      {
        fs << "\n"
           << "# Style table:\n"
           << "# table  style  N\n";
        for (auto i : tables_info)
        {
          if (i.tableStyle == TABLESTYLE::LOOKUP)
          {
            fs << "table  lookup  " << i.tableLength << std::endl;
          }
          else if (i.tableStyle == TABLESTYLE::LINEAR)
          {
            fs << "table  linear  " << i.tableLength << std::endl;
          }
          else if (i.tableStyle == TABLESTYLE::SPLINE)
          {
            fs << "table  spline  " << i.tableLength << std::endl;
          }
          else if (i.tableStyle == TABLESTYLE::BITMAP)
          {
            fs << "table  bitmap  " << i.tableLength << std::endl;
          }
        }
        fs << "#\n"
           << "# style = 'lookup' or 'linear' or 'spline' or 'bitmap' =\n"
           << "# method of interpolation\n"
           << "# N     = use N values in 'lookup', 'linear', 'spline' tables\n"
           << "# N     = use 2^N values in 'bitmap' tables\n"
           << "#\n"
           << "#\n"
           << "# 'table' style can be used multiple times. For example, interactions between\n"
           << "# I and I atoms use a 'linear' table style and interactions between I and J\n"
           << "# atoms use a 'spline' table style, then you should list the table style two\n"
           << "# times. Style indexing starts from '1'. It means, that the fisrt style is\n"
           << "# numbered '1' and the second '2' and so on so forth.\n"
           << "# Later in the pair interactions, the 'table' style must be added after the\n"
           << "# I,J atom names followed by the style number and then followed by the\n"
           << "# remaining coefficients as of 'filename', and 'keyword', and 'cutoff'.\n"
           << "#\n"
           << "#\n"
           << "# The 'filename' specifies a file containing tabulated energy\n"
           << "# and force values.\n"
           << "# The 'keyword' specifies a section of the file.\n"
           << "# The 'cutoff' is an optional coefficient in distance unit.\n"
           << "#\n"
           << "#\n"
           << "# Example 1, we use 2 styles\n"
           << "# table  linear  1000\n"
           << "# table  spline  10000\n"
           << "# I  I  table  1  II.table IIKey 4.0\n"
           << "# I  J  table  2  IJ.table IJ\n"
           << "#\n"
           << "#\n"
           << "# Example 2, we use one style\n"
           << "# table  spline  10000\n"
           << "# I  I  table  1  II.table IIKey 4.0\n"
           << "# I  J  table  1  IJ.table IJ\n"
           << "#\n"
           << "#\n"
           << "# Example 3, we use one style and have one tabulated file\n"
           << "# table  spline  10000\n"
           << "# I  I  table  1  tablefile.txt II 4.0\n"
           << "# I  J  table  1  tablefile.txt IJ 4.8\n"
           << "#\n"
           << "#\n"
           << "#- Element_1  Element_2  table  style_number  filename  keyword  [cutoff]\n";
        for (int iSpecies = 0; iSpecies < nAllSpecies; ++iSpecies)
        {
          for (int jSpecies = iSpecies; jSpecies < nAllSpecies; ++jSpecies)
          {
            if (setflag(iSpecies, jSpecies) == HYBRIDSTYLE::TABLE)
            {
              fs << allSpeciesNames[iSpecies] << "  " << allSpeciesNames[jSpecies] << "  table  " << tableNumber(iSpecies, jSpecies) << "  filename  keyword  [cutoff]" << std::endl;
            }
          }
        }
      }
      fs << "\n\n"
         << "#\n"
         << "# SNAP HYBRID style parameter file"
         << "#\n";
      fs.close();
    }
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
  bool isHybrid = (nzbls || ntables);

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
    return true;
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
  return RegisterKIMComputeArgumentsSettings(modelComputeArgumentsCreate);
}

int SNAPImplementation::ComputeArgumentsDestroy(KIM::ModelComputeArgumentsDestroy *const /*modelComputeArgumentsDestroy*/) const
{
  // Nothing to do for this case
  // Everything is good
  return false;
}

// Private member functions

#define KIM_LOGGER_OBJECT_NAME modelDriverCreate
int SNAPImplementation::OpenParameterFiles(KIM::ModelDriverCreate *const modelDriverCreate,
                                           int const numberParameterFiles,
                                           std::FILE **parameterFilePointers)
{
  std::string const *parameterFileDirectoryName;
  modelDriverCreate->GetParameterFileDirectoryName(&parameterFileDirectoryName);

  for (int i = 0; i < numberParameterFiles; ++i)
  {
    std::string const *parameterFileBasename;

    if (modelDriverCreate->GetParameterFileBasename(i, &parameterFileBasename))
    {
      LOG_ERROR("Unable to get the parameter file base name\n");
      return true;
    }

    std::string const parameterFileName = *parameterFileDirectoryName + "/" + *parameterFileBasename;

    parameterFilePointers[i] = std::fopen(parameterFileName.c_str(), "r");
    if (!parameterFilePointers[i])
    {
      HELPER_LOG_ERROR("The parameter file (" + *parameterFileBasename + ") can not be opened\n");
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

#define KIM_LOGGER_OBJECT_NAME modelDriverCreate
int SNAPImplementation::ProcessParameterFiles(KIM::ModelDriverCreate *const modelDriverCreate,
                                              int const numberParameterFiles,
                                              std::FILE *const *parameterFilePointers)
{
  char nextLine[MAXLINE];

  int endOfFileFlag;

  // keep track of known species
  std::map<KIM::SpeciesName const, int, KIM::SPECIES_NAME::Comparator> speciesMap;

  // Read SNAP coefficient file
  GetNextDataLine(parameterFilePointers[0], nextLine, MAXLINE, &endOfFileFlag);
  if (endOfFileFlag)
  {
    HELPER_LOG_ERROR("End of file in the SNAP coefficient file.\n");
    return true;
  }

  int ier = std::sscanf(nextLine, "%d %d", &nelements, &ncoeffall);
  if (ier != 2)
  {
    HELPER_LOG_ERROR("Unable to read nelements & ncoeffall from the line:\n" +
                     std::string(nextLine) +
                     "\n");
    return true;
  }

  if (nelements < 1)
  {
    HELPER_LOG_ERROR("Incorrect number of elements in the SNAP coefficient file.\n");
    return true;
  }

  if (ncoeffall < 1)
  {
    HELPER_LOG_ERROR("Incorrect number of coefficients in the SNAP coefficient file.\n");
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
      HELPER_LOG_ERROR("End of file in the SNAP coefficient file.\n");
      return true;
    }

    char spec[MAXLINE];
    ier = std::sscanf(nextLine, "%s %lf %lf", spec, &radelem[ielem], &wjelem[ielem]);
    if (ier != 3)
    {
      HELPER_LOG_ERROR("Incorrect format in the SNAP coefficient file.\n"
                       "Unable to read spec, radius & weight from the line:\n" +
                       std::string(nextLine) +
                       "\n");
      return true;
    }

    if (radelem[ielem] < 0.0)
    {
      HELPER_LOG_ERROR("Incorrect radius in the SNAP coefficient file.\n");
      return true;
    }

    KIM::SpeciesName const speciesName(spec);

    if (speciesName.ToString() == "unknown")
    {
      HELPER_LOG_ERROR("Incorrect format in the SNAP coefficient file.\n"
                       "The species name '" +
                       std::string(spec) +
                       "' is unknown.\n");
      return true;
    }

    // check for new species
    auto iter = speciesMap.find(speciesName);
    if (iter == speciesMap.end())
    {
      ier = modelDriverCreate->SetSpeciesCode(speciesName, ielem);
      if (ier)
      {
        LOG_ERROR("SetSpeciesCode failed to set the new species.\n");
        return ier;
      }

      speciesMap[speciesName] = ielem;

      elements[ielem] = speciesName.ToString();
    }
    else
    {
      HELPER_LOG_ERROR("Incorrect format in the SNAP coefficient file.\n"
                       "The Species '" +
                       std::string(spec) +
                       "' is already defined and exist.\n");
      return true;
    }

    for (int icoeff = 0; icoeff < ncoeffall; ++icoeff)
    {
      GetNextDataLine(parameterFilePointers[0], nextLine, MAXLINE, &endOfFileFlag);
      if (endOfFileFlag)
      {
        HELPER_LOG_ERROR("End of file in the SNAP coefficient file.\n");
        return true;
      }

      double coeff;
      ier = std::sscanf(nextLine, "%lf", &coeff);
      if (ier != 1)
      {
        HELPER_LOG_ERROR("Incorrect format in the SNAP coefficient file.\n");
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
    {
      break;
    }

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
    else if (!std::strcmp(keywd, "chemflag"))
    {
      chemflag = std::atoi(keyval);
      if (chemflag != 0 && chemflag != 1)
        chemflag = 1;
    }
    else if (!std::strcmp(keywd, "bnormflag"))
    {
      bnormflag = std::atoi(keyval);
      if (bnormflag != 0 && bnormflag != 1)
        bnormflag = 1;
    }
    else if (!std::strcmp(keywd, "wselfallflag"))
    {
      wselfallflag = std::atoi(keyval);
      if (wselfallflag != 0 && wselfallflag != 1)
        wselfallflag = 1;
    }
    else if (!std::strcmp(keywd, "chunksize"))
    {
      // The keyword chunksize is ignored.
    }
    else if (!std::strcmp(keywd, "diagonalstyle"))
    {
      int const diagonalstyle = std::atoi(keyval);
      if (diagonalstyle != 3)
      {
        HELPER_LOG_ERROR("Incorrect SNAP parameter file.\n"
                         "The 'diagonalstyle' keyword was removed in 2019 "
                         "since all known SNAP potentials use the default "
                         "value of '3'. (This keyword is kept only for "
                         "backward compatibility.)\n"
                         "Incorrect diagonalstyle value of " +
                         std::string(keyval) +
                         " is given\n");
        return true;
      }
    }
    else
    {
      HELPER_LOG_ERROR("Incorrect SNAP parameter file.\n"
                       "Incorrect SNAP keywords of " +
                       std::string(keywd) +
                       ", " +
                       std::string(keyval) +
                       "\n");
      return true;
    }
  }

  if (!rcutfacflag || !twojmaxflag)
  {
    HELPER_LOG_ERROR("Incorrect SNAP parameter file.\n" +
                     std::string(rcutfacflag ? "'twojmaxflag' is not given\n" : "'rcutfacflag' is not given\n"));
    return true;
  }

  // Process Hybrid parameter file
  if (numberParameterFiles > 2)
  {
    // Read the Hybrid parameter file
    GetNextDataLine(parameterFilePointers[2], nextLine, MAXLINE, &endOfFileFlag);
    if (endOfFileFlag)
    {
      HELPER_LOG_ERROR("End of file in the Hybrid parameter file.\n");
      return true;
    }

    int ier = std::sscanf(nextLine, "%d", &nHybridStyleSpecies);
    if (ier != 1)
    {
      HELPER_LOG_ERROR("Unable to read the number of unique species from the line:\n" +
                       std::string(nextLine) +
                       "\n");
      return true;
    }

    if (nHybridStyleSpecies <= 0)
    {
      HELPER_LOG_ERROR("Incorrect number of unique species '" +
                       std::to_string(nHybridStyleSpecies) +
                       "' in the Hybrid parameter file.\n");
      return true;
    }

    // Allocate memory based on the number of species

    nAllSpecies = nelements;

    allSpeciesNames.reserve(nelements + nHybridStyleSpecies);
    for (auto i : elements)
    {
      allSpeciesNames.push_back(i);
    }

    // Set up memory for the species lists
    hybridStyleSpeciesNames.resize(nHybridStyleSpecies);

    GetNextDataLine(parameterFilePointers[2], nextLine, MAXLINE, &endOfFileFlag);
    if (endOfFileFlag)
    {
      HELPER_LOG_ERROR("End of file in the Hybrid parameter file.\n");
      return true;
    }

    for (int iSpec = 0; iSpec < nHybridStyleSpecies; ++iSpec)
    {
      char *spec = std::strtok(iSpec ? NULL : nextLine, "' \t\n\r\f");
      if (!spec)
      {
        HELPER_LOG_ERROR("Unable to read the name of unique species from the line:\n" +
                         std::string(nextLine) +
                         "\n");
        return true;
      }

      KIM::SpeciesName const speciesName(spec);

      if (speciesName.ToString() == "unknown")
      {
        HELPER_LOG_ERROR("Incorrect format in the Hybrid parameter file.\n"
                         "The species name '" +
                         std::string(spec) +
                         "' is unknown.\n");
        return true;
      }

      for (int i = 0; i < iSpec; ++i)
      {
        KIM::SpeciesName const tmpName(hybridStyleSpeciesNames[i]);
        if (speciesName == tmpName)
        {
          HELPER_LOG_ERROR("Incorrect format in the Hybrid parameter file.\n"
                           "The Species '" +
                           std::string(spec) +
                           "' is already defined and exist.\n");
          return true;
        }
      }

      hybridStyleSpeciesNames[iSpec] = speciesName.ToString();

      // check for new species
      auto iter = speciesMap.find(speciesName);
      if (iter == speciesMap.end())
      {
        ier = modelDriverCreate->SetSpeciesCode(speciesName, nAllSpecies);
        if (ier)
        {
          LOG_ERROR("SetSpeciesCode failed to set the new species.\n");
          return ier;
        }

        speciesMap[speciesName] = nAllSpecies++;

        allSpeciesNames.push_back(speciesName.ToString());
      }
    } // End of loop for the number of species

    while (1)
    {
      GetNextDataLine(parameterFilePointers[2], nextLine, MAXLINE, &endOfFileFlag);
      if (endOfFileFlag)
      {
        break;
      }

      char *keywd = std::strtok(nextLine, "' \t\n\r\f");
      char *keyval1 = std::strtok(NULL, "' \t\n\r\f");
      char *keyval2 = std::strtok(NULL, "' \t\n\r\f");

      if (!std::strcmp(keywd, "zbl"))
      {
        ++nzbls;
        if (nzbls >= 2)
        {
          HELPER_LOG_ERROR("There is more than one zbl style in the Hybrid parameter file\n");
          return true;
        }

        if (!keyval1 || !keyval2)
        {
          HELPER_LOG_ERROR("Incorrect zbl style in the Hybrid parameter file.\n");
          return true;
        }

        inner = std::atof(keyval1);
        outer = std::atof(keyval2);

        if (inner <= 0.0)
        {
          HELPER_LOG_ERROR("Incorrect distance (where switching function begins) for ZBL interaction in the Hybrid parameter file.\n");
          return true;
        }

        if (inner > outer)
        {
          HELPER_LOG_ERROR("Incorrect global cutoff for ZBL interaction in the Hybrid parameter file.\n");
          return true;
        }

        continue;
      } // if zbl

      if (keyval2)
      {
        if (!std::strcmp(keyval2, "table"))
        {
          ++ntables;
          continue;
        }
      }
    } // End of while loop

    if (!nzbls && !ntables)
    {
      HELPER_LOG_ERROR("There is no zbl/table style in the Hybrid parameter file.\n");
      return true;
    }

    // Allocate memory
    setflag.resize(nAllSpecies, nAllSpecies, HYBRIDSTYLE::NONE);

    if (nzbls)
    {
      atomicNumber.resize(nAllSpecies, 0.0);
    }

    if (ntables)
    {
      tables_info.reserve(ntables);

      tables.resize(ntables);

      tableNumber.resize(nAllSpecies, nAllSpecies, -1);
    }

    // Rewind the file
    if (std::fseek(parameterFilePointers[2], 0, SEEK_SET))
    {
      HELPER_LOG_ERROR("Failed to rewind the Hybrid parameter file.\n");
      return true;
    }

    int nTableStyles = 0;
    int nTables = 0;

    while (1)
    {
      GetNextDataLine(parameterFilePointers[2], nextLine, MAXLINE, &endOfFileFlag);
      if (endOfFileFlag)
      {
        break;
      }

      char *keywd = std::strtok(nextLine, "' \t\n\r\f");
      char *keyval1 = std::strtok(NULL, "' \t\n\r\f");
      char *keyval2 = std::strtok(NULL, "' \t\n\r\f");

      if (!keyval1 || !keyval2)
      {
        continue;
      }

      if (!std::strcmp(keywd, "zbl"))
      {
        continue;
      }

      if (!std::strcmp(keywd, "table"))
      {
        if (!keyval1 || !keyval2)
        {
          HELPER_LOG_ERROR("Wrong table style in the Hybrid parameter file.\n");
          return true;
        }

        TABLESTYLE tabstyle;

        if (!std::strcmp(keyval1, "lookup"))
        {
          tabstyle = TABLESTYLE::LOOKUP;
        }
        else if (!std::strcmp(keyval1, "linear"))
        {
          tabstyle = TABLESTYLE::LINEAR;
        }
        else if (!std::strcmp(keyval1, "spline"))
        {
          tabstyle = TABLESTYLE::SPLINE;
        }
        else if (!std::strcmp(keyval1, "bitmap"))
        {
          tabstyle = TABLESTYLE::BITMAP;
        }
        else
        {
          HELPER_LOG_ERROR("Unknown table style '" +
                           std::string(keyval1) +
                           "' in the Hybrid parameter file.\n");
          return true;
        }

        int const tablength = std::atoi(keyval2);
        if (tablength < 2)
        {
          HELPER_LOG_ERROR("Illegal number of table entries " +
                           std::string(keyval2) +
                           " < 2\n");
          return true;
        }

        tables_info.resize(nTableStyles + 1);

        tables_info[nTableStyles].tableStyle = tabstyle;
        tables_info[nTableStyles].tableLength = tablength;
        ++nTableStyles;

        continue;
      }

      if ((!std::strcmp(keyval2, "zbl")) || (!std::strcmp(keyval2, "table")))
      {
        char *keyval3 = std::strtok(NULL, "' \t\n\r\f");
        char *keyval4 = std::strtok(NULL, "' \t\n\r\f");
        char *keyval5 = std::strtok(NULL, "' \t\n\r\f");
        char *keyval6 = std::strtok(NULL, "' \t\n\r\f");

        if (!keyval3 || !keyval4)
        {
          HELPER_LOG_ERROR("Incorrect zbl/table style in the Hybrid parameter file.\n");
          return true;
        }

        // convert species strings to proper type instances
        KIM::SpeciesName const specName1(keywd);
        KIM::SpeciesName const specName2(keyval1);

        // Extra check for the species
        auto iter = speciesMap.find(specName1);
        if (iter == speciesMap.end())
        {
          HELPER_LOG_ERROR("The species name '" +
                           specName1.ToString() +
                           "' is not defined\n");
          return true;
        }
        int const iSpecies = speciesMap[specName1];

        // Extra check for the species
        iter = speciesMap.find(specName2);
        if (iter == speciesMap.end())
        {
          HELPER_LOG_ERROR("The species name '" +
                           specName2.ToString() +
                           "' is not defined\n");
          return true;
        }
        int const jSpecies = speciesMap[specName2];

        if (!std::strcmp(keyval2, "zbl"))
        {
          if (keyval5 || keyval6)
          {
            HELPER_LOG_ERROR("Incorrect zbl style in the Hybrid parameter file.\n");
            return true;
          }

          // Get the species atomic numbers
          double const z_iSpecies = std::atof(keyval3);
          double const z_jSpecies = std::atof(keyval4);

          if (iSpecies == jSpecies)
          {
            if (z_iSpecies != z_jSpecies)
            {
              HELPER_LOG_ERROR("Incorrect ZBL args in the Hybrid parameter file.\n"
                               "When i == j, it is required that Z_i == Z_j.\n");
              return true;
            }
          }
          else
          {
            setflag(jSpecies, iSpecies) = HYBRIDSTYLE::ZBL;

            atomicNumber[jSpecies] = z_jSpecies;
          }

          setflag(iSpecies, jSpecies) = HYBRIDSTYLE::ZBL;

          atomicNumber[iSpecies] = z_iSpecies;

          continue;
        } // if zbl

        if (!std::strcmp(keyval2, "table"))
        {
          // specName1  specName2  table  1  filename  keyword  [cutoff]

          // Get the table style number
          int const tableStyleNumber = std::stoi(keyval3);
          if (tableStyleNumber > nTableStyles)
          {
            HELPER_LOG_ERROR("Incorrect style number in the Hybrid parameter file.\n"
                             " style number > number of defined table styles.\n");
            return true;
          }

          int const sn = tableStyleNumber - 1;
          int const nt = nTables++;

          // Construct a table object
          TABLE tt(tables_info[sn].tableStyle, tables_info[sn].tableLength);
          tables[nt] = std::move(tt);

          int tableFileNumber = 3;
          for (; tableFileNumber < numberParameterFiles; ++tableFileNumber)
          {
            std::string const *parameterFileBasename;

            if (modelDriverCreate->GetParameterFileBasename(tableFileNumber, &parameterFileBasename))
            {
              LOG_ERROR("Unable to get the parameter file base name\n");
              return true;
            }

            if (!std::strcmp(keyval4, (*parameterFileBasename).c_str()))
            {
              break;
            }
          } // End of loop throug the number of registered parameter files

          if (tableFileNumber >= numberParameterFiles)
          {
            HELPER_LOG_ERROR("Incorrect TABLE filename '" +
                             std::string(keyval4) +
                             "' in the Hybrid parameter file.\n");
            return true;
          }

          // Table keyword
          if (!keyval5)
          {
            HELPER_LOG_ERROR("Incorrect table style in the Hybrid parameter file (No table keyword).\n");
            return true;
          }

          // Read the table
          if (tables[nt].read_table(parameterFilePointers[tableFileNumber], keyval5))
          {
            HELPER_LOG_ERROR("Failed to read from the table file '" +
                             std::string(keyval4) +
                             "'.\n");
            return true;
          }

          // Check the table parameters and insure cutoff is within table
          // for BITMAP tables, file values can be in non-ascending order
          if (tables[nt].ninput <= 1)
          {
            HELPER_LOG_ERROR("Invalid table length '" +
                             std::to_string(tables[nt].ninput) +
                             "' <= 1 read from the table file '" +
                             std::string(keyval4) +
                             "'.\n");
            return true;
          }

          // Set the table cutoff
          if (keyval6)
          {
            tables[nt].cut = std::atof(keyval6);

            if (tables[nt].cut < 0.0)
            {
              HELPER_LOG_ERROR("Invalid table cutoff '" +
                               std::to_string(tables[nt].cut) +
                               "' < 0.0 in the Hybrid parameter file.\n");
              return true;
            }
          }
          else if (tables[nt].rflag != TABLEDISTANCESTYLE::NONE)
          {
            tables[nt].cut = tables[nt].rhi;
          }
          else
          {
            tables[nt].cut = tables[nt].rfile[tables[nt].ninput - 1];
          }

          double rlo;
          double rhi;

          if (tables[nt].rflag == TABLEDISTANCESTYLE::NONE)
          {
            rlo = tables[nt].rfile[0];
            rhi = tables[nt].rfile[tables[nt].ninput - 1];
          }
          else
          {
            rlo = tables[nt].rlo;
            rhi = tables[nt].rhi;
          }

          if (tables[nt].cut <= rlo || tables[nt].cut > rhi)
          {
            HELPER_LOG_ERROR("Invalid table cutoff.\n");
            return true;
          }

          if (rlo <= 0.0)
          {
            HELPER_LOG_ERROR("Invalid table rlo <= 0.0.\n");
            return true;
          }

          // match = 1, if don't need to spline read-in tables this is only
          // the case if r values needed by final tables exactly match r
          // values read from the file for tabstyle SPLINE, always need to
          // build spline tables
          tables[nt].match = 0;

          if (tables[nt].tableStyle == TABLESTYLE::LINEAR &&
              tables[nt].ninput == tables[nt].tableLength &&
              tables[nt].rflag == TABLEDISTANCESTYLE::RSQ &&
              tables[nt].rhi == tables[nt].cut)
          {
            tables[nt].match = 1;
          }

          if (tables[nt].tableStyle == TABLESTYLE::BITMAP &&
              tables[nt].ninput == 1 << tables[nt].tableLength &&
              tables[nt].rflag == TABLEDISTANCESTYLE::BITMAP &&
              tables[nt].rhi == tables[nt].cut)
          {
            tables[nt].match = 1;
          }

          if (tables[nt].rflag == TABLEDISTANCESTYLE::BITMAP &&
              tables[nt].match == 0)
          {
            HELPER_LOG_ERROR("Bitmapped table in the file does not match the requested table.\n");
            return true;
          }

          if (iSpecies != jSpecies)
          {
            setflag(jSpecies, iSpecies) = HYBRIDSTYLE::TABLE;

            tableNumber(jSpecies, iSpecies) = nt;
          }

          setflag(iSpecies, jSpecies) = HYBRIDSTYLE::TABLE;

          tableNumber(iSpecies, jSpecies) = nt;

          continue;
        } // if table
      }   // if zbl or table
    }     // End of while loop
  }       // if numberParameterFiles > 2

  // Extra check to set the flag for all species defined to have interactions
  // in the SNAP potential
  snapflag.resize(nAllSpecies ? nAllSpecies : nelements, false);
  for (auto elem : elements)
  {
    int const iSpecies = speciesMap[elem];
    snapflag[iSpecies] = true;
  }

  // everything is good
  return false;
}
#undef KIM_LOGGER_OBJECT_NAME

void SNAPImplementation::CloseParameterFiles(int const numberParameterFiles,
                                             std::FILE *const *parameterFilePointers)
{
  for (int i = 0; i < numberParameterFiles; ++i)
  {
    std::fclose(parameterFilePointers[i]);
  }
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

  // changing units of length
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
                                           ONE,
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

  // Convert to active units
  if (convertLength != ONE)
  {
    std::for_each(radelem.begin(), radelem.end(), [&](double &rc) { rc *= convertLength; });
    rmin0 *= convertLength;
    if (nzbls)
    {
      angstrom *= convertLength;
      qqr2e *= convertLength;
      inner *= convertLength;
      outer *= convertLength;
    }
    if (ntables)
    {
      for (int nt = 0; nt < ntables; ++nt)
      {
        tables[nt].rlo *= convertLength;
        tables[nt].rhi *= convertLength;
        tables[nt].fplo /= convertLength;
        tables[nt].fplo /= convertLength;
        tables[nt].fphi /= convertLength;
        tables[nt].fphi /= convertLength;
        std::for_each(tables[nt].rfile.begin(), tables[nt].rfile.end(), [&](double &r) { r *= convertLength; });
        std::for_each(tables[nt].ffile.begin(), tables[nt].ffile.end(), [&](double &f) { f /= convertLength; });
      }
    }
  }

  // Changing units of energy
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
                                       ONE,
                                       0.0,
                                       0.0,
                                       0.0,
                                       &convertEnergy);
  if (ier)
  {
    LOG_ERROR("Unable to convert energy unit");
    return ier;
  }

  // Convert to active units
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
    if (nzbls)
    {
      qqr2e *= convertEnergy;
    }
    if (ntables)
    {
      for (int nt = 0; nt < ntables; ++nt)
      {
        tables[nt].fplo *= convertEnergy;
        tables[nt].fphi *= convertEnergy;
        std::for_each(tables[nt].ffile.begin(), tables[nt].ffile.end(), [&](double &f) { f *= convertEnergy; });
        std::for_each(tables[nt].efile.begin(), tables[nt].efile.end(), [&](double &e) { e *= convertEnergy; });
      }
    }
  }

  // Changing units of charge
  double convertCharge = ONE;
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
                                       0.0,
                                       ONE,
                                       0.0,
                                       0.0,
                                       &convertCharge);
  if (ier)
  {
    LOG_ERROR("Unable to convert energy unit");
    return ier;
  }

  // Convert to active units
  if (convertCharge != ONE)
  {
    if (nzbls)
    {
      qqr2e /= convertCharge;
      qqr2e /= convertCharge;
      qelectron *= convertCharge;
    }
  }

  // Register units
  ier = modelDriverCreate->SetUnits(requestedLengthUnit,
                                    requestedEnergyUnit,
                                    requestedChargeUnit,
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
      HELPER_LOG_ERROR("Incorrect args.\n"
                       "ncoeffall = " +
                       std::to_string(ncoeffall) +
                       " ntmp = " +
                       std::to_string(ntmp) +
                       " ncoeff = " +
                       std::to_string(ncoeff) +
                       "\n");
      return true;
    }
  }
  else
  {
    ncoeff = ncoeffall - 1;
  }

  if (bnormflag != chemflag)
  {
    HELPER_LOG_WARNING("Incorrect SNAP parameter file.\n"
                       "'bnormflag' and 'chemflag' are not equal. "
                       "This is probably not what you intended\n")
  }

  // We have to creat the SNAP object and do the extra check on
  // the number of coefficients

  // Construct the SNAP object
  snap.reset(new SNA(rfac0, twojmax, rmin0,
                     switchflag, bzeroflag,
                     chemflag, bnormflag,
                     wselfallflag, nelements));

  // Extra check
  if (ncoeff != snap->ncoeff)
  {
    HELPER_LOG_ERROR("Wrong number of coefficients.\n"
                     "ncoeff = " +
                     std::to_string(ncoeff) +
                     " & SNAP ncoeff = " +
                     std::to_string(snap->ncoeff) +
                     "\nNumber of coefficients to the model and the one created in SNAP object do not match\n");
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

  if (nzbls)
  {
    // Construct the ZBL object
    zbl.reset(new ZBL(inner, outer));

    zbl->allocate(nAllSpecies);

    rcutmax = std::max(rcutmax, zbl->cut_global);

    for (int iSpecies = 0; iSpecies < nAllSpecies; ++iSpecies)
    {
      for (int jSpecies = iSpecies; jSpecies < nAllSpecies; ++jSpecies)
      {
        if (setflag(iSpecies, jSpecies) == HYBRIDSTYLE::ZBL)
        {
          zbl->set_coeff(iSpecies, jSpecies, atomicNumber[iSpecies], atomicNumber[jSpecies], angstrom, qqr2e, qelectron);
        }
      }
    }
  }

  if (ntables)
  {
    for (int nt = 0; nt < ntables; ++nt)
    {
      rcutmax = std::max(rcutmax, tables[nt].cut);

      // spline read-in values and compute r,e,f vectors within table
      if (tables[nt].match == 0)
      {
        tables[nt].spline_table();
      }

      tables[nt].compute_table();
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
                                        bool const isComputeParticleVirial,
                                        bool const isHybrid) const
{
  int const processd2E = 2;
  int const energy = 2;
  int const force = 2;
  int const particleEnergy = 2;
  int const virial = 2;
  int const particleVirial = 2;
  int const hybrid = 2;

  int const isComputeProcess_dEdr_int = static_cast<int>(isComputeProcess_dEdr);
  int const isComputeProcess_d2Edr2_int = static_cast<int>(isComputeProcess_d2Edr2);
  int const isComputeEnergy_int = static_cast<int>(isComputeEnergy);
  int const isComputeForces_int = static_cast<int>(isComputeForces);
  int const isComputeParticleEnergy_int = static_cast<int>(isComputeParticleEnergy);
  int const isComputeVirial_int = static_cast<int>(isComputeVirial);
  int const isComputeParticleVirial_int = static_cast<int>(isComputeParticleVirial);
  int const isHybrid_int = static_cast<int>(isHybrid);

  int index = 0;
  // processdE
  index += isComputeProcess_dEdr_int * processd2E * energy * force * particleEnergy * virial * particleVirial * hybrid;
  // processd2E
  index += isComputeProcess_d2Edr2_int * energy * force * particleEnergy * virial * particleVirial * hybrid;
  // energy
  index += isComputeEnergy_int * force * particleEnergy * virial * particleVirial * hybrid;
  // force
  index += isComputeForces_int * particleEnergy * virial * particleVirial * hybrid;
  // particleEnergy
  index += isComputeParticleEnergy_int * virial * particleVirial * hybrid;
  // virial
  index += isComputeVirial_int * particleVirial * hybrid;
  // particleVirial
  index += isComputeParticleVirial_int * hybrid;
  // hybrid
  index += isHybrid_int;
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
                                               "0 or 1\n "
                                               "0 = do not use switching function,\n "
                                               "1 = use switching function, 'switchflag'.");
  if (ier)
  {
    LOG_ERROR("SetParameterPointer switchflag");
    return ier;
  }

  ier = modelDriverCreate->SetParameterPointer(1, &bzeroflag, "bzeroflag",
                                               "0 or 1\n "
                                               "0 = do not subtract B0,\n "
                                               "1 = subtract B0, 'bzeroflag'.");
  if (ier)
  {
    LOG_ERROR("SetParameterPointer bzeroflag");
    return ier;
  }

  ier = modelDriverCreate->SetParameterPointer(1, &quadraticflag, "quadraticflag",
                                               "0 or 1\n "
                                               "0 = do not generate quadratic terms,\n "
                                               "1 = generate quadratic terms, 'quadraticflag'.");
  if (ier)
  {
    LOG_ERROR("SetParameterPointer quadraticflag");
    return ier;
  }

  ier = modelDriverCreate->SetParameterPointer(1, &bnormflag, "bnormflag",
                                               "0 or 1\n "
                                               "0 = do nothing,\n "
                                               "1 = if barray divided by j+1.");
  if (ier)
  {
    LOG_ERROR("SetParameterPointer bnormflag");
    return ier;
  }

  ier = modelDriverCreate->SetParameterPointer(1, &chemflag, "chemflag",
                                               "0 or 1\n "
                                               "0 = do nothing,\n "
                                               "1 = for multi-element bispectrum components.");
  if (ier)
  {
    LOG_ERROR("SetParameterPointer chemflag");
    return ier;
  }

  ier = modelDriverCreate->SetParameterPointer(1, &wselfallflag, "wselfallflag",
                                               "0 or 1\n "
                                               "0 = do nothing,\n "
                                               "1 = add wself to all element labelings.");
  if (ier)
  {
    LOG_ERROR("SetParameterPointer wselfallflag");
    return ier;
  }

  if (nzbls)
  {
    ier = modelDriverCreate->SetParameterPointer(1, &inner, "zbl_inner",
                                                 "Distance where switching function in the ZBL interaction begins.");
    if (ier)
    {
      LOG_ERROR("SetParameterPointer switchflag");
      return ier;
    }

    ier = modelDriverCreate->SetParameterPointer(1, &outer, "zbl_outer",
                                                 "Distance Global cutoff for the ZBL interaction.");
    if (ier)
    {
      LOG_ERROR("SetParameterPointer switchflag");
      return ier;
    }
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
    {
      continue;
    }

    // Get the species index for atom i
    int const iSpecies = particleSpeciesCodes[i];

    if (!snapflag[iSpecies])
    {
      continue;
    }

    {
      // Get neighbors of i
      modelComputeArguments->GetNeighborList(0, i, &numnei, &n1atom);

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

        if (snapflag[jSpecies])
        {
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
            snap->element[ninside] = jSpecies;
            ++ninside;
          }
        }
      }

      if (chemflag)
      {
        snap->compute_ui(ninside, iSpecies);
      }
      else
      {
        snap->compute_ui(ninside, 0);
      }

      snap->compute_zi();

      if (chemflag)
      {
        snap->compute_bi(iSpecies);
      }
      else
      {
        snap->compute_bi(0);
      }

      auto Bi = bispectrum.data_1D(contributing_index++);
      for (int icoeff = 0; icoeff < ncoeff; ++icoeff)
      {
        Bi[icoeff] = snap->blist[icoeff];
      }
    }
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
      {
        continue;
      }

      // Get the species index for atom i
      int const iSpecies = particleSpeciesCodes[i];

      if (!snapflag[iSpecies])
      {
        continue;
      }

      {
        // Get the 1D view to the 2D coeffelem array at row iSpecies
        auto coeffi = coeffelem.data_1D(iSpecies);

        // Get the pointer to the raw data + 1 to avoid extra sum
        double *Ci = coeffi.data() + 1;

        // Get the 1D view to the 2D beta array at row i
        auto bi = beta.data_1D(contributing_index);

        for (int icoeff = 0; icoeff < ncoeff; ++icoeff)
        {
          bi[icoeff] = Ci[icoeff];
        }

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
  }
  else
  {
    for (int i = 0, contributing_index = 0; i < cachedNumberOfParticles_; ++i)
    {
      if (!particleContributing[i])
      {
        continue;
      }

      // Get the species index for atom i
      int const iSpecies = particleSpeciesCodes[i];

      if (snapflag[iSpecies])
      {
        // Get the 1D view to the 2D coeffelem array at row iSpecies
        auto coeffi = coeffelem.data_1D(iSpecies);

        // Get the pointer to the raw data + 1 to avoid extra sum
        double *Ci = coeffi.data() + 1;

        // Get the 1D view to the 2D beta array at row i
        auto bi = beta.data_1D(contributing_index++);

        for (int icoeff = 0; icoeff < ncoeff; ++icoeff)
        {
          bi[icoeff] = Ci[icoeff];
        }
      }
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
          bool isComputeParticleVirial,
          bool isHybrid>
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
      !isComputeProcess_d2Edr2 &&
      !isHybrid)
  {
    return ier;
  }

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

  // If it is a hybrid style
  if (isHybrid)
  {
    double dEidr_ZBL;
    double dEidrByR_ZBL;
    double dEidrByR_TABLE;
    double a;
    double b;
    double fraction;
    int itable;

    TABLE const *tb;
    union_int_float_t rsq_lookup;

    int numnei = 0;
    int const *n1atom = NULL;

    // Loop over all the contributing particles
    for (int i = 0, contributing_index = 0; i < cachedNumberOfParticles_; ++i)
    {
      if (!particleContributing[i])
      {
        continue;
      }

      // Get the species index for atom i
      int const iSpecies = particleSpeciesCodes[i];

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

      if (snapflag[iSpecies])
      {
        // Get the iSpecies cutoff
        double const radi = radelem[iSpecies];

        // noutside, is the number of neighbors of I outside the cutoff

        // note Rij sign convention => dU/dRij = dU/dRj = -dU/dRi

        // Setup loop over neighbors of current particle
        for (int n = 0, noutside = numnei; n < numnei; ++n)
        {
          // Index of the neighbor atom
          int const j = n1atom[n];

          // Get the species index for atom j
          int const jSpecies = particleSpeciesCodes[j];

          double const dx = coordinates[j][0] - xi;
          double const dy = coordinates[j][1] - yi;
          double const dz = coordinates[j][2] - zi;

          if (snapflag[jSpecies])
          {
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
              snap->element[ninside] = jSpecies;
              ++ninside;
            }
            else
            {
              --noutside;
              snap->rij(noutside, 0) = dx;
              snap->rij(noutside, 1) = dy;
              snap->rij(noutside, 2) = dz;
              snap->inside[noutside] = j;
              snap->element[noutside] = jSpecies;
            }
          }
          else
          {
            --noutside;
            snap->rij(noutside, 0) = dx;
            snap->rij(noutside, 1) = dy;
            snap->rij(noutside, 2) = dz;
            snap->inside[noutside] = j;
            snap->element[noutside] = jSpecies;
          }
        }

        // compute Ui, Yi for atom i
        if (chemflag)
        {
          snap->compute_ui(ninside, iSpecies);
        }
        else
        {
          snap->compute_ui(ninside, 0);
        }

        // Get the 1D view to the 2D beta array at row i
        auto betai = beta.data_1D(contributing_index);

        // Get the pointer to the beta array of data for atom i
        double const *const bi = const_cast<double *>(betai.data());

        snap->compute_yi(bi);
      }
      else
      {
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

          snap->rij(n, 0) = dx;
          snap->rij(n, 1) = dy;
          snap->rij(n, 2) = dz;
          snap->inside[n] = j;
          snap->element[n] = jSpecies;
        }
      }

      // Compute contribution to force, etc.

      // For the neighbors of particle i within the cutoff:
      // Compute deidrj = dEi/dRj = -dEi/dRi add to Fi, subtract from Fj

      VectorOfSizeDIM deidrj;

      // Setup loop over neighbors of particle i
      for (int n = 0; n < numnei; ++n)
      {
        // Index of the neighbor atom
        int const j = snap->inside[n];

        // Get the species index for atom j
        int const jSpecies = particleSpeciesCodes[j];

        // If particle j is contributing or not
        int const jContrib = particleContributing[j];

        bool const lsnap = n < ninside;

        bool const lhalf = !(jContrib && (j < i));

        if (!lsnap)
        {
          if (!lhalf)
          {
            continue;
          }
          if (setflag(iSpecies, jSpecies) == HYBRIDSTYLE::NONE)
          {
            continue;
          }
        }

        bool lzbl = lhalf ? (setflag(iSpecies, jSpecies) == HYBRIDSTYLE::ZBL) : false;

        bool ltable = lhalf ? (setflag(iSpecies, jSpecies) == HYBRIDSTYLE::TABLE) : false;
        if (ltable)
        {
          tb = &tables[tableNumber(iSpecies, jSpecies)];
        }

        // Get the 1D view to the 2D snap->rij array at row n
        auto rij = snap->rij.data_1D(n);

        // Get the pointer to the rij_const array of data
        double const *const rij_const = const_cast<double *>(rij.data());

        double const rsq = rij_const[0] * rij_const[0] + rij_const[1] * rij_const[1] + rij_const[2] * rij_const[2];
        double const rrsq = std::sqrt(rsq);

        if (lzbl)
        {
          lzbl = (rsq < zbl->cut_globalsq);
        }

        if (ltable)
        {
          ltable = (rrsq < tb->cut) ? (rsq > tb->innersq) : false;
        }

        if (lsnap)
        {
          if (chemflag)
          {
            snap->compute_duidrj(rij_const, snap->wj[n], snap->rcutij[n], n, snap->element[n]);
          }
          else
          {
            snap->compute_duidrj(rij_const, snap->wj[n], snap->rcutij[n], n, 0);
          }

          snap->compute_deidrj(deidrj);
        }

        if (lzbl)
        {
          // Compute dphi
          if (isComputeForces ||
              isComputeProcess_dEdr ||
              isComputeVirial ||
              isComputeParticleVirial)
          {
            dEidr_ZBL = zbl->dzbldr(rrsq, iSpecies, jSpecies);

            if (rsq > zbl->cut_innersq)
            {
              double const t = rrsq - zbl->cut_inner;
              double const fswitch = t * t * (zbl->sw1(iSpecies, jSpecies) + zbl->sw2(iSpecies, jSpecies) * t);
              dEidr_ZBL += fswitch;
            }

            if (!jContrib)
            {
              dEidr_ZBL *= 0.5;
            }

            dEidrByR_ZBL = dEidr_ZBL / rrsq;
          }
          else
          {
            dEidr_ZBL = 0.0;
            dEidrByR_ZBL = 0.0;
          }
        }

        if (ltable)
        {
          int const tlm1 = tb->tableLength - 1;

          if (tb->tableStyle == TABLESTYLE::LOOKUP)
          {
            itable = static_cast<int>((rsq - tb->innersq) * tb->invdelta);
            if (itable >= tlm1)
            {
              LOG_ERROR("Pair distance (" +
                        std::to_string(rrsq) +
                        ") > table outer cutoff between types (" +
                        hybridStyleSpeciesNames[iSpecies] +
                        ", " +
                        hybridStyleSpeciesNames[jSpecies] + ")");
              return true;
            }
            dEidrByR_TABLE = tb->f[itable];
          }
          else if (tb->tableStyle == TABLESTYLE::LINEAR)
          {
            itable = static_cast<int>((rsq - tb->innersq) * tb->invdelta);
            if (itable >= tlm1)
            {
              LOG_ERROR("Pair distance (" +
                        std::to_string(rrsq) +
                        ") > table outer cutoff between types (" +
                        hybridStyleSpeciesNames[iSpecies] +
                        ", " +
                        hybridStyleSpeciesNames[jSpecies] + ")");
              return true;
            }
            fraction = (rsq - tb->rsq[itable]) * tb->invdelta;
            dEidrByR_TABLE = tb->f[itable] + fraction * tb->df[itable];
          }
          else if (tb->tableStyle == TABLESTYLE::SPLINE)
          {
            itable = static_cast<int>((rsq - tb->innersq) * tb->invdelta);
            if (itable >= tlm1)
            {
              LOG_ERROR("Pair distance (" +
                        std::to_string(rrsq) +
                        ") > table outer cutoff between types (" +
                        hybridStyleSpeciesNames[iSpecies] +
                        ", " +
                        hybridStyleSpeciesNames[jSpecies] + ")");
              return true;
            }
            b = (rsq - tb->rsq[itable]) * tb->invdelta;
            a = 1.0 - b;
            dEidrByR_TABLE = a * tb->f[itable] + b * tb->f[itable + 1] + ((a * a * a - a) * tb->f2[itable] + (b * b * b - b) * tb->f2[itable + 1]) * tb->deltasq6;
          }
          else if (tb->tableStyle == TABLESTYLE::BITMAP)
          {
            rsq_lookup.f = rsq;
            itable = rsq_lookup.i & tb->nmask;
            itable >>= tb->nshiftbits;
            fraction = (rsq_lookup.f - tb->rsq[itable]) * tb->drsq[itable];
            dEidrByR_TABLE = tb->f[itable] + fraction * tb->df[itable];
          }

          dEidrByR_TABLE *= jContrib ? -1.0 : -0.5;
        }

        // Contribution to forces
        if (isComputeForces)
        {
          if (lsnap)
          {
            forces[i][0] += deidrj[0];
            forces[i][1] += deidrj[1];
            forces[i][2] += deidrj[2];

            forces[j][0] -= deidrj[0];
            forces[j][1] -= deidrj[1];
            forces[j][2] -= deidrj[2];
          }

          if (lzbl)
          {
            forces[i][0] += dEidrByR_ZBL * rij_const[0];
            forces[i][1] += dEidrByR_ZBL * rij_const[1];
            forces[i][2] += dEidrByR_ZBL * rij_const[2];

            forces[j][0] -= dEidrByR_ZBL * rij_const[0];
            forces[j][1] -= dEidrByR_ZBL * rij_const[1];
            forces[j][2] -= dEidrByR_ZBL * rij_const[2];
          }

          if (ltable)
          {
            forces[i][0] += dEidrByR_TABLE * rij_const[0];
            forces[i][1] += dEidrByR_TABLE * rij_const[1];
            forces[i][2] += dEidrByR_TABLE * rij_const[2];

            forces[j][0] -= dEidrByR_TABLE * rij_const[0];
            forces[j][1] -= dEidrByR_TABLE * rij_const[1];
            forces[j][2] -= dEidrByR_TABLE * rij_const[2];
          }
        }

        if (isComputeEnergy || isComputeParticleEnergy)
        {
          if (lzbl)
          {
            double phi = zbl->e_zbl(rrsq, iSpecies, jSpecies);

            phi += zbl->sw5(iSpecies, jSpecies);

            if (rsq > zbl->cut_innersq)
            {
              double const t = rrsq - zbl->cut_inner;
              double const eswitch = t * t * t * (zbl->sw3(iSpecies, jSpecies) + zbl->sw4(iSpecies, jSpecies) * t);
              phi += eswitch;
            }

            if (isComputeEnergy)
            {
              *energy += jContrib ? phi : 0.5 * phi;
            }

            if (isComputeParticleEnergy)
            {
              particleEnergy[i] += 0.5 * phi;
              if (jContrib)
              {
                particleEnergy[j] += 0.5 * phi;
              }
            }
          }

          if (ltable)
          {
            double phi;

            if (tb->tableStyle == TABLESTYLE::LOOKUP)
            {
              phi = tb->e[itable];
            }
            else if (tb->tableStyle == TABLESTYLE::LINEAR)
            {
              phi = tb->e[itable] + fraction * tb->de[itable];
            }
            else if (tb->tableStyle == TABLESTYLE::SPLINE)
            {
              phi = a * tb->e[itable] +
                    b * tb->e[itable + 1] +
                    ((a * a * a - a) * tb->e2[itable] +
                     (b * b * b - b) * tb->e2[itable + 1]) *
                        tb->deltasq6;
            }
            else if (tb->tableStyle == TABLESTYLE::BITMAP)
            {
              phi = tb->e[itable] + fraction * tb->de[itable];
            }

            if (isComputeEnergy)
            {
              *energy += jContrib ? phi : 0.5 * phi;
            }

            if (isComputeParticleEnergy)
            {
              particleEnergy[i] += 0.5 * phi;
              if (jContrib)
              {
                particleEnergy[j] += 0.5 * phi;
              }
            }
          }
        }

        if (isComputeProcess_dEdr)
        {
          double dedr = lsnap ? std::sqrt(deidrj[0] * deidrj[0] + deidrj[1] * deidrj[1] + deidrj[2] * deidrj[2]) : 0.0;

          if (lzbl)
          {
            dedr += dEidr_ZBL;
          }

          if (ltable)
          {
            dedr += dEidrByR_TABLE * rrsq;
          }

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

          if (lsnap)
          {
            v[0] = rij_const[0] * deidrj[0];
            v[1] = rij_const[1] * deidrj[1];
            v[2] = rij_const[2] * deidrj[2];
            v[3] = rij_const[1] * deidrj[2];
            v[4] = rij_const[0] * deidrj[2];
            v[5] = rij_const[0] * deidrj[1];
          }
          else
          {
            v[0] = 0.0;
            v[1] = 0.0;
            v[2] = 0.0;
            v[3] = 0.0;
            v[4] = 0.0;
            v[5] = 0.0;
          }

          if (lzbl)
          {
            v[0] += rij_const[0] * rij_const[0] * dEidrByR_ZBL;
            v[1] += rij_const[1] * rij_const[1] * dEidrByR_ZBL;
            v[2] += rij_const[2] * rij_const[2] * dEidrByR_ZBL;
            v[3] += rij_const[1] * rij_const[2] * dEidrByR_ZBL;
            v[4] += rij_const[0] * rij_const[2] * dEidrByR_ZBL;
            v[5] += rij_const[0] * rij_const[1] * dEidrByR_ZBL;
          }

          if (ltable)
          {
            v[0] += rij_const[0] * rij_const[0] * dEidrByR_TABLE;
            v[1] += rij_const[1] * rij_const[1] * dEidrByR_TABLE;
            v[2] += rij_const[2] * rij_const[2] * dEidrByR_TABLE;
            v[3] += rij_const[1] * rij_const[2] * dEidrByR_TABLE;
            v[4] += rij_const[0] * rij_const[2] * dEidrByR_TABLE;
            v[5] += rij_const[0] * rij_const[1] * dEidrByR_TABLE;
          }

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

      if (snapflag[iSpecies])
      {
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
          {
            phi += Ci[icoeff] * Bi[icoeff];
          }

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
      } // snapflag[iSpecies]
    }   // End of loop over contributing particles
  }     // If it is a hybrid style
  else
  // If it is not a hybrid style
  {
    int numnei = 0;
    int const *n1atom = NULL;

    // Loop over all the contributing particles
    for (int i = 0, contributing_index = 0; i < cachedNumberOfParticles_; ++i)
    {
      if (!particleContributing[i])
      {
        continue;
      }

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
          snap->element[ninside] = jSpecies;
          ++ninside;
        }
      }

      // compute Ui, Yi for atom i
      {
        if (chemflag)
        {
          snap->compute_ui(ninside, iSpecies);
        }
        else
        {
          snap->compute_ui(ninside, 0);
        }

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

        if (chemflag)
        {
          snap->compute_duidrj(rij_const, snap->wj[n], snap->rcutij[n], n, snap->element[n]);
        }
        else
        {
          snap->compute_duidrj(rij_const, snap->wj[n], snap->rcutij[n], n, 0);
        }

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
          double const dedr = std::sqrt(deidrj[0] * deidrj[0] + deidrj[1] * deidrj[1] + deidrj[2] * deidrj[2]);

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
        {
          phi += Ci[icoeff] * Bi[icoeff];
        }

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
  }   // If it is not a hybrid style

  // everything is good
  return false;
}
#undef KIM_LOGGER_OBJECT_NAME

#undef ONE
#undef MAXLINE
#undef HELPER_LOG_ERROR
