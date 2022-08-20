/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>

#include "spxdefines.h"
#include "spxsolver.h"
#include "spxpricer.h"
#include "spxratiotester.h"
#include "spxstarter.h"
#include "slinsolver.h"
#include "slufactor.h"

namespace soplex
{
bool SPxSolver::writeState(
   const char*    filename,
   const NameSet* rowNames,
   const NameSet* colNames,
   const bool cpxFormat
   ) const
{

   std::string ofname;
   std::ofstream ofs;

   // write parameter settings
   ofname = std::string(filename) + ".set";
   ofs.open(ofname.c_str());
   if (!ofs)
      return false;

   ofs << "# SoPlex version " << SOPLEX_VERSION / 100
       << "." << (SOPLEX_VERSION / 10) % 10
       << "." << SOPLEX_VERSION % 10
       << "." << SOPLEX_SUBVERSION << std::endl << std::endl;
   ofs << "# run SoPlex as follows:" << std::endl;
   ofs << "# bin/soplex --loadset=spxcheck.set --readbas=spxcheck.bas spxcheck.mps\n" << std::endl;
   ofs << "int:representation = " << ( rep() == SPxSolver::COLUMN ? "1" : "2" ) << std::endl;
   ofs << "int:factor_update_max = " << basis().getMaxUpdates() << std::endl;
   ofs << "int:pricer = ";
   if (!strcmp(pricer()->getName(), "Auto"))
      ofs << " 0" << std::endl;
   else if (!strcmp(pricer()->getName(), "Dantzig"))
      ofs << "1" << std::endl;
   else if (!strcmp(pricer()->getName(), "ParMult"))
      ofs << "2" << std::endl;
   else if (!strcmp(pricer()->getName(), "Devex"))
      ofs << "3" << std::endl;
   else if (!strcmp(pricer()->getName(), "Steep"))
      ofs << "4" << std::endl;
   else if (!strcmp(pricer()->getName(), "SteepEx"))
      ofs << "5" << std::endl;
   ofs << "int:ratiotester = ";
   if (!strcmp(ratiotester()->getName(), "Default"))
      ofs << "0" << std::endl;
   else if (!strcmp(ratiotester()->getName(), "Harris"))
      ofs << "1" << std::endl;
   else if (!strcmp(ratiotester()->getName(), "Fast"))
      ofs << "2" << std::endl;
   else if (!strcmp(ratiotester()->getName(), "Bound Flipping"))
      ofs << "3" << std::endl;
   ofs << "real:feastol = " << feastol() << std::endl;
   ofs << "real:opttol = " << opttol() << std::endl;
   ofs << "real:epsilon_zero = " << epsilon() << std::endl;
   ofs << "real:infty = " << infinity << std::endl;
   ofs << "uint:random_seed = " << random.getSeed() << std::endl;
   ofs.close();

   // write LP
   ofname = std::string(filename) + ".mps";
   ofs.open(ofname.c_str());
   if (!ofs)
      return false;

   writeMPS(ofs, rowNames, colNames, NULL);
   ofs.close();

   // write basis
   ofname = std::string(filename) + ".bas";
   return writeBasisFile(ofname.c_str(), rowNames, colNames, cpxFormat);
}

} // namespace soplex
