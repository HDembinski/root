// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/MnParameterScan.h"
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnPrint.h"

namespace ROOT {

   namespace Minuit2 {


MnParameterScan::MnParameterScan(const FCNBase& fcn, const MnUserParameters& par) : fFCN(fcn), fParameters(par), fAmin(fcn(par.Params())) {}

MnParameterScan::MnParameterScan(const FCNBase& fcn, const MnUserParameters& par, double fval) : fFCN(fcn), fParameters(par), fAmin(fval) {}

std::vector<std::pair<double, double> > MnParameterScan::operator()(unsigned int par, unsigned int maxsteps, double low, double high) {
   // do the scan for parameter par between low and high values

   MnPrint print("MnParameterScan");

   //if(maxsteps > 101) maxsteps = 101;
   std::vector<std::pair<double, double> > result; result.reserve(maxsteps+1);
   std::vector<double> params = fParameters.Params();
   result.emplace_back(params[par], fAmin);

   if(low > high) return result;
   if(maxsteps < 2) return result;

   const auto& p = fParameters.Parameter(par);

   if(low == 0. && high == 0.) {
      low = p.HasLowerLimit() ? fParameters.Parameter(par).LowerLimit() : params[par] - 2.*fParameters.Error(par);
      high = p.HasUpperLimit() ? fParameters.Parameter(par).UpperLimit() : params[par] + 2.*fParameters.Error(par);
   }

   double x0 = low;
   double stp = (high - low)/double(maxsteps - 1);

   print.Info("Bounds", low, high, "Nsteps", maxsteps);

   for(unsigned int i = 0; i < maxsteps; i++) {
      params[par] = x0 + double(i)*stp;
      double fval = fFCN(params);
      if(fval < fAmin) {
         fParameters.SetValue(par, params[par]);
         fAmin = fval;
      }
      print.Debug([&](std::ostream& os) {
         for (const auto& p : params)
           os << p << " ";
         os << ": " << fval;
      });
      result.emplace_back(params[par], fval);
   }

   print.Info([&](std::ostream& os) {
      os << "New minimum " <<  fAmin << " at parameter values";
      for (const auto& p : fParameters.Parameters())
        os << " " << p.Value();
   });

   return result;
}

   }  // namespace Minuit2

}  // namespace ROOT
