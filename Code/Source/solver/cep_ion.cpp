// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "cep_ion.h"

#include "all_fun.h"
#include "post.h"
#include "utils.h"
#include <math.h>

namespace cep_ion {

/// @brief Modifies:
/// \code {.cpp}
///   cep_mod.Xion
/// \endcode
//
void cep_init(Simulation* simulation)
{
  using namespace consts;
  auto& com_mod = simulation->com_mod;

  #define n_debug_cep_init 
  #ifdef debug_cep_init 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  auto& cm = com_mod.cm;
  auto& cep_mod = simulation->cep_mod;
  const int nsd = com_mod.nsd;
  const int tnNo = com_mod.tnNo;
  const int nXion = cep_mod.nXion;
  #ifdef debug_cep_init 
  dmsg << "tnNo: " << tnNo;
  dmsg << "nXion: " << nXion;
  #endif

  for (auto& eq : com_mod.eq) {
    if (eq.phys != EquationType::phys_CEP) {
      continue;
    }

    if (com_mod.dmnId.size() != 0) {
      Vector<double> sA(tnNo); 
      Array<double> sF(nXion,tnNo);

      for (int a = 0; a < tnNo; a++) {
        if (!all_fun::is_domain(com_mod, eq, a, EquationType::phys_CEP)) {
          continue;
        }
        for (int iDmn = 0; iDmn < eq.nDmn; iDmn++) {
          auto cPhys = eq.dmn[iDmn].phys;
          int dID = eq.dmn[iDmn].Id;
          if ((cPhys != EquationType::phys_CEP) || !utils::btest(com_mod.dmnId(a),dID)) {
            continue;
          }
          int nX = eq.dmn[iDmn].cep.nX;
          int nG = eq.dmn[iDmn].cep.nG;
          int imyo = eq.dmn[iDmn].cep.imyo;

          Vector<double> Xl(nX); 
          Vector<double> Xgl(nG);

          cep_init_l(cep_mod, eq.dmn[iDmn].cep, nX, nG, Xl, Xgl);

          // Initialize Land model states for this node (only for TTP model)
          if (eq.dmn[iDmn].cep.cepType == ElectrophysiologyModelType::TTP) {
            cep_init_land_l(cep_mod, a);
          }

          sA(a) = sA(a) + 1.0;

          for (int i = 0; i < nX; i++) {
            sF(i,a)  = sF(i,a) + Xl(i);
          }

          for (int i = 0; i < nG; i++) {
            sF(i+nX,a) = sF(i+nX,a) + Xgl(i);
          }
        }
      }

      all_fun::commu(com_mod, sA);
      all_fun::commu(com_mod, sF);

      for (int a = 0; a < tnNo; a++) {
        if (!utils::is_zero(sA(a))) {
          for (int i = 0; i < cep_mod.Xion.nrows(); i++) {
            cep_mod.Xion(i,a) = sF(i,a) / sA(a);
          }
        }
      }

    } else {
      for (int a = 0; a < tnNo; a++) { 
        if (!all_fun::is_domain(com_mod, eq, a, EquationType::phys_CEP)) {
          continue;
        }
        int nX = eq.dmn[0].cep.nX;
        int nG = eq.dmn[0].cep.nG;
        Vector<double> Xl(nX); 
        Vector<double> Xgl(nG);

        cep_init_l(cep_mod, eq.dmn[1].cep, nX, nG, Xl, Xgl);

        // Initialize Land model states for this node (only for TTP model)
        if (eq.dmn[0].cep.cepType == ElectrophysiologyModelType::TTP) {
          cep_init_land_l(cep_mod, a);
        }

        for (int i = 0; i < nX; i++) {
          cep_mod.Xion(i,a) = Xl(i);
        }
        for (int i = 0; i < nG; i++) {
          cep_mod.Xion(i+nX,a) = Xgl(i);
        }
      }
    }
  }
}

//------------
// cep_init_l
//------------
//
void cep_init_l(CepMod& cep_mod, cepModelType& cep, int nX, int nG, Vector<double>& X, Vector<double>& Xg)
{
  switch (cep.cepType) {

    case ElectrophysiologyModelType::AP:
      cep_mod.ap.init(nX, X);
    break;

    case ElectrophysiologyModelType::BO:
      cep_mod.bo.init(nX, X);
    break;

    case ElectrophysiologyModelType::FN:
      cep_mod.fn.init(nX, X);
    break;

    case ElectrophysiologyModelType::TTP:
      cep_mod.ttp.init(cep.imyo, nX, nG, X, Xg);
    break;
  }
}

//-------------
// cep_init_land_l
//-------------
// Initialize Land model state variables for a specific node
//
void cep_init_land_l(CepMod& cep_mod, const int nodeId)
{
  // Initialize Land model state variables to reasonable resting values
  // Based on Land 2016 model, these are typical initial values for resting state
  cep_mod.Y_land(0, nodeId) = 0.99;   // XB - bound cross-bridges (resting = 1)
  cep_mod.Y_land(1, nodeId) = 0.0;   // XW - weakly bound cross-bridges (resting = 0)
  cep_mod.Y_land(2, nodeId) = 0.0;   // CaTRPN - troponin bound to calcium (resting = 0)
  cep_mod.Y_land(3, nodeId) = 0.0;   // XS - strongly bound cross-bridges (resting = 0)
  cep_mod.Y_land(4, nodeId) = 0.0;   // ZETAW - weakly bound cross-bridge distortion (resting = 0)
  cep_mod.Y_land(5, nodeId) = 0.0;   // ZETAS - strongly bound cross-bridge distortion (resting = 0)

}

//-----------
// cep_integ
//-----------
// State variable integration.
//
void cep_integ(Simulation* simulation, const int iEq, const int iDof, const Array<double>& Dg)
{
  static bool IPASS = true;

  using namespace consts;

  auto& com_mod = simulation->com_mod;
  auto& cm_mod = simulation->cm_mod;

  #define n_debug_cep_integ 
  #ifdef debug_cep_integ
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  auto& cm = com_mod.cm;
  int tnNo = com_mod.tnNo;
  double dt = com_mod.dt;
  double time = com_mod.time;

  auto& cep_mod = simulation->cep_mod;
  auto& cem = cep_mod.cem;
  auto& eq = com_mod.eq[iEq];

  auto& Yo = com_mod.Yo;
  auto& Dn = com_mod.Dn;
  auto& Xion = cep_mod.Xion;
  int nXion = cep_mod.nXion;

  Vector<double> I4f(tnNo);

  #ifdef debug_cep_integ
  dmsg << "cem.cpld: " << cem.cpld;
  dmsg << "time: " << time;
  #endif

  // Electromechanics: get fiber stretch for stretch activated currents
  // Compute lambda at both time n (old) and time n+1 (new) for stabilization
  //
  // std::cout << "Finding structural equation" << std::endl;
  if (cem.cpld) {
    for (int iM = 0; iM < com_mod.nMsh; iM++) {
      auto& msh = com_mod.msh[iM];

      if (msh.nFn != 0) {
        Vector<double> sA(msh.nNo);
        post::fib_strech(simulation, iEq, msh, Dg, sA);
        for (int a = 0; a < msh.nNo; a++) {
          int Ac = msh.gN(a);
          I4f(Ac) = sA(a);
        }
      }
    }
  }

  //  Ignore first pass as Xion is already initialized
  if (IPASS) {
    IPASS = false;

  // Copy action potential after diffusion as first state variable
  } else {
    for (int Ac = 0; Ac < tnNo; Ac++) {
      Xion(0,Ac) = Yo(iDof,Ac);
    }
  }

  // Integrate electric potential based on cellular activation model
  //
  if (com_mod.dmnId.size() != 0) {
    Vector<double> sA(tnNo); 
    Array<double> sF(nXion,tnNo); 
    Vector<double> sY(tnNo);
    Vector<double> sTa(tnNo);
    Vector<double> sKa(tnNo);
    sA = 0.0;
    sF = 0.0;
    sY = 0.0;
    sTa = 0.0;
    sKa = 0.0;

    for (int Ac = 0; Ac < tnNo; Ac++) {
      if (!all_fun::is_domain(com_mod, eq, Ac, Equation_CEP)) {
        continue;
      }

      for (int iDmn = 0; iDmn < eq.nDmn; iDmn++) {
        auto& dmn = eq.dmn[iDmn];
        auto cPhys = dmn.phys;
        int dID = dmn.Id;

        if (cPhys != Equation_CEP || !utils::btest(com_mod.dmnId(Ac),dID)) {
          continue;
	}

        int nX = dmn.cep.nX;
        int nG = dmn.cep.nG;
        #ifdef debug_cep_integ
        dmsg << "nX: " << nX ;
        dmsg << "nG: " << nG ;
        #endif

        auto Xl = Xion.rows(0,nX-1,Ac); // extract column vector of ionic state variables per node

        // [NOTE] nG can be 0.
        Vector<double> Xgl;
        if (nG != 0) {
          Xgl.resize(nG);
          for (int i = 0; i < nG; i++) {
            Xgl(i) = Xion(i+nX,Ac);
          }
        }

        // Extract Land model states for this node (only for TTP model)
        Vector<double> Y_land_node;
        if (dmn.cep.cepType == ElectrophysiologyModelType::TTP) {
          Y_land_node = cep_mod.Y_land.col(Ac);
        }
        double yl = 0.0;
        if (cem.cpld) {
          yl = cem.Ya(Ac);
        }

        cep_integ_l(cep_mod, dmn.cep, nX, nG, Xl, Xgl, time-dt, yl, I4f(Ac), dt);

        sA(Ac) = sA(Ac) + 1.0;
        for (int i = 0; i < nX; i++) {
          sF(i,Ac) += Xl(i);
        }

        for (int i = 0; i < nG; i++) {
          sF(nX+i,Ac) += Xgl(i);
        }

        if (cem.cpld) {
          double Ta_dyne = Ta_node *1.0e-3;
          double Ka_dyne = Ka_node *1.0e-3;
          yl = Ta_dyne;
          sTa(Ac) += Ta_dyne;
          sKa(Ac) += Ka_dyne;
          sY(Ac) = sY(Ac) + yl;
#ifdef debug_land_actv
          if (Ac < 5 && cm.mas(cm_mod)) {
            std::cout << "[cep_integ] node=" << Ac
                      << ", lambda_prev=" << cem.lambda_prev(Ac)
                      << ", lambda_curr=" << cem.lambda_curr(Ac)
                      << ", dlambda_dt=" << dlambda_dt_vec(Ac)
                      << ", Ta=" << Ta_dyne
                      << ", Ka=" << Ka_dyne
                      << ", tension=" << yl
                      << std::endl;
          }
#endif
        }

        // Copy updated Land model states back to global array (only for TTP model)
        if (dmn.cep.cepType == ElectrophysiologyModelType::TTP && Y_land_node.size() > 0) {
          cep_mod.Y_land.set_col(Ac, Y_land_node);
        }
      }
    }

    all_fun::commu(com_mod, sA);
    all_fun::commu(com_mod, sF);

    if (cem.cpld) {
      all_fun::commu(com_mod, sY);
      all_fun::commu(com_mod, sTa);
      all_fun::commu(com_mod, sKa);
    }

    // Debug: Print accumulated values for node 0 from all processors
    // Get values from the global arrays after communication
    double global_c_Ca = 0.0;
    double global_I4f = 0.0;
    double global_I4fRate = 0.0;
    double global_ActiveTension = 0.0;
    
    if (com_mod.dmnId.size() != 0) {
      // For domain-based case, get from global arrays
      if (com_mod.msh[0].gN(0) < com_mod.tnNo) {
        int global_node_0 = com_mod.msh[0].gN(0);
        global_c_Ca = Xion(3, global_node_0);
        global_I4f = cem.lambda_curr(global_node_0);  // Current lambda for diagnostics
        global_I4fRate = dlambda_dt_vec(global_node_0);
        if (cem.cpld) {
          global_ActiveTension = cem.Ya(global_node_0);
        }
      }
    } else {
      // For non-domain case, get directly
      global_c_Ca = Xion(3, 0);
      global_I4f = cem.lambda_curr(0);  // Current lambda for diagnostics
      global_I4fRate = dlambda_dt_vec(0);
      if (cem.cpld) {
        global_ActiveTension = cem.Ya(0);
      }
    }
    
    // Only print from master processor
    // if (cm.mas(cm_mod)) {
    //   std::cout << "[cep_integ] Node 0 (global): c_Ca=" << global_c_Ca 
    //             << ", I4f=" << global_I4f << ", I4fRate=" << global_I4fRate 
    //             << ", ActiveTension=" << global_ActiveTension << " dyne/cm²" << std::endl;
    // }

    for (int Ac = 0; Ac < tnNo; Ac++) {
      if (!utils::is_zero(sA(Ac))) {
        Xion.set_col(Ac, sF.col(Ac) / sA(Ac));
        if (cem.cpld) {
          cem.Ta(Ac) = sTa(Ac) / sA(Ac);
          cem.Ka(Ac) = sKa(Ac) / sA(Ac);
          cem.Ya(Ac) = cem.Ta(Ac);
        }
      }
    }

  } else {
    for (int Ac = 0; Ac < tnNo; Ac++) {
      if (!all_fun::is_domain(com_mod, eq, Ac, Equation_CEP)) {
        continue;
      }

      int nX = eq.dmn[0].cep.nX;
      int nG = eq.dmn[0].cep.nG;
      auto Xl = Xion.rows(0,nX-1,Ac);
      auto Xgl = Xion.rows(nX,nX+nG-1,Ac);

      // Extract Land model states for this node (only for TTP model)
      Vector<double> Y_land_node;
      if (eq.dmn[0].cep.cepType == ElectrophysiologyModelType::TTP) {
        Y_land_node = cep_mod.Y_land.col(Ac);
      }

      double yl = 0.0;
      if (cem.cpld) {
        yl = cem.Ya(Ac);
      }

      cep_integ_l(cep_mod, eq.dmn[0].cep, nX, nG, Xl, Xgl, time-dt, yl, I4f(Ac), dt);

      for (int i = 0; i < nX; i++) {
        Xion(i,Ac) = Xl(i);
      }

      for (int i = 0; i < nG; i++) {
        Xion(nX+i,Ac) = Xgl(i);
      }

      // Copy updated Land model states back to global array (only for TTP model)
      if (eq.dmn[0].cep.cepType == ElectrophysiologyModelType::TTP && Y_land_node.size() > 0) {
        cep_mod.Y_land.set_col(Ac, Y_land_node);
      }

      if (cem.cpld) {
        double Ta_dyne = Ta_node *1.0e-3;
        double Ka_dyne = Ka_node *1.0e-3;
        cem.Ta(Ac) = Ta_dyne;
        cem.Ka(Ac) = Ka_dyne;
        cem.Ya(Ac) = Ta_dyne;
        yl = Ta_dyne;
#ifdef debug_land_actv
        if (Ac < 5 && cm.mas(cm_mod)) {
          std::cout << "[cep_integ] node=" << Ac
                    << ", lambda_prev=" << cem.lambda_prev(Ac)
                    << ", lambda_curr=" << cem.lambda_curr(Ac)
                    << ", dlambda_dt=" << dlambda_dt_vec(Ac)
                    << ", Ta=" << Ta_dyne
                    << ", Ka=" << Ka_dyne
                    << ", tension=" << tension
                    << std::endl;
        }
#endif
      }
    }

    // Debug: Print values for node 0 in non-domain case
    // Only print from master processor
      if (cm.mas(cm_mod)) {
        double global_c_Ca = Xion(3, 0);
        double global_I4f = cem.lambda_curr(0);  // Current lambda for diagnostics
        double global_I4fRate = dlambda_dt_vec(0);
        double global_ActiveTension = cem.cpld ? cem.Ya(0) : 0.0;
        
        std::cout << "[cep_integ] Node 0 (non-domain): c_Ca=" << global_c_Ca 
                  << ", I4f=" << global_I4f << ", I4fRate=" << global_I4fRate 
                  << ", ActiveTension=" << global_ActiveTension << " dyne/cm²" << std::endl;
      }
  }

  for (int Ac = 0; Ac < tnNo; Ac++) {
    Yo(iDof,Ac) = Xion(0,Ac);   // copy back to global solution vector (action potential only)
  }
}

//-------------
// cep_integ_l
//-------------
// Integrate local electrophysiology variables from t1 to t1+dt. Also
// integrate excitation-activation variables form coupled electro-
// mechanics. The equations are integrated at domain nodes.
//
// Integrates electrophysiology variables and computes active tension using analytical dlambda_dt
//
void cep_integ_l(CepMod& cep_mod, cepModelType& cep, int nX, int nG, Vector<double>& X, Vector<double>& Xg, 
    const double t1, double& yl, const double I4f, const double dt)
{
  using namespace consts;

  #define n_debug_cep_integ_l
  #ifdef debug_cep_integ_l
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  auto& cem = cep_mod.cem;

  // Feedback coefficient for stretch-activated-currents
  // Use lambda_new (predicted state) for SAC coupling
  double Ksac = 0.0;
  if (lambda_new > 1.0) {
     Ksac = cep.Ksac * (lambda_new - 1.0);
  } else {
     Ksac = 0.0;
  }

  // Total time steps
  int nt = static_cast<int>(dt/cep.dt); // for pure EP simulations, cep.dt = dt, so nt = 1
  // but if coupled with mechanics, cep.dt can be smaller than dt (which comes from com_mod.dt) 

  // External stimulus duration
  int icl = static_cast<int>(fmax(floor(t1/cep.Istim.CL),0.0));
  double Ts = cep.Istim.Ts + static_cast<double>(icl)*cep.Istim.CL;
  double Te = Ts + cep.Istim.Td;
  double eps = std::numeric_limits<double>::epsilon();

  #ifdef debug_cep_integ_l
  dmsg << "nt: " << nt;
  dmsg << "Ksac: " << Ksac;
  dmsg << "icl: " << icl;
  dmsg << "Ts: " << Ts;
  dmsg << "Te: " << Te;
  dmsg << "cep.cepType: " << cep.cepType;
  dmsg << "cep.odes.tIntTyp: " << cep.odes.tIntType;
  #endif

  Ta_out = 0.0;
  Ka_out = 0.0;

  switch (cep.cepType) {
    case ElectrophysiologyModelType::AP: {
      Vector<int> IPAR(2); 
      Vector<double> RPAR(2);
      IPAR(0) = cep.odes.maxItr;
      IPAR(1) = 0;
      RPAR(0) = cep.odes.absTol;
      RPAR(1) = cep.odes.relTol;

      switch (cep.odes.tIntType) {
        case TimeIntegratioType::FE: {
          for (int i = 0; i < nt; i++) {
            double t = t1 + static_cast<double>(i) * cep.dt;
            double Istim;
            if (t >= Ts-eps &&  t <= Te+eps) {
              Istim = cep.Istim.A;
            } else {
              Istim = 0.0;
            }
            cep_mod.ap.integ_fe(nX, X, t, cep.dt, Istim, Ksac);
  
            // Electromechanics excitation-activation
            if (cem.aStress) {
              double epsX;
              cep_mod.ap.actv_strs(X(0), cep.dt, yl, epsX);
            }
          }
        } break; 

        case TimeIntegratioType::RK4: {
          for (int i = 0; i < nt; i++) {
            double t = t1 + static_cast<double>(i) * cep.dt;
            double Istim;
            if (t >= Ts-eps &&  t <= Te+eps) {
              Istim = cep.Istim.A;
            } else {
              Istim = 0.0;
            }
            cep_mod.ap.integ_rk(nX, X, t, cep.dt, Istim, Ksac);

            //  Electromechanics excitation-activation
            if (cem.aStress) {
              double epsX;
              cep_mod.ap.actv_strs(X(0), cep.dt, yl, epsX);
            }
          }
        } break; 

        case TimeIntegratioType::CN2: {
          for (int i = 0; i < nt; i++) {
            double t = t1 + static_cast<double>(i) * cep.dt;
            double Istim;
            if (t >= Ts-eps &&  t <= Te+eps) {
              Istim = cep.Istim.A;
            } else {
              Istim = 0.0;
            }

            cep_mod.ap.integ_cn2(nX, X, t, cep.dt, Istim, Ksac, IPAR, RPAR);

            // Electromechanics excitation-activation
            if (cem.aStress) {
              double epsX;
              cep_mod.ap.actv_strs(X(0), cep.dt, yl, epsX);
            }
          }
        } break; 
      } 
    } break; 

    case ElectrophysiologyModelType::BO: {
      Vector<int> IPAR(2); 
      Vector<double> RPAR(5);
      IPAR(0) = cep.odes.maxItr;
      IPAR(1) = 0.0;
      RPAR(0) = cep.odes.absTol;
      RPAR(1) = cep.odes.relTol;

      switch (cep.odes.tIntType) {
        case TimeIntegratioType::FE: {
          for (int i = 0; i < nt; i++) {
            double t = t1 + static_cast<double>(i) * cep.dt;
            double Istim;
            if (t >= Ts-eps &&  t <= Te+eps) {
              Istim = cep.Istim.A;
            } else {
              Istim = 0.0;
            }

            cep_mod.bo.integ_fe(cep.imyo, nX, X, t, cep.dt, Istim, Ksac, RPAR);

            // Electromechanics excitation-activation
            if (cem.aStress) {
              double epsX;
              cep_mod.bo.actv_strs(X(0), cep.dt, yl, epsX);
            } else if (cem.aStrain) {
              cep_mod.bo.actv_strn(X(3), lambda_new*lambda_new, cep.dt, yl);  // actv_strn expects I4f = lambda^2
            }
          }
        } break;

        case TimeIntegratioType::RK4: {
          for (int i = 0; i < nt; i++) {
            double t = t1 + static_cast<double>(i) * cep.dt;
            double Istim;
            if (t >= Ts-eps &&  t <= Te+eps) {
              Istim = cep.Istim.A;
            } else {
              Istim = 0.0;
            }

            cep_mod.bo.integ_rk(cep.imyo, nX, X, t, cep.dt, Istim, Ksac, RPAR);

            // Electromechanics excitation-activation
            if (cem.aStress) {
              double epsX;
              cep_mod.bo.actv_strs(X(0), cep.dt, yl, epsX);
            } else if (cem.aStrain) {
              cep_mod.bo.actv_strn(X(3), lambda_new*lambda_new, cep.dt, yl);  // actv_strn expects I4f = lambda^2
            }
          }
        } break;

        case TimeIntegratioType::CN2: {
          for (int i = 0; i < nt; i++) {
            double t = t1 + static_cast<double>(i) * cep.dt;
            double Istim;
            if (t >= Ts-eps &&  t <= Te+eps) {
              Istim = cep.Istim.A;
            } else {
              Istim = 0.0;
            }

            cep_mod.bo.integ_cn2(cep.imyo, nX, X, t, cep.dt, Istim, Ksac, IPAR, RPAR);

            // Electromechanics excitation-activation
            if (cem.aStress) {
              double epsX;
              cep_mod.bo.actv_strs(X(0), cep.dt, yl, epsX);
            } else if (cem.aStrain) {
              cep_mod.bo.actv_strn(X(3), lambda_new*lambda_new, cep.dt, yl);  // actv_strn expects I4f = lambda^2
            }
          }
        } break;
      } 
    } break; 

    case ElectrophysiologyModelType::FN: {
      Vector<int> IPAR(2); 
      Vector<double> RPAR(2);
      IPAR(0) = cep.odes.maxItr;
      IPAR(0) = cep.odes.maxItr;
      IPAR(1) = 0;
      RPAR(0) = cep.odes.absTol;
      RPAR(1) = cep.odes.relTol;

      switch (cep.odes.tIntType) {
        case TimeIntegratioType::FE: {
          for (int i = 0; i < nt; i++) {
            double t = t1 + static_cast<double>(i) * cep.dt;
            double Istim;
            if (t >= Ts-eps &&  t <= Te+eps) {
              Istim = cep.Istim.A;
            } else {
              Istim = 0.0;
            }
            cep_mod.fn.integ_fe(nX, X, t, cep.dt, Istim);
          }
        } break;

        case TimeIntegratioType::RK4: {
          for (int i = 0; i < nt; i++) {
            double t = t1 + static_cast<double>(i) * cep.dt;
            double Istim;
            if (t >= Ts-eps &&  t <= Te+eps) {
              Istim = cep.Istim.A;
            } else {
              Istim = 0.0;
            }
            cep_mod.fn.integ_rk(nX, X, t, cep.dt, Istim);
          }
        } break;

        case TimeIntegratioType::CN2: {
          for (int i = 0; i < nt; i++) {
            double t = t1 + static_cast<double>(i) * cep.dt;
            double Istim;
            if (t >= Ts-eps &&  t <= Te+eps) {
              Istim = cep.Istim.A;
            } else {
              Istim = 0.0;
            }
            cep_mod.fn.integ_cn2(nX, X, t, cep.dt, Istim, IPAR, RPAR);
           }
        } break;
      }
    } break; 

    case ElectrophysiologyModelType::TTP: {
      Vector<int> IPAR(2); 
      Vector<double> RPAR(18);
      IPAR(0) = cep.odes.maxItr;
      IPAR(1) = 0;
      RPAR(0) = cep.odes.absTol;
      RPAR(1) = cep.odes.relTol;
      double Ta_kPa = 0.0;
      double Ka_kPa = 0.0;
      switch (cep.odes.tIntType) {
        case TimeIntegratioType::FE: {
          for (int i = 0; i < nt; i++) {
            double t = t1 + static_cast<double>(i) * cep.dt;
            double Istim;
            if (t >= Ts-eps &&  t <= Te+eps) {
              Istim = cep.Istim.A;
            } else {
              Istim = 0.0;
            }
            cep_mod.ttp.integ_fe(cep.imyo, nX, nG, X, Xg, t, cep.dt, Istim, Ksac, RPAR);

            // Electromechanics excitation-activation
            if (cem.aStress) {
              double epsX;
              cep_mod.ttp.actv_strs(X(3), cep.dt, yl, epsX);
            } else if (cem.aStrain) {
              cep_mod.ttp.actv_strn(X(3), lambda_new*lambda_new, cep.dt, yl);  // actv_strn expects I4f = lambda^2
            }
          }
        } break;

        case TimeIntegratioType::RK4: {
          for (int i = 0; i < nt; i++) {
            double t = t1 + static_cast<double>(i) * cep.dt;
            double Istim;
            if (t >= Ts-eps &&  t <= Te+eps) {
              Istim = cep.Istim.A;
            } else {
              Istim = 0.0;
            }

            cep_mod.ttp.integ_rk(cep.imyo, nX, nG, X, Xg, t, cep.dt, Istim, Ksac, RPAR);

            // Electromechanics excitation-activation
            if (cem.aStress) {
              double epsX;
              cep_mod.ttp.actv_strs(X(3), cep.dt, yl, epsX);
            } else if (cem.aStrain) {
              cep_mod.ttp.actv_strn(X(3), lambda_new*lambda_new, cep.dt, yl);  // actv_strn expects I4f = lambda^2
            }
          }
        } break;

        case TimeIntegratioType::CN2: {
          for (int i = 0; i < nt; i++) {
            double t = t1 + static_cast<double>(i) * cep.dt;
            double Istim;
            if (t >= Ts-eps &&  t <= Te+eps) {
              Istim = cep.Istim.A;
            } else {
              Istim = 0.0;
            }

            cep_mod.ttp.integ_cn2(cep.imyo, nX, nG, X, Xg, t, cep.dt, Istim, Ksac, IPAR, RPAR);

            // Electromechanics excitation-activation
            if (cem.aStress) {
              double epsX;
              cep_mod.ttp.actv_strs(X(3), cep.dt, yl, epsX);
            } else if (cem.aStrain) {
              cep_mod.ttp.actv_strn(X(3), lambda_new*lambda_new, cep.dt, yl);  // actv_strn expects I4f = lambda^2
            }
          }
        } break;
      }

    } break; 
  } 

  if (isnan(X(0)) ||  isnan(yl)) {
    throw std::runtime_error("[cep_integ_l] A NaN has been computed during time integration of electrophysiology variables.");
  }
}

//-------------
// cep_integ_active_tension
//-------------
// Update active tension by integrating Land-Niederer ODE system during each Newton iteration
// This is called every Newton iteration before structural equations are solved
//
void cep_integ_active_tension(Simulation* simulation, const int iEq, const Array<double>& Dg, const Array<double>& Yg)
{
  using namespace consts;
  
  auto& com_mod = simulation->com_mod;
  auto& cep_mod = simulation->cep_mod;
  auto& cem = cep_mod.cem;
  auto& eq = com_mod.eq[iEq];
  
  if (!cem.cpld || !cem.aStress) {
    return;  // No active stress coupling, nothing to do
  }
  
  const int tnNo = com_mod.tnNo;
  const double dt = com_mod.dt;
  
  // Ensure arrays are allocated
  if (cem.Ta.size() != tnNo) {
    cem.Ta.resize(tnNo);
    cem.Ta = 0.0;
  }
  if (cem.Ka.size() != tnNo) {
    cem.Ka.resize(tnNo);
    cem.Ka = 0.0;
  }
  if (cem.Ya.size() != tnNo) {
    cem.Ya.resize(tnNo);
    cem.Ya = 0.0;
  }
  if (cem.lambda_prev.size() != tnNo) {
    cem.lambda_prev.resize(tnNo);
    cem.lambda_prev = 1.0;
  }
  if (cem.lambda_curr.size() != tnNo) {
    cem.lambda_curr.resize(tnNo);
    cem.lambda_curr = 1.0;
  }

  // Find structural equation to compute fiber stretch
  int structEq = -1;
  for (int a = 0; a < com_mod.nEq; a++) { 
    if (com_mod.eq[a].phys == EquationType::phys_struct || 
        com_mod.eq[a].phys == EquationType::phys_ustruct) {
      structEq = a;
      break;
    }
  }

  // Compute current lambda from displacement Dg
  Vector<double> lambda_curr_vec(tnNo);
  Vector<double> dlambda_dt_vec(tnNo);
  lambda_curr_vec = 1.0;
  dlambda_dt_vec = 0.0;

  if (structEq != -1) {
    for (int iM = 0; iM < com_mod.nMsh; iM++) {
      auto& msh = com_mod.msh[iM];
      if (msh.nFn != 0) {
        Vector<double> sA_new(msh.nNo);
        Vector<double> sDlambda_dt(msh.nNo);
        post::fib_strech(simulation, structEq, msh, Dg, Yg, sA_new, sDlambda_dt);
        for (int a = 0; a < msh.nNo; a++) {
          int Ac = msh.gN(a);
          lambda_curr_vec(Ac) = sA_new(a);
          dlambda_dt_vec(Ac) = sDlambda_dt(a);
        }
      }
    }
  }

  // Update lambda history
  cem.lambda_prev = cem.lambda_curr;
  cem.lambda_curr = lambda_curr_vec;

  // Find CEP equation to get calcium concentration
  int cepEq = -1;
  for (int a = 0; a < com_mod.nEq; a++) {
    if (com_mod.eq[a].phys == EquationType::phys_CEP) {
      cepEq = a;
      break;
    }
  }

  if (cepEq == -1) {
    std::cout << "[cep_integ_active_tension] WARNING: No CEP equation found, cannot get calcium concentration" << std::endl;
    return;
  }

  // Get domain information for CEP equation
  auto& cep_eq = com_mod.eq[cepEq];
  
  // Integrate Land model for each node
  if (com_mod.dmnId.size() != 0) {
    // Domain-based case
    for (int Ac = 0; Ac < tnNo; Ac++) {
      if (!all_fun::is_domain(com_mod, cep_eq, Ac, EquationType::phys_CEP)) {
        continue;
      }

      // Find domain for this node
      for (int iDmn = 0; iDmn < cep_eq.nDmn; iDmn++) {
        auto& dmn = cep_eq.dmn[iDmn];
        auto cPhys = dmn.phys;
        int dID = dmn.Id;

        if (cPhys != EquationType::phys_CEP || !utils::btest(com_mod.dmnId(Ac), dID)) {
          continue;
        }

        // Check if TTP model is used
        if (dmn.cep.cepType != ElectrophysiologyModelType::TTP) {
          continue;
        }

        // Get calcium concentration (index 3 in TTP model)
        double c_Ca = cep_mod.Xion(3, Ac);

        // Get Land model state vector for this node
        // Always start from states at time n, not from previous Newton iteration
        Vector<double> Y_land_node;
        if (cep_mod.Y_land_n.size() > 0) {
          Y_land_node = cep_mod.Y_land_n.col(Ac);  // Start from time n
        } else {
          Y_land_node = cep_mod.Y_land.col(Ac);  // Fallback if not saved
        }

        // print out Y_land_node(0) for debugging for Ac = 0
        // if (Ac == 0) {
        //   std::cout << "[cep_integ_active_tension] Y_land_node(0) = " << Y_land_node(0) << std::endl;
        //   std::cout << "[cep_integ_active_tension] Y_land_node(1) = " << Y_land_node(1) << std::endl;
        // }

        // Get lambda values
        double lambda_old = cem.lambda_prev(Ac);
        double lambda_new = cem.lambda_curr(Ac);
        double dlambda_dt = dlambda_dt_vec(Ac);

        // Integrate Land model
        double Ta_kPa = 0.0;
        double Ka_kPa = 0.0;
        cep_mod.ttp.actv_strs_land(Y_land_node, c_Ca, lambda_old, lambda_new, dlambda_dt, dt, Ta_kPa, Ka_kPa);

        // Store updated Land model states back to global array using set_col
        cep_mod.Y_land.set_col(Ac, Y_land_node);

        // Convert to dyne/cm² and store
        cem.Ta(Ac) = Ta_kPa *1.0e-3;
        cem.Ka(Ac) = Ka_kPa *1.0e-3;
        cem.Ya(Ac) = cem.Ta(Ac);
      }
    }
  } else {
    // Non-domain case
    for (int Ac = 0; Ac < tnNo; Ac++) {
      if (!all_fun::is_domain(com_mod, cep_eq, Ac, EquationType::phys_CEP)) {
        continue;
      }

      // Check if TTP model is used
      if (cep_eq.dmn[0].cep.cepType != ElectrophysiologyModelType::TTP) {
        continue;
      }

      // Get calcium concentration (index 3 in TTP model)
      double c_Ca = cep_mod.Xion(3, Ac);

      // Get Land model state vector for this node
      Vector<double> Y_land_node;
      if (cep_mod.Y_land_n.size() > 0) {
        Y_land_node = cep_mod.Y_land_n.col(Ac);
      } else {
        Y_land_node = cep_mod.Y_land.col(Ac);
      }

      if (Y_land_node.size() == 0) {
        continue;  // Skip if Land states not initialized
      }

      if (Ac == 0) {
          std::cout << "[cep_integ_active_tension] Y_land_node(0) = " << Y_land_node(0) << std::endl;
        }

      // Get lambda values
      double lambda_old = cem.lambda_prev(Ac);
      double lambda_new = cem.lambda_curr(Ac);
      double dlambda_dt = dlambda_dt_vec(Ac);

      // Integrate Land model
      double Ta_kPa = 0.0;
      double Ka_kPa = 0.0;
      cep_mod.ttp.actv_strs_land(Y_land_node, c_Ca, lambda_old, lambda_new, dlambda_dt, dt, Ta_kPa, Ka_kPa);

      // Store updated Land model states back to global array using set_col
      cep_mod.Y_land.set_col(Ac, Y_land_node);

      // Convert to dyne/cm² and store
      cem.Ta(Ac) = Ta_kPa *1.0e-3;
      cem.Ka(Ac) = Ka_kPa *1.0e-3;
      cem.Ya(Ac) = cem.Ta(Ac);
    }
  }
}
};


