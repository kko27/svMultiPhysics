/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "celec_mech.h"
using namespace std; 
#include <algorithm> 
#include "mat_fun.h"
#include "utils.h"
#include <math.h>

celec_mech::celec_mech()
{
}

celec_mech::~celec_mech()
{
}

void celec_mech::land_model(const Vector<double>& Y, Vector<double>& dY, const double Cai, 
                    const double nominal_lambda, const double dlambda_dt, double& Ta, double& Ka){

    // Land 2016 Model Parameters (Local variables for now)
    double CaRef = 0.805;
    double eta_Tm = 5.0; 
    double k_uw = 0.182; 
    double k_ws = 0.012; 
    double Tref = 120; 
    double k_TRPN = 0.1;
    double eta_TRPN = 2.0; 
    double k_u = 1.0; 
    double TRPN50 = 0.35; 
    double rw = 0.5; 
    double rs = 0.25; 
    double gamma_s = 0.0085; 
    double gamma_w = 0.615; 
    double phi = 2.23; 
    double Aeff = 25.0; 
    double beta0 = 2.3; 
    double beta1 = -2.4; 

    // State Variables 
    double XB = Y(0); 
    double XW = std::max(0.0, Y(1)); 
    double CaTRPN = std::max(0.0, Y(2)); 
    double XS = std::max(0.0, Y(3)); 
    double ZETAW = Y(4);  // Swapped to match MATLAB
    double ZETAS = Y(5);  // Swapped to match MATLAB

    // Cross-bridging Model 
    double lambda = std::min(1.2, nominal_lambda); 
    double LFac = std::max(0.0, 1.0 + beta0 * (lambda + std::min(0.87,lambda) - 1.87)); 
    double XU   = (1 - XB) - XS - XW; 

    double k_b = k_u * pow(TRPN50, eta_Tm)/ (1 - rs - (1 - rs)*rw);
    dY(0)  = k_b * std::min(100.0, pow(CaTRPN, -(eta_Tm/2.0))) * XU  - k_u * pow(CaTRPN, (eta_Tm/2.0)) *  XB;

    double k_wu = (k_uw * (1.0/rw - 1.0) - k_ws);
    double gamma_wu = gamma_w * std::abs(ZETAW); 
    dY(1)  = k_uw*XU - k_wu*XW - k_ws*XW - gamma_wu*XW;

    double CaT50 = CaRef + beta1*std::min(0.2, lambda - 1.0); 
    dY(2) = k_TRPN * (pow((Cai/CaT50), eta_TRPN) * (1.0 - CaTRPN) - CaTRPN); 

    double k_su = k_ws * rw * (1.0/rs - 1.0); 
    double term1, term2;
    if (ZETAS > 0.0) {term1 = ZETAS;} else {term1 = 0.0;}
    if (ZETAS < -1.0) {term2 = -ZETAS - 1.0;} else {term2 = 0.0;}
    double gamma_su = gamma_s * std::max(term1, term2);
    dY(3) = k_ws*XW - k_su*XS - gamma_su*XS; 

    double Aw = Aeff * rs/((1.0 - rs) * rw + rs); 
    double cw = phi * k_uw * (1.0 - rs)*(1.0 - rw) / ((1.0 - rs)*rw);
    dY(4) = Aw * dlambda_dt - cw*ZETAW;  // Y(4) = ZETAW

    double As = Aw; 
    double cs = phi * k_ws * (1.0 - rs)*rw/rs; 
    dY(5) = As * dlambda_dt - cs*ZETAS;  // Y(5) = ZETAS 
    // Active and Total Tension
    Ta = LFac * (Tref/rs) * ((ZETAS+1.0) * XS + (ZETAW) * XW);
    Ka = LFac * (Tref/rs) * (As*XS + Aw*XW); 
}

void celec_mech::LandJacobian(Array<double>& J, double t, double Cai, double nominal_lambda, const Vector<double>& y)
{
    //-------------------------------------------------------------------------------
    // Parameters
    //-------------------------------------------------------------------------------
    // Land 2016 Model Parameters (Local variables for now)
    double CaRef = 0.805;
    double eta_Tm = 5.0; 
    double k_uw = 0.182; 
    double k_ws = 0.012; 
    double Tref = 120; 
    double k_TRPN = 0.1;
    double eta_TRPN = 2.0; 
    double k_u = 1.0; 
    double TRPN50 = 0.35; 
    double rw = 0.5; 
    double rs = 0.25; 
    double gamma_s = 0.0085; 
    double gamma_w = 0.615; 
    double phi = 2.23; 
    double beta1 = -2.4; 

    // State Variables 
    double XB = y(0); 
    double XW = std::max(0.0, y(1)); 
    double CaTRPN = std::max(0.0, y(2)); 
    double XS = std::max(0.0, y(3)); 
    double ZETAW = y(4);  // Swapped to match MATLAB
    double ZETAS = y(5);  // Swapped to match MATLAB
    double lambda = std::min(1.2, nominal_lambda); 
    double XU = (1 - XB) - XS - XW; 

    // Row 1: dXB/dt 
    double k_b = k_u * (std::pow(TRPN50, eta_Tm))/(1.0 - rs - (1.0 - rs)*rw); 
    double TRPN_power = std::pow(CaTRPN, -eta_Tm/2.0);
    double min_factor = std::min(100.0, TRPN_power);
    J(0,0) = -k_b * min_factor - k_u * std::pow(CaTRPN, eta_Tm/2.0);  // d/dXB
    J(0,1) = -k_b * min_factor;                            // d/dXW
    J(0,3) = -k_b * min_factor;                            // d/dXS    
    // d/dCaTRPN 
    if (TRPN_power <= 100) {
      J(0,2) = k_b * (-eta_Tm/2.0) * std::pow(CaTRPN, -eta_Tm/2.0 - 1.0) * XU - k_u * (eta_Tm/2.0) * std::pow(CaTRPN, eta_Tm/2.0 - 1.0) * XB;
    } else {
      J(0,2) = - k_u * (eta_Tm/2.0) * std::pow(CaTRPN, eta_Tm/2.0 - 1.0) * XB;
    }

    // Row 2: dXW/dt 
    double k_wu = (k_uw * (1.0/rw - 1.0) - k_ws);
    double gamma_wu = gamma_w * std::abs(ZETAW); 
    J(1,0) = -k_uw;
    J(1,1) = -k_uw - k_wu - k_ws - gamma_wu; 
    J(1,3) = -k_uw;
    double sign_val = (ZETAW > 0.0) ? 1.0 : (ZETAW < 0.0) ? -1.0 : 0.0;
    J(1,4) = -gamma_w * sign_val * XW;

    // Row 3: dCaTRPN/dt 
    double CaT50 = CaRef + beta1*std::min(0.2, lambda - 1.0); 
    J(2,2) = k_TRPN*(-1*std::pow((Cai/CaT50), eta_TRPN) - 1.0); 

    // Row 4: dXS/dt 
    double k_su = k_ws * rw * (1.0/rs - 1.0); 
    // Compute gamma_su using piecewise logic
    double gamma_su = 0.0;
    if (ZETAS > 0.0) {
        gamma_su = gamma_s * ZETAS;
    } else if (ZETAS < -1.0) {
        gamma_su = gamma_s * (-ZETAS - 1.0);
    }
    // Derivative of gamma_su with respect to ZETAS
    double d_gamma_su_rate = 0.0;
    if (ZETAS > 0.0) {
        d_gamma_su_rate = gamma_s;
    } else if (ZETAS < -1.0) {
        d_gamma_su_rate = -gamma_s;
    }
    J(3,1) = k_ws;
    J(3,3) = -k_su - gamma_su;
    J(3,5) = -d_gamma_su_rate * XS;

    // Row 5: dZETAW/dt (Y(4) = ZETAW)
    double cw = phi*k_uw*(1.0 - rs)*(1.0 - rw)/((1.0 - rs)*rw); 
    J(4,4) = -cw; 

    // Row 6: dZETAS/dt (Y(5) = ZETAS)
    double cs = phi*k_ws*(1.0 - rs)*rw/rs;
    J(5,5) = -cs; 
    }

void celec_mech::integ_rk(const int nY, Vector<double>& Y, double& Ta, double& Ka,
    const double dt, const double Cai, const double nominal_lambda, const double dlambda_dt)
{
  double dt6 = dt / 6.0;
  Vector<double> frk1(nY), frk2(nY), frk3(nY), frk4(nY);

  // RK4: 1st pass
  auto Yrk = Y;
  land_model(Yrk, frk1, Cai, nominal_lambda, dlambda_dt, Ta, Ka); 

  // RK4: 2nd pass
  Yrk = Y + 0.5*dt*frk1;
  land_model(Yrk, frk2, Cai, nominal_lambda, dlambda_dt, Ta, Ka); 

  // RK4: 3rd pass
  Yrk = Y + 0.5*dt*frk2;
  land_model(Yrk, frk3, Cai, nominal_lambda, dlambda_dt, Ta, Ka);

  // RK4: 4th pass
  Yrk = Y + dt*frk3;
  land_model(Yrk, frk4, Cai, nominal_lambda, dlambda_dt, Ta, Ka);

  Y = Y + dt6 * (frk1 + 2.0*(frk2 + frk3) + frk4);
}

/// @brief Time integration performed using Crank-Nicholson method
void celec_mech::integ_cn2(const int nY, Vector<double>& Yn, double& Ta, double& Ka, const double Ts, 
    const double dt, const double Cai, const double nominal_lambda, const double dlambda_dt, Vector<int>& IPAR, Vector<double>& RPAR)
{
  int itMax = IPAR(0);
  double atol  = RPAR(0);
  double rtol  = RPAR(1);
  auto Im = mat_fun::mat_id(nY);

  Vector<double> fn(nY);
  land_model(Yn, fn, Cai, nominal_lambda, dlambda_dt, Ta, Ka);

  int k  = 0;
  auto Yk = Yn;
  bool l1 = false;
  bool l2 = false;
  bool l3 = false;
  double t = Ts + dt;
  double eps = std::numeric_limits<double>::epsilon();

  while (true) {
    k = k + 1;
    Vector<double> fk(nY);
    land_model(Yk, fk, Cai, nominal_lambda, dlambda_dt, Ta, Ka);

    auto rK = Yk - Yn - 0.5*dt*(fk + fn);

    double rmsA = 0.0;
    double rmsR = 0.0;

    for (int i = 0; i < nY; i++) {
      rmsA = rmsA + pow(rK(i),2.0);
      rmsR = rmsR + pow(rK(i) / (Yk(i)+eps), 2.0);
    }

    rmsA = sqrt(rmsA / static_cast<double>(nY));
    rmsR = sqrt(rmsR / static_cast<double>(nY));

    l1 = (k > itMax);
    l2 = (rmsA <= atol);
    l3  = (rmsR <= rtol);
    if (l1 || l2 || l3) {
      break;
    }

    Array<double> JAC(nY, nY);
    LandJacobian(JAC, t, Cai, nominal_lambda, Yk);
    JAC = Im - 0.5*dt*JAC;
    JAC = mat_fun::mat_inv(JAC, nY);
    rK  = mat_fun::mat_mul(JAC, rK);
    Yk  = Yk - rK;
  }
  Yn = Yk;

  land_model(Yn, fn, Cai, nominal_lambda, dlambda_dt,Ta, Ka);

  if (!l2 && !l3) {
    IPAR(1) = IPAR(1) + 1;
  }
}