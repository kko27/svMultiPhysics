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

 void celec_mech::land_model(const Vector<double>& Y, Vector<double>& dY, const double Cai, 
                    const double lambda, const double dlambda_dt, double& T, double& Ta, double& Tp){

    // Land 2016 Model Parameters (Local variables for now)
    double perm50 = 0.35;
    double TRPN_n = 2;
    double koff = 0.1;
    double dr = 0.25;
    double wfrac = 0.5;
    double TOT_A = 25;
    double ktm_unblock = 1;
    double beta_1 = -2.4;
    double beta_0 = 2.3;
    double gamma = 0.0085;
    double gamma_wu = 0.615;
    double phi = 2.23;
    double nperm = 5;
    double ca50 = 0.805;
    double Tref = 120;
    double nu = 7;
    double mu = 3;

    double k_ws = 0.004 * mu;
    double k_uw = 0.026  * nu;
    double cdw = phi * k_uw * (1-dr)*(1-wfrac) /  ((1-dr)*wfrac);
    double cds = phi * k_ws * (1-dr)*wfrac / dr;
    double k_wu = k_uw * (1/wfrac - 1) - k_ws;
    double k_su = k_ws * (1/dr - 1) * wfrac; 
    double A = (0.25 * TOT_A) / ((1-dr)*wfrac + dr) * (dr/0.25);

    // State Variables 
    double XS = std::max(0.0, Y(1));
    double XW = std::max(0.0, Y(2));
    double TRPN = std::max(0.0, Y(3)); 
    double TmBlocked = Y(4);
    double ZETAS = Y(5);
    double ZETAW = Y(6); 

    // Cross-bridging Model 
    double lambda = std::min(1.2, lambda); 
    double LFac = std::max(0.0, 1 + beta_0 * (lambda + std::min(0.87,lambda) - 1.87)); 
    double XU   = (1 - TmBlocked) - XW - XS;
    double xb_ws = k_ws * XW;
    double xb_uw = k_uw * XU ;
    double xb_wu = k_wu * XW;
    double xb_su = k_su * XS;

    double term1, term2;
    if (ZETAS > 0.0) {term1 = ZETAS;} else {term1 = 0.0;}
    if (ZETAS < -1.0) {term2 = -ZETAS - 1.0;} else {term2 = 0.0;}
    double gamma_rate = gamma * std::max(term1, term2);

    double xb_su_gamma = gamma_rate * XS;
    double gamma_rate_w  = gamma_wu * std::abs(ZETAW);
    double xb_wu_gamma = gamma_rate_w * XW;  

    dY(1)  = xb_ws - xb_su - xb_su_gamma;
    dY(2)  = xb_uw - xb_wu - xb_ws - xb_wu_gamma;
    ca50 = ca50 + beta_1*min(0.2,lambda - 1);
    dY(3)  = koff * (pow((Cai/ca50), TRPN_n) * (1-TRPN) - TRPN); 

    double XSSS = dr * 0.5;
    double XWSS = (1-dr) * wfrac * 0.5;
    double ktm_block = ktm_unblock * pow(perm50, nperm) *  0.5 / (0.5 - XSSS - XWSS);
    dY(4)  = ktm_block * min(100.0, pow(TRPN, -(nperm/2))) * XU  - ktm_unblock * pow(TRPN, (nperm/2)) *  TmBlocked;

    // Velocity Dependence
    dY(5) = A * dlambda_dt - cds * ZETAS;
    dY(6) = A * dlambda_dt - cdw * ZETAW;

    // Passive Model
    double par_k = 7;
    double b = 9.1;
    double eta_l = 200;
    double eta_s = 20;
    double a = 2.1;

    double Cd = Y(7);
    double C = lambda - 1;
    double eta; 

    if (C - Cd < 0)
        {eta = eta_s;}
    else
        {eta = eta_l;}

    double dCd_dt = par_k * (C - Cd)/eta; 
    dY(7) = dCd_dt;

    double Fd = eta * dCd_dt;
    double F1 = (exp( b * C) - 1);
    Tp = a * (F1 + Fd);

    // Active and Total Tension
    Ta = LFac * (Tref/dr) * ((ZETAS+1) * XS + (ZETAW) * XW);
    T = Ta + Tp;
}

void celec_mech::integ_rk(const int imyo, const int nY, Vector<double>& Y, double& T,  double& Ta, double& Tp, 
    const double dt, const double Cai, const double lambda, const double dlambda_dt)
{
  double dt6 = dt / 6.0;
  Vector<double> frk1(nY), frk2(nY), frk3(nY), frk4(nY);

  // RK4: 1st pass
  auto Yrk = Y;
  land_model(Yrk, frk1, Cai, lambda, dlambda_dt, T, Ta, Tp); 

  // RK4: 2nd pass
  Yrk = Y + 0.5*dt*frk1;
  land_model(Yrk, frk2, Cai, lambda, dlambda_dt, T, Ta, Tp); 

  // RK4: 3rd pass
  Yrk = Y + 0.5*dt*frk2;
  land_model(Yrk, frk3, Cai, lambda, dlambda_dt, T, Ta, Tp);
  
  // RK4: 4th pass
  Yrk = Y + dt*frk3;
  land_model(Yrk, frk4, Cai, lambda, dlambda_dt, T, Ta, Tp);

  Y = Y + dt6 * (frk1 + 2.0*(frk2 + frk3) + frk4);
}
