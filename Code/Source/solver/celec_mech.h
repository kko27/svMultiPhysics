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

#ifndef CELEC_MECH_H
#define CELEC_MECH_H
#include "Vector.h"
#include "Array.h"

class celec_mech {
public:
    celec_mech(); 
    ~celec_mech();
    void land_model(const Vector<double>& Y, Vector<double>& dY, const double Cai, const double nominal_lambda, 
        const double dlambda_dt, double& Ta, double& Ka);

    void LandJacobian(Array<double>& J, double t, double Cai, 
            double lambda, const Vector<double>& y);

    void integ_rk(const int nY, Vector<double>& Y, double& Ta, double& Ka, 
    const double dt, const double Cai, const double nominal_lambda, const double dlambda_dt);

    void integ_cn2(const int nY, Vector<double>& Yn, double& Ta, double& Ka, const double Ts, 
    const double dt, const double Cai, const double nominal_lambda, const double dlambda_dt, Vector<int>& IPAR, Vector<double>& RPAR);
};

#endif