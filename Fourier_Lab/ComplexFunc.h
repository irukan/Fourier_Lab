//
//  ComplexFunc.h
//  Fourier_Lab
//
//  Created by kayama on 2016/03/05.
//  Copyright © 2016年 kayama. All rights reserved.
//

#ifndef ComplexFunc_h
#define ComplexFunc_h
#include <complex>
using namespace std;

inline
void ComplexDFT(const vector<double>& rSrc, const vector<double>&iSrc,
                vector<double>& rDest, vector<double>& iDest, vector<double>& spec)
{
    size_t dataN = rSrc.size();
    double dx = M_2PI / dataN;
    
    rDest.resize(dataN);
    iDest.resize(dataN);
    spec.resize(dataN);
    
    complex<double>* Src = new complex<double>[dataN];
    for (size_t i = 0; i < dataN; i++)
        Src[i] = complex<double>(rSrc[i], iSrc[i]);

    for (size_t i = 0; i < dataN; i++)
    {
        complex<double> sum(0, 0);
        for (size_t j = 0; j < dataN; j++)
        {
            double theta = dx * i * j;
            sum += Src[j] * exp(complex<double>(0.0, theta));
        }
        rDest[i] = sum.real();
        iDest[i] = sum.imag();
        spec[i] = sqrt(rDest[i]*rDest[i] + iDest[i]*iDest[i]);
    }
    
    delete [] Src;
}
#endif /* ComplexFuncr_ */
