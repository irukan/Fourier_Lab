//
//  MyTaylorFunc_h
//  Fourier_Lab
//
//  Created by kayama on 2016/02/14.
//  Copyright © 2016年 kayama. All rights reserved.
//

#ifndef MyTaylorFunc_h
#define MyTaylorFunc_h
#include <vector>
#include <math.h>

using namespace std;

#define M_2PI M_PI*2.0
#define M_1_2PI 1/M_PI/2

inline
void rotateTheta(double* theta, int* signSin, int* signCos)
{
    *theta = *theta - (int)(*theta * M_1_2PI)*M_2PI;
    *signSin = 1;
    *signCos = 1;
    if (*theta > M_PI)
    {
        *theta = M_2PI- *theta;
        *signSin = -1; // sin
    }
    
    if (*theta > M_PI_2)
    {
        *theta = M_PI - *theta;
        *signCos = -1; // cos
    }
}

#define coef2 1.0/2
#define coef4 1.0/24
#define coef6 1.0/720
#define coef8 1.0/40320
#define coef10 1.0/3628800
#define coef12 1.0/479001600

inline double
myTaylorCos(double theta, int sign)
{    
    double thetaMul2 = theta * theta;
    double keisu = thetaMul2;
    double ret = 1.0;
    ret -= coef2 * keisu;
    keisu *= thetaMul2;
    ret += coef4 * keisu;
    keisu *= thetaMul2;
    ret -= coef6 * keisu;
    keisu *= thetaMul2;
    ret += coef8 * keisu;
    //ret += coef8  * theta * theta * theta * theta * theta * theta * theta * theta;
    //ret -= coef10 * theta * theta * theta * theta * theta * theta * theta * theta * theta * theta;
    //ret += coef12 * theta * theta * theta * theta * theta * theta * theta * theta * theta * theta * theta * theta;
    
    return ret * sign;
}

///#define myTaylorSin(theta,sign) (myTaylorCos(theta-M_PI_2,sign))
#define coef3 1.0/6
#define coef5 1.0/120
#define coef7 1.0/5040
#define coef9 1.0/362880

inline double
myTaylorSin(double theta, int sign)
{
    double thetaMul2 = theta * theta;

    double ret = theta;
    double keisu = theta * thetaMul2;
    ret -= coef3 * keisu;
    keisu *= thetaMul2;
    ret += coef5 * keisu;
    keisu *= thetaMul2;
    ret -= coef7 * keisu;
    keisu *= thetaMul2;
    ret += coef9 * keisu;
    
    return ret * sign;
}

inline
void MyTaylorDFT(const vector<double>& rSrc, const vector<double>&iSrc,
            vector<double>& rDest, vector<double>& iDest, vector<double>& spec)
{
    size_t dataN = rSrc.size();
    double dx = M_2PI / dataN;
    
    rDest.resize(dataN);
    iDest.resize(dataN);
    spec.resize(dataN);
    
    for (size_t i = 0; i < dataN; i++)
    {
        double ReSum = 0.0;
        double ImSum = 0.0;
        
        for (size_t k = 0; k< dataN; k++)
        {
            double theta = dx * i * k;
            int signSin, signCos;
            rotateTheta(&theta, &signSin, &signCos);
            double sin = myTaylorSin(theta, signSin);
            double cos = myTaylorCos(theta, signCos);
            
            ReSum += rSrc[k] * cos + iSrc[k] * sin;
            ImSum += -rSrc[k] * sin + iSrc[k] * cos;
        }
        rDest[i] = ReSum;
        iDest[i] = ImSum;
        spec[i] = sqrt(ReSum * ReSum + ImSum * ImSum);
    }
}

// 回転因子計算をN^2ループ外で行う
inline
void MyTaylorDFT2(const vector<double>& rSrc, const vector<double>&iSrc,
                 vector<double>& rDest, vector<double>& iDest, vector<double>& spec)
{
    int dataN = (int)rSrc.size();
    double dx = M_2PI / dataN;
    
    rDest.resize(dataN);
    iDest.resize(dataN);
    spec.resize(dataN);
    
    // 回転因子の計算
    double *w1 = new double[dataN];
    double *w2 = new double[dataN];
    for (size_t i = 0; i < dataN; i++)
    {
        double theta = dx * i;
        int signSin, signCos;
        rotateTheta(&theta, &signSin, &signCos);

        w1[i] = myTaylorCos(theta, signCos);
        w2[i] = myTaylorSin(theta, signSin);
    }
    
    int index = -1;
    for (int i = 0; i < dataN; i++)
    {
        double ReSum = 0.0;
        double ImSum = 0.0;
        
        for (int k = 0; k< dataN; k++)
        {
            index = i * k;
            index -= (index/dataN) * dataN;
            
            ReSum += rSrc[k] * w1[index] + iSrc[k] * w2[index];
            ImSum += -rSrc[k] * w2[index] + iSrc[k] * w1[index];
        }
        rDest[i] = ReSum;
        iDest[i] = ImSum;
        spec[i] = sqrt(ReSum * ReSum + ImSum * ImSum);
    }
}

#endif /* MyTaylorFunc_h */
