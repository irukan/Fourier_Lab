//
//  MyFunc.h
//  Fourier_Lab
//
//  Created by kayama on 2016/02/14.
//  Copyright © 2016年 kayama. All rights reserved.
//

#ifndef MyFunc_h
#define MyFunc_h
#include <vector>
#include <math.h>

using namespace std;

static double* cosTable;
static int cosTableNum;
void initTable(int num)
{
    cosTable = new double[num];
    cosTableNum = num;
    double dx = M_PI * 2 / num;
    for (int i=0; i< num; i++)
        cosTable[i] = cos(i * dx);
}

inline double myCos(double x)
{
    double sign=1;
    int i;
    
    if (x<0) x=-x;              /*  cos(x)=cos(-x)     */
    x/=M_PI;
    x=(x-(int)(x/2)*2);         /*  cos(x)=cos(x+2*PI) */
    if (x>1) x=2-x;             /*  cos(x)=cos(x-PI)   */
    if (x>.5) {x=1-x; sign=-1;} /* -cos(x)=cos(PI-x)   */
    i=(int)(x*2*cosTableNum);
    return (sign*cosTable[i]);
}
#define mySin(x) (myCos(x-M_PI/2))


#define coef2 1.0/2
#define coef4 1.0/24
#define coef6 1.0/720
#define coef8 1.0/40320
#define coef10 1.0/3628800
#define coef12 1.0/479001600

#define M_2PI M_PI*2.0
#define M_1_2PI 1/M_PI/2

inline double
myTaylorCos(double sita)
{
    int sign = 1;
    sita = sita - (int)(sita * M_1_2PI)*M_2PI;
    if (sita > M_PI)
        sita = M_2PI- sita;
    if (sita> M_PI_2)
    {
        sita = M_PI - sita;
        sign = -1;
    }
    double sitaMul2 = sita * sita;
    double keisu = sitaMul2;
    double ret = 1.0;
    ret -= coef2 * keisu;
    keisu *= sitaMul2;
    ret += coef4 * keisu;
    keisu *= sitaMul2;
    ret -= coef6 * keisu;
    //ret += coef8  * sita * sita * sita * sita * sita * sita * sita * sita;
    //ret -= coef10 * sita * sita * sita * sita * sita * sita * sita * sita * sita * sita;
    //ret += coef12 * sita * sita * sita * sita * sita * sita * sita * sita * sita * sita * sita * sita;
    
    return ret*sign;
}

#define myTaylorSin(sita) (myTaylorCos(sita-M_PI_2))
//#define coef3 1.0/6
//#define coef5 1.0/120
//#define coef7 1.0/5040
//#define coef9 1.0/362880
//
//inline double
//myTaylorSin(double sita)
//{
//    int sign = 1;
//    sita = sita - (int)(sita * M_1_2PI)*M_2PI;
//    if (sita > M_PI)
//    {
//        sita = M_2PI- sita;
//        sign = -1;
//    }
//    if (sita> M_PI_2)
//    {
//        sita = M_PI - sita;
//    }
//    double sitaMul2 = sita * sita;
//
//    double ret = sita;
//    double keisu = sita * sitaMul2;
//    ret -= coef3 * keisu;
//    keisu *= sitaMul2;
//    ret += coef5 * keisu;
//    keisu *= sitaMul2;
//    ret -= coef7 * keisu;
//    
//    return ret*sign;
//}


void MyDFT1(const vector<double>& rSrc, const vector<double>&iSrc,
               vector<double>& rDest, vector<double>& iDest, vector<double>& spec)
{
    size_t dataN = rSrc.size();
    double dx = 2 * M_PI / dataN;
    
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
            
            ReSum += rSrc[k] * myCos(theta) + iSrc[k] * mySin(theta);
            ImSum += -rSrc[k] * mySin(theta) + iSrc[k] * myCos(theta);
        }
        rDest[i] = ReSum;
        iDest[i] = ImSum;
        spec[i] = sqrt(ReSum * ReSum + ImSum * ImSum);
    }
}

inline
void MyDFT2(const vector<double>& rSrc, const vector<double>&iSrc,
            vector<double>& rDest, vector<double>& iDest, vector<double>& spec)
{
    size_t dataN = rSrc.size();
    double dx = 2 * M_PI / dataN;
    
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
            double sin = myTaylorSin(theta);
            double cos = myTaylorCos(theta);
            
            ReSum += rSrc[k] * cos + iSrc[k] * sin;
            ImSum += -rSrc[k] * sin + iSrc[k] * cos;
        }
        rDest[i] = ReSum;
        iDest[i] = ImSum;
        spec[i] = sqrt(ReSum * ReSum + ImSum * ImSum);
    }
}

#endif /* MyFunc_h */
