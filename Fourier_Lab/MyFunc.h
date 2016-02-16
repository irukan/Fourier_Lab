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
    
    //if (x<0) x=-x;              /*  cos(x)=cos(-x)     */
    x/=M_PI;
    x=(x-(int)(x/2)*2);         /*  cos(x)=cos(x+2*PI) */
    if (x>1) x=2-x;             /*  cos(x)=cos(x-PI)   */
    if (x>.5) {x=1-x; sign=-1;} /* -cos(x)=cos(PI-x)   */
    i=(int)(x*2*cosTableNum);
    return (sign*cosTable[i]);
}
#define mySin(x) (myCos(x-M_PI/2))

#define coef3 1.0/6
#define coef5 1.0/120
#define coef7 1.0/5040
#define coef9 1.0/362880

inline double
myTaylorSin(double sita)
{
    double ret = sita;
    //double sitaMul = sita;
    double sitaMul2 = sita * sita;
    //ret += sita;
    
    //sitaMul *= sitaMul2;
    ret -= coef3 * sitaMul2;

    //sitaMul *= sitaMul2;
    ret += coef5 * sitaMul2;

    //sitaMul *= sitaMul2;
    ret -= coef7 * sitaMul2;

    //sitaMul *= sitaMul2;
    ret += coef9 * sitaMul2;
    
    return ret;
}
#define myTaylorCos(sita) (myTaylorSin(sita-M_PI/2))

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
            
            ReSum += rSrc[k] * myTaylorCos(theta) + iSrc[k] * myTaylorSin(theta);
            ImSum += -rSrc[k] * myTaylorSin(theta) + iSrc[k] * myTaylorCos(theta);
        }
        rDest[i] = ReSum;
        iDest[i] = ImSum;
        spec[i] = sqrt(ReSum * ReSum + ImSum * ImSum);
    }
}

#endif /* MyFunc_h */
