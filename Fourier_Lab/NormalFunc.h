//
//  NormalFunc.h
//  Fourier_Lab
//
//  Created by kayama on 2016/02/14.
//  Copyright © 2016年 kayama. All rights reserved.
//

#ifndef NormalFunc_h
#define NormalFunc_h
#include<vector>
#include<iostream>
#include<math.h>
#include "Plotter.h"
using namespace std;

inline
void NormalDFT(const vector<double>& rSrc, const vector<double>&iSrc,
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

            ReSum += rSrc[k] * cos(theta) + iSrc[k] * sin(theta);
            ImSum += -rSrc[k] * sin(theta) + iSrc[k] * cos(theta);
        }
        rDest[i] = ReSum;
        iDest[i] = ImSum;
        spec[i] = sqrt(ReSum * ReSum + ImSum * ImSum);
    }
}

// 回転因子計算をN^2ループ外で行う
inline
void NormalDFT2(const vector<double>& rSrc, const vector<double>&iSrc,
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
        w1[i] = cos(dx * i);
        w2[i] = sin(dx * i);
    }
    
    int index = -1;
    for (int i = 0; i < dataN; i++)
    {
        double ReSum = 0.0;
        double ImSum = 0.0;
        
        for (int k = 1; k< dataN; k++)
        {
            //index = (index + k) % dataN;
            index = i * k;
            index -= (index/dataN) * dataN;
            
            ReSum += rSrc[k] * w1[index] + iSrc[k] * w2[index];
            ImSum += -rSrc[k] * w2[index] + iSrc[k] * w1[index];
        }
        rDest[i] = ReSum;
        iDest[i] = ImSum;
        spec[i] = sqrt(ReSum * ReSum + ImSum * ImSum);
    }
    
    delete [] w1;
    delete [] w2;
}


#endif /* NormalFunc_h */
