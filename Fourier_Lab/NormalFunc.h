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
#include<math.h>
using namespace std;

void NormalDFT(const vector<double>& rSrc, const vector<double>&iSrc,
                vector<double>& rDest, vector<double>& iDest, vector<double>& spec)
{
    size_t dataN = rSrc.size();
    double dx = 2 * M_PI / dataN;
    
    rDest.resize(dataN);
    iDest.resize(dataN);
    spec.resize(dataN);
    
    for (size_t i =0; i < dataN; i++)
    {
        double ReSum = 0.0;
        double ImSum = 0.0;
        
        for (size_t k=0; k< dataN; k++)
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
#endif /* NormalFunc_h */
