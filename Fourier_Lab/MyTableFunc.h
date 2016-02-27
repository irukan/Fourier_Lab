//
//  MyTableFunc.h
//  Fourier_Lab
//
//  Created by kayama on 2016/02/26.
//  Copyright © 2016年 kayama. All rights reserved.
//

#ifndef MyTableFunc_h
#define MyTableFunc_h
#include "MyTaylorFunc.h"

using namespace std;

static double* cosTable;
static double* sinTable;
static int TableNum;
void initTable(int num)
{
    cosTable = new double[num];
    sinTable = new double[num];
    TableNum = num;
    
    double dx = M_2PI / num;
    
    for (int i = 0; i< num; i++)
    {
        sinTable[i] = sin(i * dx);
        cosTable[i] = cos(i * dx);
    }
}

inline
void convIdx(int* idx)
{
    *idx -= (*idx / TableNum) * TableNum;
}
inline
double myTableCos(size_t idx)
{
    return cosTable[idx];
}

inline
double myTableSin(size_t idx)
{
    return sinTable[idx];
}

inline
void MyTableDFT(const vector<double>& rSrc, const vector<double>&iSrc,
                vector<double>& rDest, vector<double>& iDest, vector<double>& spec)
{
    size_t dataN = rSrc.size();
    
    rDest.resize(dataN);
    iDest.resize(dataN);
    spec.resize(dataN);
    
    for (int i = 0; i < dataN; i++)
    {
        double ReSum = 0.0;
        double ImSum = 0.0;
        
        for (int k = 0; k< dataN; k++)
        {
            int index = i * k;
            convIdx(&index);
            
            double Sin = myTableSin(index);
            double Cos = myTableCos(index);
            
            ReSum += rSrc[k] * Cos + iSrc[k] * Sin;
            ImSum += -rSrc[k] * Sin + iSrc[k] * Cos;
        }
        rDest[i] = ReSum;
        iDest[i] = ImSum;
        spec[i] = sqrt(ReSum * ReSum + ImSum * ImSum);
    }
}

#endif /* MyTableFunc_h */
