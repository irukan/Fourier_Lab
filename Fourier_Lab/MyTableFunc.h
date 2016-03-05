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
#include <tmmintrin.h>
#include <xmmintrin.h>

using namespace std;

static double* cosTable;
static double* sinTable;
static double* sseTableSin;
static double* sseTableCos;
static int TableNum;
static int SSETableNum;

void initTable(int num)
{
    cosTable = new __attribute__((aligned(16))) double[num];
    sinTable = new __attribute__((aligned(16))) double[num];
    TableNum = num;
    
    sseTableSin = new __attribute__((aligned(16))) double[num * 2];
    sseTableCos = new __attribute__((aligned(16))) double[num * 2];
    
    SSETableNum = num * 2;
    double dx = M_2PI / num;
    
    for (int i = 0; i< num; i++)
    {
        sinTable[i] = sin(i * dx);
        cosTable[i] = cos(i * dx);
    }
    
    for (int i = 0; i< num; i++)
    {
        sseTableSin[i*2] = sin(i * dx);//Real
        sseTableSin[i*2 + 1] = -sin(i * dx);//Imag
        
        sseTableCos[i*2] = cos(i * dx);//Real
        sseTableCos[i*2 + 1] = cos(i * dx);//Imag
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
void convIdxSSE(int *idx)
{
    *idx -= (*idx / TableNum) * TableNum;
    *idx *= 2;
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
            
            ReSum +=  iSrc[k] * Sin + rSrc[k] * Cos;
            ImSum += -rSrc[k] * Sin + iSrc[k] * Cos;
        }
        rDest[i] = ReSum;
        iDest[i] = ImSum;
        spec[i] = sqrt(ReSum * ReSum + ImSum * ImSum);
    }
}

inline
void MyTableDFT_SSE(const vector<double>& rSrc, const vector<double>&iSrc,
                vector<double>& rDest, vector<double>& iDest, vector<double>& spec)
{
    size_t dataN = rSrc.size();
    
    rDest.resize(dataN);
    iDest.resize(dataN);
    spec.resize(dataN);
    
    for (int i = 0; i < dataN; i ++)
    {
        __m128d Re_sum = _mm_setzero_pd();
        __m128d Im_sum = _mm_setzero_pd();

        for (int k = 0; k< dataN; k += 2)
        {
            __m128d Re_src = _mm_loadu_pd(&rSrc[k]);
            __m128d Im_src_ = _mm_loadu_pd(&iSrc[k]);
            
            int index1 = i * k;
            convIdx(&index1);
            int index2 = i * (k+1);
            convIdx(&index2);
            
//            __m128d Sin = _mm_load_pd(&sinTable[index]);
//            __m128d Cos = _mm_load_pd(&cosTable[index]);
            __m128d Sin = _mm_set_pd(sinTable[index2], sinTable[index1]);
            __m128d Cos = _mm_set_pd(cosTable[index2], cosTable[index1]);
            
            __m128d temp_Re = _mm_add_pd(_mm_mul_pd(Re_src, Cos), _mm_mul_pd(Im_src_, Sin));
            __m128d temp_Im = _mm_sub_pd(_mm_mul_pd(Im_src_, Cos), _mm_mul_pd(Re_src, Sin));
            
            Re_sum = _mm_add_pd(Re_sum, temp_Re);
            Im_sum = _mm_add_pd(Im_sum, temp_Im);
        }
        
        rDest[i] = Re_sum[0] + Re_sum[1];
        iDest[i] = Im_sum[0] + Im_sum[1];
        spec[i] = sqrt(rDest[i] * rDest[i] + iDest[i] * iDest[i]);
//        _mm_storeu_pd(&rDest[i], Re_sum);
//        _mm_storeu_pd(&iDest[i], Im_sum);
//        _mm_storeu_pd(&spec[i], _mm_sqrt_pd(_mm_add_pd(_mm_mul_pd(Re_sum, Re_sum), _mm_mul_pd(Im_sum, Im_sum))));
    }
}

inline
void MyTableDFT_SSE2(const vector<double>& rSrc, const vector<double>&iSrc,
                    vector<double>& rDest, vector<double>& iDest, vector<double>& spec)
{
    size_t dataN = rSrc.size();
    
    rDest.resize(dataN);
    iDest.resize(dataN);
    spec.resize(dataN);
    
    double* complex = new __attribute__((aligned(16))) double[dataN * 2];
    for (int i = 0; i < dataN; i++)
    {
        complex[2*i] = rSrc[i];
        complex[2*i + 1] = iSrc[i];
    }
    
    __m128d minus = {1, -1};
    for (int i = 0; i < dataN; i++)
    {
        __m128d Sum = _mm_setzero_pd();
        
        for (int k = 0; k< dataN; k++)
        {
            int index = i * k;
            convIdx(&index);
            __m128d Cos = _mm_set1_pd(cosTable[index]);
            __m128d Sin = _mm_set1_pd(sinTable[index]);
            
            __m128d src = _mm_load_pd(&complex[k*2]);
            
            Sum = _mm_add_pd(Sum, _mm_mul_pd(src, Cos));
            Sum = _mm_add_pd(Sum, _mm_mul_pd(_mm_shuffle_pd(src, src, _MM_SHUFFLE2(0, 1)), _mm_mul_pd(Sin, minus)));
            //ReSum +=  rSrc[k] * Cos + iSrc[k] * Sin;
            //ImSum +=  iSrc[k] * Cos + rSrc[k] * -Sin;
        }
        rDest[i] = Sum[0];
        iDest[i] = Sum[1];
        Sum = _mm_mul_pd(Sum, Sum);
        spec[i] = Sum[0] + Sum[1];
    }
}


#endif /* MyTableFunc_h */
