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

inline
void NormalFFT(const vector<double>& rSrc, const vector<double>&iSrc,
                vector<double>& rDest, vector<double>& iDest, vector<double>& spec)
{
    int dataN = (int)rSrc.size();
    int N_2 = dataN / 2;
    double dx = M_2PI / dataN;
    
    // 2のべき常数
    int ln2 = log(dataN);

    // データのコピー
    double *real = new double[dataN]; // rSrc
    double *imag = new double[dataN]; // iSrc
    std::copy(rSrc.begin(), rSrc.end(), real);
    std::copy(iSrc.begin(), iSrc.end(), imag);
    
    // 回転因子の計算
    double *w = new double[dataN];
    for (size_t i = 0; i < N_2; i++)
    {
        w[i] = cos(dx * i);
        w[i + N_2] = -sin(dx * i);
    }
    
    // ビット逆順への入力データ並び替え
    double t1, t2;
    int j = 0;
    int k = 0;
    for (size_t i = 0; i < N_2; i++)
    {
        if (i < j)
        {
            t1 = real[j];
            real[j] = real[i];
            real[i] = t1;
            
            t2 = imag[j];
            imag[j] = imag[i];
            imag[i] = t2;
        }
        k = N_2;
        while (k <= j)
        {
            j = j - k;
            k = k / 2;
        }
        
        j = j + k;
    }
    
    // バタフライ演算
    for (int i = 1; i <= ln2; i++)
    {
        int m = pow(2, i);
        int h = m / 2;
        
        for (int j = 0; j < h; j++)
        {
            double w1 = w[j * (dataN/m)]; // sin
            double w2 = w[j * (dataN/m) + N_2]; // cos
            
            for (int k = 0; k< dataN; k += m)
            {
                int kp = k + h;
                
                double v1 = real[kp] * w1 - imag[kp] * w2;
                double v2 = real[kp] * w2 + imag[kp] * w1;
                
                // ビットの戻し
                t1 = real[k] + v1;
                real[kp] = real[k] - v1;
                real[k] = t1;
                
                t2 = imag[k] + v2;
                imag[kp] = imag[k] - v2;
                imag[k] = t2;
            }
        }
    }

    // スペクトル
    spec.resize(dataN);
    for (size_t i = 0; i < dataN; i++)
        spec[i] = sqrt(real[i] * real[i] + imag[i] * imag[i]);
    
    // Real, Imag
    rDest = std::vector<double>(real, &real[dataN]);
    iDest = std::vector<double>(imag, &imag[dataN]);
}

#endif /* NormalFunc_h */
