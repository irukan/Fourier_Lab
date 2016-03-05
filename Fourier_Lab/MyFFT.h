//
//  MyFFT.h
//  Fourier_Lab
//
//  Created by kayama on 2016/03/01.
//  Copyright © 2016年 kayama. All rights reserved.
//

#ifndef MyFFT_h
#define MyFFT_h
#include <complex>

typedef std::complex<double> complex_t;

inline
void butterfly(int n, int q, complex_t* x) // バタフライ演算
// n : 系列長
// q : ブロックの開始位置
// x : バタフライ演算する系列(入出力)
{
    const int m = n/2;
    const double theta0 = 2*M_PI/n;
    
    if (n > 1) {
        for (int p = 0; p < m; p++) {
            const complex_t wp = complex_t(cos(p*theta0), -sin(p*theta0));
            const complex_t a = x[q + p + 0];
            const complex_t b = x[q + p + m];
            x[q + p + 0] =  a + b;
            x[q + p + m] = (a - b) * wp;
        }
        butterfly(n/2, q + 0, x);
        butterfly(n/2, q + m, x);
    }
}

inline
void bit_reverse(int n, complex_t* x) // ビット反転並べ替え
// n : 系列長
// x : ビット反転並べ替えする系列(入出力)
{
    for (int i = 0, j = 1; j < n-1; j++) {
        for (int k = n >> 1; k > (i ^= k); k >>= 1);
        if (i < j) std::swap(x[i], x[j]); // x[i]とx[j]を交換する
    }
}
inline
void fft(int n, complex_t* x) // フーリエ変換
// n : 系列長
// x : フーリエ変換する系列(入出力)
{
    butterfly(n, 0, x);
    bit_reverse(n, x);
    //for (int k = 0; k < n; k++) x[k] /= n;
}


inline
void FFT(const vector<double>& rSrc, const vector<double>&iSrc,
               vector<double>& rDest, vector<double>& iDest, vector<double>& spec)
{
    int dataN = rSrc.size();
    
    complex<double>* Src = new complex<double>[dataN];
    for (size_t i = 0; i < dataN; i++)
        Src[i] = complex<double>(rSrc[i], iSrc[i]);
    
    fft(dataN, Src);
    rDest.resize(dataN);
    iDest.resize(dataN);
    spec.resize(dataN);
    
    for (size_t i = 0 ;i < dataN; i++)
    {
        rDest[i] = Src[i].real();
        iDest[i] = Src[i].imag();
        spec[i] = sqrt(rDest[i]*rDest[i] + iDest[i]*iDest[i]);
    }
    delete [] Src;
}








inline
void FFTCore(const vector<double>& rSrc, const vector<double>&iSrc,
             vector<double>& rDest, vector<double>& iDest, vector<double>& spec, size_t start, size_t end);


inline
void NormalFFT(const vector<double>& rSrc, const vector<double>&iSrc,
                   vector<double>& rDest, vector<double>& iDest, vector<double>& spec)
{
    size_t dataN = rSrc.size();
    
    rDest.resize(dataN);
    iDest.resize(dataN);
    spec.resize(dataN);

    FFTCore(rSrc, iSrc, rDest, iDest, spec,  0, dataN);
    for (size_t i = 0; i < dataN; i++)
    {
        rDest[i] = rDest[dataN - i];
        iDest[i] = iDest[dataN - i];
        spec[i] = spec[dataN - i];
    }

}

inline
void FFTCore(const vector<double>& rSrc, const vector<double>&iSrc,
             vector<double>& rDest, vector<double>& iDest, vector<double>& spec, size_t start, size_t end)
{
    size_t dataN = (end - start);
    size_t dataN_1_4 = dataN / 4;
    size_t dataN_1_2 = dataN / 2;
    size_t dataN_3_4 = dataN_1_2 + dataN_1_4;
    
    double dx = M_2PI / dataN;
    
    
    for (size_t i = 0; i < dataN_1_2; i++)
    {
        double ReSum = 0.0;
        double ImSum = 0.0;
        
        for (size_t k = 0; k< dataN_1_2; k++)
        {
            double theta = dx * i * k;
            
            ReSum += rSrc[k] * cos(theta) + iSrc[k] * sin(theta);
            ImSum += -rSrc[k] * sin(theta) + iSrc[k] * cos(theta);
        }
        rDest[i] = ReSum;
        iDest[i] = ImSum;
        spec[i] = sqrt(ReSum * ReSum + ImSum * ImSum);
    }
    for (size_t i = dataN_1_2; i < end; i++)
    {
        double ReSum = 0.0;
        double ImSum = 0.0;
        
        for (size_t k = dataN_1_2; k< end; k++)
        {
            double theta = dx * i * k;
            
            ReSum += rSrc[k] * cos(theta) + iSrc[k] * sin(theta);
            ImSum += -rSrc[k] * sin(theta) + iSrc[k] * cos(theta);
        }
        rDest[i] = ReSum;
        iDest[i] = ImSum;
        spec[i] = sqrt(ReSum * ReSum + ImSum * ImSum);
    }
//    for (size_t i = dataN_1_2; i < dataN_3_4; i++)
//    {
//        double ReSum = 0.0;
//        double ImSum = 0.0;
//        
//        for (size_t k = dataN_1_2; k< dataN_3_4; k++)
//        {
//            double theta = dx * i * k;
//            
//            ReSum += rSrc[k] * cos(theta) + iSrc[k] * sin(theta);
//            ImSum += -rSrc[k] * sin(theta) + iSrc[k] * cos(theta);
//        }
//        rDest[i] = ReSum;
//        iDest[i] = ImSum;
//        spec[i] = sqrt(ReSum * ReSum + ImSum * ImSum);
//    }
//    for (size_t i = dataN_3_4; i < end; i++)
//    {
//        double ReSum = 0.0;
//        double ImSum = 0.0;
//        
//        for (size_t k = dataN_3_4; k< end; k++)
//        {
//            double theta = dx * i * k;
//            
//            ReSum += rSrc[k] * cos(theta) + iSrc[k] * sin(theta);
//            ImSum += -rSrc[k] * sin(theta) + iSrc[k] * cos(theta);
//        }
//        rDest[i] = ReSum;
//        iDest[i] = ImSum;
//        spec[i] = sqrt(ReSum * ReSum + ImSum * ImSum);
//    }
    

}


//inline
//void NormalFFT(const vector<double>& rSrc, const vector<double>&iSrc,
//               vector<double>& rDest, vector<double>& iDest, vector<double>& spec)
//{
//    int dataN = (int)rSrc.size();
//    int N_2 = dataN / 2;
//    double dx = M_2PI / dataN;
//    
//    // 2のべき常数
//    int ln2 = log(dataN);
//    
//    // データのコピー
//    double *real = new double[dataN]; // rSrc
//    double *imag = new double[dataN]; // iSrc
//    std::copy(rSrc.begin(), rSrc.end(), real);
//    std::copy(iSrc.begin(), iSrc.end(), imag);
//    
//    // 回転因子の計算
//    double *w = new double[dataN];
//    for (size_t i = 0; i < N_2; i++)
//    {
//        w[i] = cos(dx * i);
//        w[i + N_2] = -sin(dx * i);
//    }
//    
//    // ビット逆順への入力データ並び替え
//    double t1, t2;
//    int j = 0;
//    int k = 0;
//    for (size_t i = 0; i < N_2; i++)
//    {
//        if (i < j)
//        {
//            t1 = real[j];
//            real[j] = real[i];
//            real[i] = t1;
//            
//            t2 = imag[j];
//            imag[j] = imag[i];
//            imag[i] = t2;
//        }
//        k = N_2;
//        while (k <= j)
//        {
//            j = j - k;
//            k = k / 2;
//        }
//        
//        j = j + k;
//    }
//    
//    // バタフライ演算
//    for (int i = 1; i <= ln2; i++)
//    {
//        int m = pow(2, i);
//        int h = m / 2;
//        
//        for (int j = 0; j < h; j++)
//        {
//            double w1 = w[j * (dataN/m)]; // sin
//            double w2 = w[j * (dataN/m) + N_2]; // cos
//            
//            for (int k = 0; k< dataN; k += m)
//            {
//                int kp = k + h;
//                //cout << "m:" << m << " kp:" << kp << endl;
//                
//                double v1 = real[kp] * w1 - imag[kp] * w2;
//                double v2 = real[kp] * w2 + imag[kp] * w1;
//                
//                // ビットの戻し
//                t1 = real[k] + v1;
//                real[kp] = real[k] - v1;
//                real[k] = t1;
//                
//                t2 = imag[k] + v2;
//                imag[kp] = imag[k] - v2;
//                imag[k] = t2;
//            }
//        }
//    }
//    
//    // スペクトル
//    spec.resize(dataN);
//    for (size_t i = 0; i < dataN; i++)
//        spec[i] = sqrt(real[i] * real[i] + imag[i] * imag[i]);
//    
//    // Real, Imag
//    rDest = std::vector<double>(real, &real[dataN]);
//    iDest = std::vector<double>(imag, &imag[dataN]);
//    
//    delete [] w;
//}
#endif /* MyFFT_h */
