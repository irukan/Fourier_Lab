//
//  MyTableFunc.h
//  Fourier_Lab
//
//  Created by kayama on 2016/02/26.
//  Copyright © 2016年 kayama. All rights reserved.
//

#ifndef MyTableFunc_h
#define MyTableFunc_h
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

inline double myTableCos(double x)
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
#define myTableSin(x) (myTableCos(x-M_PI/2))


void MyTableDFT(const vector<double>& rSrc, const vector<double>&iSrc,
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
            
            ReSum += rSrc[k] * myTableCos(theta) + iSrc[k] * myTableSin(theta);
            ImSum += -rSrc[k] * myTableSin(theta) + iSrc[k] * myTableCos(theta);
        }
        rDest[i] = ReSum;
        iDest[i] = ImSum;
        spec[i] = sqrt(ReSum * ReSum + ImSum * ImSum);
    }
}

#endif /* MyTableFunc_h */
