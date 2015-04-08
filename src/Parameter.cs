using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;

namespace IBD
{
    class Parameter
    {
        public double SITA = 0;
        public double K = 0;
        public double factorialK = 0;
        public double POWSK;
        public double h;
        public void creatPar(ArrayList controlP)
        {
            double s;
            double F1 = 1.0;
            double F2 = 1.0;
            double Sum = 0;
            double Lnsum = 0;
            for (int i = 0; i < controlP.Count; i++)
            {
                Sum = Sum + (double)controlP[i];
                Lnsum = Lnsum + Math.Log((double)controlP[i], Math.E);
            }
            s = Math.Log(Sum / (controlP.Count), Math.E) - Lnsum / controlP.Count;
            K = (3 - s + Math.Sqrt(Math.Pow((s - 3), 2) + 24 * s)) / (12 * s);
            SITA = Sum / (controlP.Count * K);
            h = K * SITA;
            POWSK = Math.Pow(SITA, K);
            int vv = (int)K + 4;
            int jj = (int)(K + 5);
            for (int i = jj - 1; i > 0; i--)
            {
                F1 = F1 * i;

            }
            for (int i = vv - 1; i > 0; i--)
            {
                F2 = F2 * i;
            }
            factorialK = ((F1 - F2) * (K + 4 - vv) + F2) / ((K + 3) * (K + 2) * (K + 1) * K);
        }
    }
}
