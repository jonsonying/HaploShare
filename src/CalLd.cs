using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
namespace IBD
{

    class CalLd
    {
        public double par;
        public double[,] Dx;
        public double r2;
       
        public void CalComGroup(Inputdata x)
        {
            Dx = new double[x.controlData.GetUpperBound(1), x.controlData.GetUpperBound(1)];
            double[,] linkageData = new double[x.Tem - 6, 6];
            for (int i = 0; i < x.controlData.GetUpperBound(1); i++)
            {
                for (int j = i + 1; j < x.controlData.GetUpperBound(1); j++)
                {
                    int AABB = 0;
                    int AABb = 0;
                    int AaBB = 0;
                    int AaBb = 0;
                    int AAbb = 0;
                    int aaBB = 0;
                    int aaBb = 0;
                    int aabb = 0;
                    int Aabb = 0;
                    for (int t = 0; t < x.PeoNum; t++)
                    {
                        if (x.controlData[t, i] == 0 && x.controlData[t, j] == 0)
                        {
                            AABB++;
                        }
                        if (x.controlData[t, i] == 0 && x.controlData[t, j] == 1)
                        {
                            AABb++;
                        }
                        if (x.controlData[t, i] == 1 && x.controlData[t, j] == 0)
                        {
                            AaBB++;
                        }
                        if (x.controlData[t, i] == 1 && x.controlData[t, j] == 1)
                        {
                            AaBb++;
                        }
                        if (x.controlData[t, i] == 0 && x.controlData[t, j] == 2)
                        {
                            AAbb++;
                        }
                        if (x.controlData[t, i] == 2 && x.controlData[t, j] == 0)
                        {
                            aaBB++;
                        }
                        if (x.controlData[t, i] == 2 && x.controlData[t, j] == 1)
                        {
                            aaBb++;
                        }
                        if (x.controlData[t, i] == 2 && x.controlData[t, j] == 2)
                        {
                            aabb++;
                        }
                        if (x.controlData[t, i] == 1 && x.controlData[t, j] == 2)
                        {
                            Aabb++;
                        }
                    }
                    double LinkAB = AABB + AABb + AaBB + AaBb;
                    double PA = (double)(2 * AABB + 2 * AABb + AaBB + 2 * AAbb + Aabb + AaBb) / (2 * x.PeoNum);
                    double PB = (double)(2 * AABB + AABb + 2 * AaBB + 2 * aaBB + aaBb + AaBb) / (2 * x.PeoNum);
                    double PAB = 0.25;
                    double Pab = 0.25;
                    double PAb = 0.25;
                    double PaB = 0.25;
                    double LastPar = 0;
                    double AB = 0;
                    double Ab = 0;
                    par = 0.5;
                    double DM=0;
                    do{
                         LastPar = par;
                         PAB = (2 * AABB + AABb + AaBB + AaBb * par) / (2 * x.PeoNum);
                         PAb = (2 * AAbb + AABb + Aabb + AaBb * (1 - par)) / (2 * x.PeoNum);
                         PaB = (2 * aaBB + AaBB + aaBb + AaBb * (1 - par)) / (2 * x.PeoNum);
                         Pab = (2 * aabb + Aabb + aaBb + AaBb * par) / (2 * x.PeoNum);
                         par = PAB * Pab / (PAB * Pab + PAb * PaB);
                         AB = 2 * AABB + AABb + AaBB + AaBb * par;
                         Ab = 2 * AAbb + AABb + Aabb + AaBb * (1 - par);
                    } while (Math.Abs(par-LastPar) > 0.000001);
                   double D = PAB * Pab - PAb * PaB;
                    if (D > 0)
                    {
                        if (PA*(1-PB)> (1 - PA) * PB)
                        {
                            DM = (1 - PA) * PB;
                        }
                        else
                        {
                            DM = PA * (1 - PB);

                        }
                    }
                    if (D < 0)
                    {
                        if (-PA* PB > -(1-PA)*(1-PB))
                        {
                            DM = -PA * PB;
                        }
                        else
                        {
                            DM = -(1 - PA) * (1 - PB);
                        } 
                    }
                    Dx[i,j] = D / DM;
                    r2 = D * D / (PAB * PAb * PaB * Pab);
                    if (r2 < 0.1)
                        continue;
                    else
                    {
                        linkageData[i, 0] = i;
                        linkageData[i, 1] = j;
                        linkageData[i, 2] = PAB;
                        linkageData[i, 3] = PAb;
                        linkageData[i, 4] = PaB;
                        linkageData[i, 4] = Pab;
                    }       
                }
            }
        }
    }
}
