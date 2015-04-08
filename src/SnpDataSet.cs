using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;
using System.IO;

namespace IBD
{
    class SnpDataSet  //Two data set: Case and Control.
    {
        public int[] individual;
        private string[] individualName;
        private int snpCount;
        private int[,] snp;
        private double[] maf; //Minor allele frequency
        private int[] physDistance;
        private int[] geneDistance;

        public void InitialIndividual(int size)
        {
            this.individual = new int[size];
        }
        public void SetIndividual(int[] newIndividual)
        {
            this.individual = new int[newIndividual.Length];
            for (int i = 0; i < newIndividual.Length; i++)
                this.individual[i] = newIndividual[i];
        }
        public void SetIndividual(int index, int value)
        {
            this.individual[index] = value;
        }
        public int[] GetIndividual()
        {
            int[] returnIndividuals = new int[this.individual.Length];
            for (int i = 0; i < this.individual.Length; i++)
                returnIndividuals[i] = this.individual[i];
            return returnIndividuals;
        }
        public int GetIndividual(int index)
        {
            return this.individual[index];
        }
        public void SetSnpCount(int count)
        {
            this.snpCount = count;
        }
        public int GetSnpCount()
        {
            return this.snpCount;
        }
        public void InitialSnp(int individualSize, int snpSize)
        {
            this.snp = new int[individualSize, snpSize];
        }
        public void SetSnp(int individualId, int snpId, int value)
        {
            this.snp[individualId, snpId] = value;
        }
        public void SetSnp(int[,] newSnpData)
        {
            this.snp = new int[newSnpData.GetUpperBound(0)+1, newSnpData.GetUpperBound(1)+1];
            for (int i = 0; i < newSnpData.GetUpperBound(0)+1; i++)
                for (int j = 0; j < newSnpData.GetUpperBound(1); j++)
                    this.snp[i, j] = newSnpData[i, j];
        }
        public int GetSnp(int individualId, int snpId)
        {
            return this.snp[individualId, snpId];
        }
        public void SetIndividualName(string[] newIndividualName)
        {
            this.individualName = new string[newIndividualName.Length];
            for (int i = 0; i < newIndividualName.Length; i++)
                this.individualName[i] = newIndividualName[i];
        }
        public void InitialMAF()
        {
            this.maf = new double[this.GetSnpCount()];
            for (int i = 0; i < this.GetSnpCount(); i++)
            { 
                int majorAlleleCount = 0;
                int minorAlleleCount = 0;
                for (int j = 0; j < this.GetIndividual().Length; j++)
                {
                    if (this.GetSnp(j, i) == 2)
                        minorAlleleCount += 2;
                    else if (this.GetSnp(j, i) == 1)
                    {
                        minorAlleleCount++;
                        majorAlleleCount++;
                    }
                    else
                        majorAlleleCount++;
                }
                double frequency = (double)minorAlleleCount / ((double)majorAlleleCount + (double)minorAlleleCount);
                this.maf[i] = frequency;
            }
        }
        public double GetMAF(int snpId)
        {
            return this.maf[snpId];
        }
        
      
        public void Linkage(double[,] linkageData)
        {
            double par;
            double[,] lD;
            double r2;
            
            for (int i = 0; i < this.GetSnpCount(); i++)
            {
                for (int j = i + 1; j < this.GetSnpCount(); j++)
                {
                    int AABB = 0;
                    int AABb = 0;
                    int AAbb = 0;
                    int AaBB = 0;
                    int AaBb = 0;
                    int Aabb = 0;
                    int aaBB = 0;
                    int aaBb = 0;
                    int aabb = 0;
                    for (int t = 0; t < this.GetIndividual().Length; t++)
                    {
                        if (this.GetSnp(t,i)== 0)
                        {
                            if (this.GetSnp(t, j) == 0)
                                AABB++;
                            else if (this.GetSnp(t, j) == 1)
                                AABb++;
                            else
                                AAbb++;
                        }
                        else if (this.GetSnp(t, i) == 1)
                        {
                            if (this.GetSnp(t,j) == 0)
                                AaBB++;
                            else if (this.GetSnp(t, j) == 1)
                                AaBb++;
                            else
                                Aabb++;
                        }
                        else
                        {
                            if (this.GetSnp(t, j) == 0)
                                aaBB++;
                            else if (this.GetSnp(t, j) == 1)
                                aaBb++;
                            else
                                aabb++;                            
                        }
                    }
                    double linkAB = AABB + AABb + AaBB + AaBb;
                    double pA = (double)(2 * AABB + 2 * AABb + AaBB + 2 * AAbb + Aabb + AaBb) / (2 * this.individual.Length);
                    double pB = (double)(2 * AABB + AABb + 2 * AaBB + 2 * aaBB + aaBb + AaBb) / (2 * this.individual.Length);
                    double pAB = 0.25;
                    double pab = 0.25;
                    double pAb = 0.25;
                    double paB = 0.25;
                    double lastPar = 0;
                    double AB = 0;
                    double Ab = 0;
                    par = 0.5;
                    double dPrime = 0;
                    double zeroExpectedProb = 1 - Math.Pow(0.95, (double)1/this.GetIndividual().Length*2);
                    do
                    {
                        lastPar = par;
                        pAB = (2 * AABB + AABb + AaBB + AaBb * par) / (2 * this.individual.Length);
                        pAb = (2 * AAbb + AABb + Aabb + AaBb * (1 - par)) / (2 * this.individual.Length);
                        paB = (2 * aaBB + AaBB + aaBb + AaBb * (1 - par)) / (2 * this.individual.Length);
                        pab = (2 * aabb + Aabb + aaBb + AaBb * par) / (2 * this.individual.Length);
                        par = pAB * pab / (pAB * pab + pAb * paB);
                        AB = 2 * AABB + AABb + AaBB + AaBb * par;
                        Ab = 2 * AAbb + AABb + Aabb + AaBb * (1 - par);
                    } while (Math.Abs(par - lastPar) > zeroExpectedProb);
                    double dValue = pAB * pab - pAb * paB;
                    if (dValue >= 0)
                    {
                        if (pA * (1 - pB) > (1 - pA) * pB)
                        {
                            dPrime = (1 - pA) * pB;
                        }
                        else
                        {
                            dPrime = pA * (1 - pB);
                        }
                    }
                    else
                    {
                        if (-pA * pB > -(1 - pA) * (1 - pB))
                        {
                            dPrime = -pA * pB;
                        }
                        else
                        {
                            dPrime = -(1 - pA) * (1 - pB);
                        }
                    }
                   

                    r2 = dValue * dValue / (pAB * pAb * paB * pab);
                    if (true)
                    {
                        linkageData[i, 0] = i;
                        linkageData[i, 1] = j;
                        if (pAB + pAb > (double)0)
                        {
                            linkageData[i, 2] = pAB / (pAB + pAb);
                            linkageData[i, 3] = pAb / (pAB + pAb);
                            if (pAB / (pAB + pAb) < zeroExpectedProb)
                            {
                                linkageData[i, 2] = zeroExpectedProb;
                                linkageData[i, 3] = 1 - zeroExpectedProb;
                            }
                            if (pAb / (pAB + pAb) < zeroExpectedProb)
                            {
                                linkageData[i, 3] = zeroExpectedProb;
                                linkageData[i, 2] = 1 - zeroExpectedProb;
                            }
                        }
                        else
                        {
                            linkageData[i, 2] = pB;
                            linkageData[i, 3] = 1 - pB;
                        }
                        if (paB + pab > (double)0)
                        {
                            linkageData[i, 4] = paB / (paB + pab);
                            linkageData[i, 5] = pab / (paB + pab);
                            linkageData[i, 6] = par;

                            if (paB / (paB + pab) < zeroExpectedProb)
                            {
                                linkageData[i, 4] = zeroExpectedProb;
                                linkageData[i, 5] = 1 - zeroExpectedProb;
                            }
                            if (pab / (paB + pab) < zeroExpectedProb)
                            {
                                linkageData[i, 5] = zeroExpectedProb;
                                linkageData[i, 4] = 1 - zeroExpectedProb;
                            }                        
                        }
                        else
                        {
                            linkageData[i, 4] = pB;
                            linkageData[i, 5] = 1 - pB;
                        }
                        break;
                    }
                }
            }
        }
        public void SetIndividualList(ArrayList individualList)
        {
            for (int id = 0; id < this.individualName.Length; id++)
            {
               
                int[] newGenotype = new int[this.snpCount];
                for (int snpID = 0; snpID < this.snpCount; snpID++)
                {
                    newGenotype[snpID] = this.snp[id, snpID];
                }
                Individual newIndividual = new Individual(newGenotype, id);
                individualList.Add(newIndividual);
            }
        }
    }
}
