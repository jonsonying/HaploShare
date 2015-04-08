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
        public string[] individualName;
        public int snpCount;
        public int[,] snp;
        public void WindowSlide(ArrayList shareWindowList)  
        {
            
            int windowId = 0;
            ArrayList allIndividuals = new ArrayList();     //put all the individuals into the arraylist
            foreach (int individual in this.individual)
            {
                allIndividuals.Add(individual);
            }

            for (int i = 0; i <= this.snpCount; i = i + 50)      //Window sliding from 0 to end
            {
                Window newWindow = new Window();
                newWindow.start = i;
                if (i + 10 > this.snpCount)
                {
                    newWindow.end = this.snpCount;
                }
                else 
                {
                    newWindow.end = newWindow.start + 50;
                }
                newWindow.id = windowId;
                windowId ++;
                newWindow.inWindowExtend(allIndividuals, newWindow.start, this, shareWindowList);
                Console.Write("{0}\r\n", i);
                //Start the window slide for each window.
            }
        } //failed function
        public void Simulation(int length, int groupSize, double[,] linkageData)
        {
            Random randomSeed = new Random();
            int[] group = new int[groupSize];
            int start = (int)randomSeed.Next(0, this.snpCount - length);
            for (int i = 0; i < group.Length; i++)
            {
                int newIndividual = (int)(randomSeed.NextDouble() * (double)this.individual.Length);
                bool flag = true;
                for (int j = 0; j < i; j++)
                {
                    if (group[j] == newIndividual)
                    {
                        flag = false;
                        i--;
                    }
                }
                if (flag)
                    group[i] = newIndividual;
            }
            // generate a group of individuals to do the simulation
            int[] sharedHaplotype = new int[length];
            for (int i = 0; i < group.Length; i++)
            {
                int[] haplotype1 = new int[length];
                int[] haplotype2 = new int[length];
                for (int j = start; j < start + length; j++)
                {
                    if (j == start)
                    {
                        if (this.snp[group[i], j] == 0)
                        {
                            haplotype1[0] = 0;
                            haplotype2[0] = 0;
                        }
                        else if (this.snp[group[i], j] == 2)
                        {
                            haplotype1[0] = 2;
                            haplotype2[0] = 2;
                        }
                        else
                        {
                            int tmp = (int)randomSeed.Next(0, 2);
                            if (tmp == 0)
                            {
                                haplotype1[0] = 2;
                                haplotype2[0] = 0;
                            }
                            else
                            {
                                haplotype1[0] = 0;
                                haplotype2[0] = 2;
                            }
                        }
                    }
                    else
                    {
                        if (this.snp[group[i], j] == 0)   // two major: BB
                        {
                            haplotype1[j - start] = 0;
                            haplotype2[j - start] = 0;
                        }
                        else if (this.snp[group[i], j] == 2) //two nimor: bb
                        {
                            haplotype1[j - start] = 2;
                            haplotype2[j - start] = 2;
                        }
                        else                                // major+minor: Bb
                        {
                            if (this.snp[group[i], j - 1] != 1)  //AABb or aaBb, half chance
                            {
                                int tmp = (int)randomSeed.Next(0, 2);
                                if (tmp == 0)
                                    haplotype1[j - start] = 0;
                                else
                                    haplotype2[j - start] = 2;
                            }
                            else                       //AaBb, the percentage for AB+ab is stroed in linkage data [snp, 6];
                            {
                                double tmp = randomSeed.NextDouble();
                                if (tmp < linkageData[i, 6])
                                {
                                    haplotype1[j - start] = haplotype1[j - start - 1];
                                    haplotype2[j - start] = haplotype2[j - start - 1];
                                }
                                else
                                {
                                    haplotype1[j - start] = haplotype2[j - start - 1];
                                    haplotype2[j - start] = haplotype1[j - start - 1];                                   
                                }
                            }
                        }
                    }
                }
                if (i == 0)
                    haplotype2.CopyTo(sharedHaplotype, 0);
                else
                {
                    sharedHaplotype.CopyTo(haplotype1, 0);
                    for (int j = start; j < start + length; j++)
                    {
                        if (haplotype1[j - start] == 0)
                        {
                            if (haplotype2[j - start] == 0)
                                this.snp[group[i], j] = 0;
                            else
                                this.snp[group[i], j] = 1;
                        }
                        else
                        {
                            if (haplotype2[j - start] == 0)
                                this.snp[group[i], j] = 1;
                            else
                                this.snp[group[i], j] = 2;
                        }
                    }
                }
            }
            
            FileStream bFile = new FileStream(@"F:\test1\result1.txt", FileMode.Create);
            StreamWriter sw = new StreamWriter(bFile);
            
            for (int i = 0; i < 5; i++)
            {
                sw.Write("{0}\t", group[i]);
            }
            sw.Write("\r\n{0}\r\n",start);
            for (int t = 0; t < this.snpCount; t++)
            {
                if (t < start|| t >= start+200)
                    sw.Write(0);
                else
                    sw.Write("{0}", sharedHaplotype[t-start]);
            }
            sw.Write("\r\n");
            for (int i = 0; i < 5; i++)
            {
                //sw.Write(individualName[i]);
                for (int t = 0; t < this.snpCount; t++)
                {
                    sw.Write("{0}", this.snp[group[i], t]);
                }
                sw.Write("\r\n");
            }
            sw.Write("\r\n");
            for (int i = 0; i < this.individual.Length; i++)
            {
                //sw.Write(individualName[i]);
                for (int t = 0; t < this.snpCount; t++)
                {
                    sw.Write("{0}", this.snp[i, t]);
                }
                sw.Write("\r\n");
            }
            sw.Close();
            
            //Console.ReadKey();

        }
        public void SearchRare(ArrayList rareBlockList)
        {
            for (int s = 0; s < snpCount-2; s++)
            {
                ArrayList group = new ArrayList();
                int end = s;
                foreach (int i in this.individual)
                {
                    if (this.snp[i, s] * this.snp[i, s + 1] * this.snp[i, s + 2] > 0)
                    {
                        group.Add(i);
                    }
                }
                if (group.Count >= 3)
                {
                   // Console.Write("{0} {1}\r\n", s, group.Count);
                    end +=2;
                    Boolean flag = true;

                    while (flag)
                    {
                        if (end == this.snpCount-1) break;
                        foreach (int j in group)
                        {
                            if (this.snp[j,end+1] == 0)
                            {
                                //end++;
                                flag = false;
                                break;                                
                            }
                            
                        }
                        end++;
                        
                    }

                    RareBlock newRareBlock = new RareBlock();
                    newRareBlock.SetStart(s);
                    newRareBlock.SetEnd(end-1);
                    newRareBlock.InitialSequence(newRareBlock.end-newRareBlock.start+1);
                    for (int i = 0; i < newRareBlock.sequence.Length; i++)
                        newRareBlock.SetSequence(i, 1);
                    newRareBlock.group = (int[])group.ToArray(typeof(int));
                    RareBlock lastRareBlock = new RareBlock();
                    rareBlockList.Add(newRareBlock);
                    /*
                    if (newRareBlock.group.Length < this.individual.Length * 100 / 100)
                    {
                        if (rareBlockList.Count < 1)
                            rareBlockList.Add(newRareBlock);
                        else
                        {
                            lastRareBlock = (RareBlock)rareBlockList[rareBlockList.Count - 1];
                            //lastRareBlock2 = (RareBlock)rareBlockList[rareBlockList.Count - 2];
                            if (lastRareBlock.end == newRareBlock.end && lastRareBlock.group.Length == newRareBlock.group.Length)
                            {
                                int i = 0;
                                int j = 0;
                                for (i = 0; i < newRareBlock.group.Length; i++)
                                {
                                    if (lastRareBlock.group[j] == newRareBlock.group[i])
                                        j++;
                                    else break;
                                }
                                if (i != newRareBlock.group.Length)
                                    rareBlockList.Add(newRareBlock);

                            }
                            else
                                rareBlockList.Add(newRareBlock);
                        }
                    }
                     */

                }
            }
        }
        public void Linkage(double[,] linkageData)
        {
            double par;
            double[,] lD;
            double r2;
            //lD = new double[this.snpCount, this.snpCount];
            for (int i = 0; i < this.snpCount; i++)
            {
                for (int j = i + 1; j < this.snpCount; j++)
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
                    for (int t = 0; t < this.individual.Length; t++)
                    {
                        if (this.snp[t, i] == 0)
                        {
                            if (this.snp[t, j] == 0)
                                AABB++;
                            else if (this.snp[t, j] == 1)
                                AABb++;
                            else
                                AAbb++;
                        }
                        else if (this.snp[t, i] == 1)
                        {
                            if (this.snp[t, j] == 0)
                                AaBB++;
                            else if (this.snp[t, j] == 1)
                                AaBb++;
                            else
                                Aabb++;
                        }
                        else
                        {
                            if (this.snp[t, j] == 0)
                                aaBB++;
                            else if (this.snp[t, j] == 1)
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
                    } while (Math.Abs(par - lastPar) > 0.000001);
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
                   //lD[i, j] = dValue / dPrime;
                    r2 = dValue * dValue / (pAB * pAb * paB * pab);
                    if (r2 >= 0)
                    {
                        linkageData[i, 0] = i;
                        linkageData[i, 1] = j;
                        linkageData[i, 2] = pAB / (pAB + pAb);
                        linkageData[i, 3] = pAb / (pAB + pAb);
                        linkageData[i, 4] = paB / (paB + pab);
                        linkageData[i, 5] = pab / (paB + pab);
                        linkageData[i, 6] = par;
                        break;
                    }
                }
            }
        }
    }
}
