using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;
using System.IO;

namespace IBD
{
    public struct HaplotypePair
    {
        public int firstHaplotypeID;
        public int secondHaplotypeID;
        public double firstFrequency;
        public double secondFrequency;
        public HaplotypePair(int newFirst, int newSecond, double newFirstF, double newSecondF)
        {
            this.firstHaplotypeID = newFirst;
            this.firstFrequency = newFirstF;
            this.secondHaplotypeID = newSecond;
            this.secondFrequency = newSecondF;
        }
    }
    class Block  //Haplotype information in the block
    {
        private int start; // start position of the block
        private int end;  //end position of the block
        private int length; //length of the block
        private int[,] haplotype; //haplotype string of all haplotypes in the block
                                  //[i, ] stores the ith haplotype.
        private double[] frequency;  //haplotype frequency of all the haplotype.
        private int startSnpID;
        private int endSnpID;
        public int GetHaplotypeCount()
        {
            return this.frequency.Length;
        }
        public void SetStartSnpID(int newStart)
        {
            this.startSnpID = newStart;
        }
        public void SetEndSnpID(int newEnd)
        {
            this.endSnpID = newEnd;
        }
        public int GetStartSnpID()
        {
            return this.startSnpID;
        }
        public int GetEndSnpID()
        {
            return this.endSnpID;
        }
        public Block GetCopy()
        {
            Block newBlock = new Block();
            newBlock.start = this.GetStart();
            newBlock.end = this.GetEnd();
            newBlock.startSnpID = this.startSnpID;
            newBlock.endSnpID = this.endSnpID;
            newBlock.length = this.length;
            newBlock.haplotype = (int[,])this.haplotype.Clone();
            newBlock.frequency = (double[])this.frequency.Clone();
            return newBlock;
        }
        public int GetStart()
        {
            return this.start;
        }
        public void SetStart(int newStart)
        {
            this.start = newStart;
        }
        public int GetEnd()
        {
            return this.end;    
        }
        public void SetEnd(int newEnd)
        {
            this.end = newEnd;
            this.length = this.end - this.start + 1;
        }
        public int Gethaplotype(int haplotypeIndex, int snpIndex)
        {
            if (haplotypeIndex > this.haplotype.GetUpperBound(0))
                haplotypeIndex = this.haplotype.GetUpperBound(0);
            if (snpIndex > this.haplotype.GetUpperBound(1))
                haplotypeIndex = this.haplotype.GetUpperBound(1);
            return this.haplotype[haplotypeIndex, snpIndex];
        }
        public void SetHaplotype (ArrayList newHaplotype)
        {
            this.haplotype = new int[newHaplotype.Count, ((char[])newHaplotype[0]).Length];
            for (int i = 0; i < newHaplotype.Count; i++)
            { 
                char[] singleHaplotye = (char[])newHaplotype[i];
                for (int j = 0; j < singleHaplotye.Length; j++)
                {
                    if (singleHaplotye[j].Equals('2'))
                        this.haplotype[i, j] = 0;
                    else // '1'
                        this.haplotype[i, j] = 1;
                }
            }
        }
        public void SetFrequency(double[] newFrequency)
        {
            this.frequency = new double[newFrequency.Length];
            for (int i = 0; i < newFrequency.Length; i++)
                this.frequency[i] = newFrequency[i];
        }
        public double GetFrequency(int index)
        {
            if (index > this.frequency.GetUpperBound(0))
                return this.frequency[this.frequency.GetUpperBound(0)];
            return this.frequency[index];
        }
        public ArrayList GenotypeToHaplotype(int[] genotype, double[,] lingkageData)
        {
            int[,] homoPosition = new int[genotype.Length, 2]; //store the homozygosity position [position, allele]
            for (int i = 0; i < genotype.Length; i++)
                homoPosition[i,0] = -1;
            int homoIndex = 0;
            for (int i = 0; i < genotype.Length; i++)
            {
                if (genotype[i] != 1)
                {
                    homoPosition[homoIndex, 0] = i;
                    homoPosition[homoIndex, 1] = genotype[i];
                    homoIndex++;
                }
            }
            ArrayList candidateHaplotype = new ArrayList();
            for (int i = 0; i < this.haplotype.GetUpperBound(0)+1; i++) //for each haplotype
            {
                for (int j = 0; j < homoPosition.GetUpperBound(0)+1; j++) //for each homoposition
                {
                    if (homoPosition[j, 0] == -1)
                    {
                        candidateHaplotype.Add(i);
                        break;
                    }
                    else
                    {
                        if (homoPosition[j, 1] != this.haplotype[i, homoPosition[j, 0]] * 2)
                        {
                            break;
                        }
                        else if (j == this.haplotype.GetUpperBound(1))
                        {
                            candidateHaplotype.Add(i);
                        }
                    }
                }
            }

            bool[,] combination = new bool[candidateHaplotype.Count, candidateHaplotype.Count]; //store the combination possability for each haplotype combination, false means impossible.
            for (int i = 0; i < candidateHaplotype.Count; i++)
                for (int j = 0; j < candidateHaplotype.Count; j++)
                    combination[i, j] = false; //all set to false at the beginning
            for (int i = 0; i < candidateHaplotype.Count; i++)//fill in the combination table
                                                              //for each candidate haplotype
            {
                for (int j = 0; j < candidateHaplotype.Count; j++) //for each second haplotype
                {
                    for (int k = 0; k < genotype.Length; k++) //for each position
                    {
                        if (haplotype[(int)candidateHaplotype[i], k] + haplotype[(int)candidateHaplotype[j], k] != genotype[k])
                        {
                            combination[i, j] = false;
                            break;
                        }
                        else if (k == genotype.Length - 1)
                            combination[i, j] = true;
                    }
                }
            }
            ArrayList haplotypePairList = new ArrayList();


            for (int i = 0; i < candidateHaplotype.Count; i++)
            {
                for (int j = i; j < candidateHaplotype.Count; j++)
                {
                    if (combination[i, j] == true)
                    {
                        for (int n = 0; n < genotype.Length; n++)
                        {
                            if (haplotype[(int)candidateHaplotype[i], n] + haplotype[(int)candidateHaplotype[j], n] != genotype[n])
                                Console.WriteLine("here\n");
                        }
                        double newFrequency = this.frequency[(int)candidateHaplotype[i]] * this.frequency[(int)candidateHaplotype[j]];
                        HaplotypePair newHaplotypePair = new HaplotypePair((int)candidateHaplotype[i], (int)candidateHaplotype[j], this.frequency[(int)candidateHaplotype[i]], this.frequency[(int)candidateHaplotype[j]]);
                        haplotypePairList.Add(newHaplotypePair);
                    }
                }
            }
            HaplotypePairCompare haplotypePairCompare = new HaplotypePairCompare();
            haplotypePairList.Sort(haplotypePairCompare);
            haplotypePairList.Reverse();

            if (haplotypePairList.Count == 0) //no haplotype pair can be combined to this genotype.
            {

                if (candidateHaplotype.Count != 0) // can find a haplotype from the genotype
                {
                    HaplotypePair newHaplotypePair = new HaplotypePair((int)candidateHaplotype[0], -1, this.frequency[(int)candidateHaplotype[0]], 0.001);
                    haplotypePairList.Add(newHaplotypePair);
                }
                else //can not find one haplotype from the geontype.
                {
                    int[] haplotype1 = new int[genotype.Length];
                    int[] haplotype2 = new int[genotype.Length];
                    Random rdSeed = new Random();
                    for (int snpIndex = 0; snpIndex < this.endSnpID - this.startSnpID + 1; snpIndex++)
                    {
                        if (genotype[snpIndex] == 0)
                        {
                            haplotype1[snpIndex] = 0;
                            haplotype2[snpIndex] = 0;
                        }
                        else if (genotype[snpIndex] == 2)
                        {
                            haplotype1[snpIndex] = 1;
                            haplotype2[snpIndex] = 1;
                        }
                        else //genotype == 0
                        {
                            if (snpIndex == 0 || genotype[snpIndex - 1] != 1)
                            {

                                if (rdSeed.NextDouble() <= 0.5)
                                {
                                    haplotype1[snpIndex] = 0;
                                    haplotype2[snpIndex] = 1;
                                }
                                else
                                {
                                    haplotype2[snpIndex] = 0;
                                    haplotype1[snpIndex] = 1;
                                }
                            }
                            else
                            {
                                if (rdSeed.NextDouble() <= lingkageData[this.startSnpID + snpIndex - 1, 6])
                                {
                                    haplotype1[snpIndex] = haplotype1[snpIndex - 1];
                                    haplotype2[snpIndex] = haplotype2[snpIndex - 1];
                                }
                                else
                                {
                                    haplotype1[snpIndex] = haplotype1[snpIndex - 1];
                                    haplotype2[snpIndex] = haplotype2[snpIndex - 1];
                                }
                            }
                        }
                    }

                    int[,] newHaplotype = new int[this.haplotype.GetUpperBound(0)+2, this.haplotype.GetUpperBound(1)+1];
                    double[] newFrequency = new double[this.frequency.GetUpperBound(0)+2];
                    for (int i = 0; i < newHaplotype.GetUpperBound(0) + 1; i++)
                    {
                        for (int j = 0; j < newHaplotype.GetUpperBound(1) + 1; j++)
                        {
                            if (i == newHaplotype.GetUpperBound(0))
                            {
                                newHaplotype[i, j] = haplotype1[j];
                            }
                            else
                            {
                                newHaplotype[i, j] = this.haplotype[i, j];
                            }
                        }
                        if (i == newHaplotype.GetUpperBound(0))
                            newFrequency[i] = Math.Max(newFrequency[i - 1], 0.001);
                        else
                            newFrequency[i] = this.frequency[i];
                    }
                    this.haplotype = newHaplotype;
                    this.frequency = newFrequency;
                    HaplotypePair newHaplotypePair = new HaplotypePair(this.haplotype.GetUpperBound(0), -1, this.frequency[this.frequency.GetUpperBound(0)], 0.001);
                    haplotypePairList.Add(newHaplotypePair);
                }
                
            }
            return haplotypePairList;
        }
        private class HaplotypePairCompare : IComparer
        {
            public int Compare(Object x, Object y)
            {                                                                       //the order of the blocks
                return (int)(((HaplotypePair)x).firstFrequency * ((HaplotypePair)x).secondFrequency * 100000 - ((HaplotypePair)y).firstFrequency * ((HaplotypePair)y).secondFrequency * 100000);
            }
        }
    }
}
