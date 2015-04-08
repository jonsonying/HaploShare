using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace IBD
{
    class RareBlock
    {
        private int start;
        private int end;
        private int[] sequence;
        private int[] haplotype;
        //public ArrayList shareRarePosition;
        private double[] p;
        private int[] group;
        private double cumulativeP;

        public int GetStart()
        {
            return this.start;
        }  //get the start position of the region
        public void SetStart(int newStart)
        {
            this.start = newStart;
        }   //set the start position of the region
        public int GetEnd()
        {
            return this.end;
        } //get the end position of the region
        public void SetEnd(int newEnd)
        {
            this.end = newEnd;
        } // set the end position of the region
        public void InitialSequence(int length)
        {
            this.sequence = new int[length];        
        }  //initial the sequence of the region
        public int[] GetSequence()
        {
            int[] returnSequence = new int[this.sequence.Length];
            for (int i = 0; i < this.sequence.Length; i++)
                returnSequence[i] = this.sequence[i];
            return returnSequence;
        }
        public int GetSequence(int index)
        {
            if (index > this.sequence.Length)
                return -1;
            else
                return this.sequence[index];
        }
        public void SetSequence(int[] newSequence)
        {
            this.sequence = new int[newSequence.Length];
            for (int i = 0; i < newSequence.Length; i++)
            {
                this.sequence[i] = newSequence[i];
            }
        }
        public void SetSequence(int index, int value)
        {
            this.sequence[index] = value;
        }
        public void InitialHaplotype(int length)
        {
            this.haplotype = new int[length];
        }
        public int[] GetHaplotype()
        {
            int[] returnHaplotype = new int[this.haplotype.Length];
            for (int i = 0; i < this.haplotype.Length; i++)
                returnHaplotype[i] = this.haplotype[i];
            return returnHaplotype;
        }
        public int GetHaplotype(int index)
        {
            if (index > this.haplotype.Length)
                return -1;
            else
                return this.haplotype[index];
        }
        public void SetHaplotype(int[] newHaplotype)
        {
            this.haplotype = new int[newHaplotype.Length];
            for (int i = 0; i < newHaplotype.Length; i++)
            {
                this.haplotype[i] = newHaplotype[i];
            }
        }
        public void SetHaplotype(int index, int value)
        {
            this.haplotype[index] = value;
        }
        public int[] GetGroup()
        {
            int[] returnGroup = new int[this.group.Length];
            for (int i = 0; i < this.group.Length; i++)
            {
                returnGroup[i] = this.group[i];
            }
            return returnGroup;
        }
        public int GetGroup(int index)
        {
            return this.group[index];
        }
        public void SetGroup(int[] newGroup)
        {
            this.group = new int[newGroup.Length];
            for (int i = 0; i < newGroup.Length; i++)
            {
                this.group[i] = newGroup[i];
            }
        }
        public void InitialGroup(int length)
        {
            this.group = new int[length];
        }
        public void SetGroup(int index, int value)
        {
            this.group[index] = value;
        }
        public void InitialP()
        {
            p = new double[this.end - this.start + 1];
        }
        public void SetP(int index, double p)
        {
            this.p[index] = p;
        }
        public double[] GetP()
        {
            double[] returnP = new double[this.p.Length];
            for (int i = 0; i < this.p.Length; i++)
                returnP[i] = this.p[i];
            return returnP;
        }
        public double GetP(int index)
        {
            return this.p[index];        
        }
        public RareBlock GetCopy()
        {
            return (RareBlock)MemberwiseClone();
        }
        public double GetMinusLogP()
        {
            return cumulativeP;
        }
        public void SetMinusLogP()
        {
            double logP = 0;
            for (int i = 0; i < this.GetSequence().Length; i++)
                logP -= Math.Log10(this.GetP(i));
            this.cumulativeP = logP;
        }
        public void Coverage(ArrayList rareBlockList, SnpDataSet dataSet)
        {
            int[] coverage = new int[dataSet.GetSnpCount()];
            int total = 0;
            foreach (RareBlock rareBlock in rareBlockList)
            {
                for (int i = rareBlock.start; i <= rareBlock.end; i++)
                    coverage[i] = 1;
            }
            foreach (int i in coverage)
                if (i == 1)
                    total++;
            Console.WriteLine(total);
        }    
        public bool CheckSame(RareBlock rareBlockToCheck)   //check if two rare blocks are the same or not.
        {
            if (this.GetStart() != rareBlockToCheck.GetStart()) return false;
            if (this.GetEnd() != rareBlockToCheck.GetEnd()) return false;
            if (this.GetGroup().Length != rareBlockToCheck.GetGroup().Length) return false;
            for (int i = 0; i < this.GetGroup().Length; i++)
            {
                if (this.GetGroup(i) != rareBlockToCheck.GetGroup(i))
                    return false;
            }
            return true;
        }
        public bool CheckContain(RareBlock rareBlockToCheck) //check if RareBlockToCheck is in the region of "this" with same group.
        {
            if (this.GetGroup().Length != rareBlockToCheck.GetGroup().Length) return false;
            if (this.GetEnd() > rareBlockToCheck.GetEnd()) return false;
            if (this.GetStart() < rareBlockToCheck.GetStart()) return false;
            for (int i = 0; i < this.GetGroup().Length; i++)
            {
                if (this.GetGroup(i) != rareBlockToCheck.GetGroup(i))
                    return false;
            }
            return true;
        }
        public bool CkeckSubGroup(RareBlock rareBlockToCheck)
        {
            if (this.GetStart() != rareBlockToCheck.GetStart() || this.GetEnd() != rareBlockToCheck.GetEnd()||this.GetGroup().Length <= rareBlockToCheck.GetGroup().Length)
                return false;
            
            int i = 0;
            int j = 0;
            while (j < rareBlockToCheck.GetGroup().Length)
            {
                if (this.GetGroup(i) == rareBlockToCheck.GetGroup(j))
                {
                    i++;
                    j++;
                }
                else if (this.GetGroup(i) < rareBlockToCheck.GetGroup(j))
                    i++;
                else
                    return false;
                if (i == this.GetGroup().Length && j != rareBlockToCheck.GetGroup().Length)
                    return false;
            }
            return true;
        }
        private bool CheckCommon(SnpDataSet dataSet) //check snp data between rare blocks
        { 
            int[] error = new int[this.group.Length];   //error data.
            this.sequence = new int[this.end - this.start + 1];

            for( int e = 0; e < error.Length; e++ )     // all is 0 at the begining.
                error[e] = 0;
            for (int i = start ; i < end; i++)  //loop from start to end for snp
            {
                int count0 = 0;                                     //count how many 0
                int count2 = 0;                                     //count how many 2    
                int[] dataTmp = new int[this.group.Length];         //store the data, avoid multi data accessing
                //for (int d = 0; d < dataTmp.Length; d++)            //
                //    d = 0;
                int j = 0;
                foreach (int ind in this.group)                     // for each guy in group ,check this snp data
                {
                    int data = dataSet.GetSnp(ind, i);
                    if (data == 2)
                        count2++;
                    else if (data == 0)
                        count0++;
                    dataTmp[j] = data;                              //store the data.        
                    j++;
                }
                if (count2 <= count0)                               //two minor allele is the outlier in the group
                {
                    for (int m = 0; m < this.group.Length; m++)
                    {
                        if (dataTmp[m] == 2)
                            error[m]++;
                    }
                    sequence[i - start] = 0;
                }
                else                                                //two major allele is the outlier in the group
                {
                    for (int m = 0; m < this.group.Length; m++)     
                    {
                        if (dataTmp[m] == 0)
                            error[m]++;
                    }
                    sequence[i - start] = 2;
                }
                if (count2 == 0 && count0 == 0)
                    sequence[i - start] = 1;
            }

            ArrayList newGroup = new ArrayList();
            for (int n = 0; n < this.group.Length; n++)
            {
                if (!(error[n] > 0))   //no error or little error, keep this guy.
                    newGroup.Add(this.group[n]);
            }
            if (newGroup.Count < 3) return false;
            else
            {
                this.group = (int[])newGroup.ToArray(typeof(int));
                return true;
            }
        }
        public void Extend(SnpDataSet dataSet)
        {
            int newStart = this.GetStart();
            int newEnd = this.GetEnd();
            ArrayList newEndingSequence = new ArrayList();
            while (newEnd < dataSet.GetSnpCount() - 1)
            {
                int data = -1;
                bool flag = false;
                for (int i = 1; i < this.group.Length; i++)
                {

                    if (data == -1)
                    {
                        if (dataSet.GetSnp(this.group[i], newEnd + 1) != 1)
                            data = dataSet.GetSnp(this.group[i], newEnd + 1);
                    }
                    else if (Math.Abs(dataSet.GetSnp(this.group[i], newEnd + 1) - data) == 2)
                    {
                        flag = true;
                        break;
                    }
                    
                }
                if (flag) break;
                if (data == -1)
                    newEndingSequence.Add(1);
                else if (data == 0)
                    newEndingSequence.Add(0);
                else
                    newEndingSequence.Add(2);

                newEnd++;
            }
            ArrayList newStartingSequence = new ArrayList();
            while (newStart > 0)
            {
                int data = -1;
                bool flag = false;
                for (int i = 1; i < this.GetGroup().Length; i++)
                {

                    if (data == -1)
                    {
                        if (dataSet.GetSnp(this.group[i], newStart - 1) != 1)
                            data = dataSet.GetSnp(this.group[i], newStart - 1);
                    }
                    else if (Math.Abs(dataSet.GetSnp(this.group[i], newStart - 1) - data) == 2)
                    {
                        flag = true;
                        break;
                    }
                }
                if (flag) break;
                if (data == -1)
                    newStartingSequence.Add(1);
                else if (data == 0)
                    newStartingSequence.Add(0);
                else
                    newStartingSequence.Add(2);
                
                newStart--;
            }
            newStartingSequence.Reverse();
            int[] newWholeSequence = new int[this.sequence.Length + newStartingSequence.Count+ newEndingSequence.Count];
            newStartingSequence.CopyTo(newWholeSequence);
            this.GetSequence().CopyTo(newWholeSequence, newStartingSequence.Count);
            newEndingSequence.CopyTo(newWholeSequence, newStartingSequence.Count + this.sequence.Length);
            this.SetSequence(newWholeSequence);
            this.SetStart(newStart);
            this.SetEnd(newEnd);
            if (this.GetEnd() - this.GetEnd() > 170)
                ;
            while (newEnd < dataSet.GetSnpCount())
            {
                ArrayList group0 = new ArrayList();
                ArrayList group1 = new ArrayList();
                ArrayList group2 = new ArrayList();
                for (int i = 0; i < this.GetGroup().Length; i++)
                {
                    if (this.GetGroup(i) == 0)
                        group0.Add(this.GetGroup(i));
                    else if (this.GetGroup(i) == 2)
                        group2.Add(this.GetGroup(i));
                    else
                        group1.Add(this.GetGroup(i));
                }
                
            }
        }
        public void Combine(int i, ArrayList rareBlockList, ArrayList combinedBlockList, SnpDataSet dataSet)
        {
            
            if (i == rareBlockList.Count - 1) return; //reach the end of the rareBlockList.
            ArrayList newGroup = new ArrayList();
            ArrayList remainGroup = new ArrayList();
            bool flag = true;
            for (int j = i + 1; j < i + 50 && j < rareBlockList.Count; j++) // find combined block from this block to 100 blocks away.
            {                
                RareBlock block = new RareBlock();
                //Console.WriteLine("{0} {1}", i, j);
                block = (RareBlock)((RareBlock)rareBlockList[j]).MemberwiseClone();

                if (block.GetStart() < this.GetEnd()) continue;
                int m = 0;
                int n = 0;
                newGroup.Clear();
                remainGroup.Clear();
                for (m = 0; m < this.GetGroup().Length; m++ )  //find the intersection of the two group of individuals.
                {
                    if (n == block.GetGroup().Length)
                    {
                        remainGroup.Add(this.GetGroup(m));
                        break;

                    }
                    if (this.group[m] == block.GetGroup(n))
                    {
                        newGroup.Add(this.group[m]);
                        n++;
                    }
                    else if (this.group[m] > block.GetGroup(n))
                    {
                        n++;
                        m--;
                    }
                    else
                        remainGroup.Add(this.GetGroup(m));
                }

                if (newGroup.Count > 2) //more than 1 guys appare in the two rare blocks
                {
                    //Console.Write("{0} {1} {2}\t", i, j, newGroup.Count);
                    //foreach (int ind in newGroup)
                    //    Console.Write("{0} ",ind);
                    //Console.WriteLine();
                    //Console.ReadKey();
                    RareBlock newRareBlock = new RareBlock();
                    newRareBlock.SetStart(this.GetStart());
                    newRareBlock.SetEnd(block.GetEnd());
                    newRareBlock.SetGroup((int[])newGroup.ToArray(typeof(int)));
                    if (newRareBlock.CheckCommon(dataSet)) //check the snp data between the rare blocks
                    //if (true)
                    {
                        if (newRareBlock.GetGroup().Length == 3)
                        {
                            newRareBlock.Extend(dataSet);
                            if (newRareBlock.GetEnd() - newRareBlock.GetStart() > 170)
                            {
                                if (combinedBlockList.Count > 0)
                                {
                                    RareBlock lastRareBlock = (RareBlock)((RareBlock)combinedBlockList[combinedBlockList.Count - 1]).MemberwiseClone();

                                    if (!lastRareBlock.CheckSame(this))
                                    {
                                        if (this.CheckContain(lastRareBlock))
                                            combinedBlockList[combinedBlockList.Count - 1] = this.GetCopy();
                                        else
                                            combinedBlockList.Add(this.GetCopy());
                                    }
                                }
                                else
                                    combinedBlockList.Add(this.GetCopy());
                            }
                        }
                        else
                            newRareBlock.Combine(j, rareBlockList, combinedBlockList, dataSet);//if pass, continue combining.   
                        //break;       //When we find another block to combine to this block ,stop searching for this block.
                                     //Time saving a lot, can miss some, discuss later on.
                        if (this.group.Length - newRareBlock.group.Length == 0)
                        {
                            flag = false;
                            break;
                        }
                        //else if (newRareBlock.group.Length / (this.group.Length - newRareBlock.group.Length) >= 3) break;
                    }
                    /*
                    if (this.group.Length - newRareBlock.group.Length >= 3)
                    {
                        RareBlock remainGroupBlock = new RareBlock();
                        remainGroupBlock.start = this.start;
                        remainGroupBlock.end = this.end;
                        remainGroupBlock.sequence = this.sequence;
                        remainGroupBlock.group = (int[])remainGroup.ToArray(typeof(int));
                        remainGroupBlock.Combine(j, rareBlockList, combinedBlockList, dataSet);
                    }
                     */
                    
                }
            }
            if (flag && (this.GetEnd() - this.GetEnd()) > 70)
            {
                if (combinedBlockList.Count > 0)
                {
                    RareBlock lastRareBlock = (RareBlock)combinedBlockList[combinedBlockList.Count - 1];
                 
                    if (!lastRareBlock.CheckSame(this))
                    {
                        if (this.CheckContain(lastRareBlock))
                            combinedBlockList[combinedBlockList.Count - 1] = this.GetCopy();
                        else
                            combinedBlockList.Add(this.GetCopy());
                    }
                }
                else
                    combinedBlockList.Add(this.GetCopy());
            }    
        }
        //private void MarcovModel(double p,int start, int end, int[] haplotype, int[] linkageData)
        public bool Pvalue(double[,] linkageData, ArrayList pairwiseIBDList, int index)
        {
            this.InitialP();
            this.InitialHaplotype(this.GetSequence().Length);
            this.SetHaplotype(this.GetSequence());
            this.SetP(0, 0);
            for (int i = 0; i < this.GetP().Length; i++)
            {
                if (i == 0)  //first snp of the region, set the shared haplotype to be the major allele.
                {
                    if (this.GetSequence(i) > 1)
                        this.SetHaplotype(i, 2);
                    else
                        this.SetHaplotype(i, 0);
                    this.SetP(0, 1);
                    continue;
                }
                if (this.haplotype[i] != 1) //if the haplotype is determined (0 for share common, 2 for share rare, 1 for uncertain)
                {
                    int currentHap = this.GetHaplotype(i);
                    int lastHap = this.GetHaplotype(i - 1);
                    if (lastHap == 1)
                    {
                        this.SetP(i, 1);
                        continue;
                    }
                    if (currentHap == 0) //for cases to check, major->major, major->minor, minor->major,minor->minor;
                    {
                        if (lastHap == 0)
                            this.SetP(i, linkageData[this.GetStart() + i - 1, 2]);
                        else
                            this.SetP(i, linkageData[this.GetStart() + i - 1, 4]);
                    }
                    else
                    {
                        if (lastHap == 0)
                            this.SetP(i, linkageData[this.GetStart() + i - 1, 3]);
                        else
                            this.SetP(i, linkageData[this.GetStart() + i - 1, 5]);
                    }
                }
                else  //if the haplotype at position i is not determined (the value is 1), but in this case,
                // the previous snp haplotype can determine this haplotype (say, last one is major and major->minor is 0%, then determine that this position is major)
                {
                    this.SetP(i, 1);
                }

            }
     
 //*** now  we should check the break point of the haplotype. (at which position the conversion P is almost zero.(which is simulated to be 1-power(0.95, sampleSize)))    
            int lastBreak = 0;  //last break point position.
            bool flag = false;  //false means that no break point. true means there is at least one break point.
            //int lastBreak = this.GetSequence().Length;
            //return true;
            for (int i = 0; i < this.GetP().Length; i++)
            {
                if (this.GetP(i) <= 0.001)//1.0 / (linkageData.GetUpperBound(0) + 1))
                {
                    if (i - lastBreak < 100)
                    {
                        lastBreak = i;
                    }
                    else
                    {
                        RareBlock newBlock = new RareBlock();
                        newBlock.InitialSequence(i - lastBreak);
                        for (int j = 0; j < newBlock.GetSequence().Length; j++)
                            newBlock.SetSequence(j, this.GetSequence(j + lastBreak));
                        newBlock.SetStart(lastBreak);
                        newBlock.SetEnd(i - 1);
                        newBlock.SetGroup(this.GetGroup());
                        pairwiseIBDList.Insert(index + 1, newBlock);
                        lastBreak = i;
                        
                    }
                    flag = true;
                }
            }
            if (flag)
            {
                pairwiseIBDList.RemoveAt(index);
                return false;
            }
            return true;
            
        }
        public bool PvalueWithMM(double[,] linkageData, ArrayList pairwiseIBDList, int index)
        { 
            this.InitialP();
            this.InitialHaplotype(this.GetSequence().Length);
            this.SetHaplotype(this.GetSequence());
            this.SetP(0, 0);
            for (int i = 0; i < this.GetP().Length; i++)
            {
                if (i == 0)  //first snp of the region, set the shared haplotype to be the major allele.
                {
                    if (this.GetSequence(i) > 1)
                        this.SetHaplotype(i, 2);
                    else
                        this.SetHaplotype(i, 0);
                    this.SetP(0, 1);
                    continue;
                }
                if (this.haplotype[i] != 1) //if the haplotype is determined (0 for share common, 2 for share rare, 1 for uncertain)
                {
                    int currentHap = this.GetHaplotype(i);
                    int lastHap = this.GetHaplotype(i - 1);
                    if (lastHap == 1)
                    {
                        this.SetP(i, 1);
                        continue;
                    }
                    if (currentHap == 0) //for cases to check, major->major, major->minor, minor->major,minor->minor;
                    {
                        if (lastHap == 0)
                            this.SetP(i, linkageData[this.GetStart() + i - 1, 2]);
                        else
                            this.SetP(i, linkageData[this.GetStart() + i - 1, 4]);
                    }
                    else
                    {
                        if (lastHap == 0)
                            this.SetP(i, linkageData[this.GetStart() + i - 1, 3]);
                        else
                            this.SetP(i, linkageData[this.GetStart() + i - 1, 5]);
                    }
                }
                else  //if the haplotype at position i is not determined (the value is 1), but in this case,
                // the previous snp haplotype can determine this haplotype (say, last one is major and major->minor is 0%, then determine that this position is major)
                {
                    if (this.haplotype[i - 1] == 0)
                    {
                        if (linkageData[this.GetStart() + i - 1, 2] == (double)1)
                        {
                            this.SetSequence(i, 0);
                            this.SetHaplotype(i, 0);
                            this.SetP(i, 1);
                            continue;
                        }
                        else if (linkageData[this.GetStart() + i - 1, 3] == (double)1)
                        {
                            this.SetSequence(i, 2);
                            this.SetHaplotype(i, 2);
                            this.SetP(i, 1);
                            continue;
                        }
                    }
                    else
                    {
                        if (linkageData[this.start + i - 1, 4] == (double)1)
                        {
                            this.SetSequence(i, 0);
                            this.SetHaplotype(i, 0);
                            this.SetP(i, 1);
                            continue;
                        }
                        else if (linkageData[this.start + i - 1, 5] == (double)1)
                        {
                            this.SetSequence(i, 2);
                            this.SetHaplotype(i, 2);
                            this.SetP(i, 1);
                            continue;
                        }
                    }
                    if (i == this.GetP().Length - 1)
                    {
                        if (this.GetHaplotype(i - 1) == 0)
                        {
                            if (linkageData[this.GetStart() + i - 1, 2] >= 0.5)
                            {
                                this.SetHaplotype(i, 0);
                                this.SetP(i, linkageData[this.GetStart() + i - 1, 2]);
                            }
                            else
                            {
                                this.SetHaplotype(i, 2);
                                this.SetP(i, linkageData[this.GetStart() + i - 1, 3]);
                            }
                        }
                        else
                        {
                            if (linkageData[this.GetStart() + i - 1, 4] >= 0.5)
                            {
                                this.SetHaplotype(i, 0);
                                this.SetP(i, linkageData[this.GetStart() + i - 1, 4]);
                            }
                            else
                            {
                                this.SetHaplotype(i, 2);
                                this.SetP(i, linkageData[this.GetStart() + i - 1, 5]);
                            }       
                        }
                        continue;
                    }
                    //if still here, means that we should decide the haplotype by p value.
                    int j = i;  //j is the start positon of uncertain haplotype region
                    while (i < this.p.Length - 1)
                    {
                        if (this.haplotype[i + 1] == 1)
                        {
                            //this.p[i + 1] = 1;
                            i++;
                        }
                        else
                        {
                            i++;
                            break;
                        }
                        //determine the length of the uncertain haplotype region.
                    }
                    int length = i - j; //the length, now i is the end of the region+1
                    /*
                                        if (length >= 3)  // if the region is >=3, it should be the rareblock which we start the searching.
                                        { 
                                            for (int k = 0; k < length; k++)
                                                this.haplotype[j+k] = 2;
                                            i =j -1;  //go back to the start position of the uncertain region.
                                            continue;
                                        }
                    */
                    //this will cause many impossible combination, omitted now. 7th,Nov,2008
                    //if not, we should calculate which haplotype is most likely to be shared by these individuals.
                    double[,] markovP = new double[length + 1, 4];  //marcov chain data, for each row stores: 
                    //[,0]:major of this snp will choose major(0) or minor(2), [,1]the p value of this chosen
                    //[,2]:major of this snp will choose major(0) or minor(2), [,3]the p value of this chosen                                        
                    markovP[length, 0] = haplotype[i];
                    markovP[length, 2] = haplotype[i];
                    if (i < this.GetP().Length - 1)
                    {
                        if (markovP[length, 0] == 0)
                        {
                            markovP[length, 1] = linkageData[j + length - 1, 2];
                            markovP[length, 3] = linkageData[j + length - 1, 4];
                        }
                        else
                        {
                            markovP[length, 1] = linkageData[j + length - 1, 3];
                            markovP[length, 3] = linkageData[j + length - 1, 5];
                        }
                    }//fill the last row first, coz which one should follow is determined at the end of the uncertain region.
                    else
                    {
                        markovP[length, 1] = 0.5;
                        markovP[length, 3] = 0.5;
                    }
                    for (int m = length - 1; m >= 0; m--) //then fill the row from bottom to top.
                    {
                        double AB = linkageData[j + m - 1 + this.start, 2] * markovP[m + 1, 1];
                        double Ab = linkageData[j + m - 1 + this.start, 3] * markovP[m + 1, 3];
                        if (AB >= Ab)
                        {
                            markovP[m, 0] = 0;
                            markovP[m, 1] = AB;
                        }
                        else
                        {
                            markovP[m, 0] = 2;
                            markovP[m, 1] = Ab;
                        }
                        double aB = linkageData[j + m - 1 + this.start, 4] * markovP[m + 1, 1];
                        double ab = linkageData[j + m - 1 + this.start, 5] * markovP[m + 1, 3];
                        if (aB >= ab)
                        {
                            markovP[m, 2] = 0;
                            markovP[m, 3] = aB;
                        }
                        else
                        {
                            markovP[m, 2] = 2;
                            markovP[m, 3] = ab;
                        }
                    }// finished filling the marcov table.
                    /* if (this.sequence[j - 1] == 0)
                       {
                           this.p[i] = (int)markovP[0, 1];
                       }
                       else
                       {
                           this.p[i] = (int)markovP[0, 3];
                       }
                       */
                    int n = 0;
                    while (n < length) //read from the marcov table and assign all the haplotype with shared region.
                    {
                        if (this.haplotype[j - 1 + n] == 0)
                            this.SetHaplotype(j + n, (int)markovP[n, 0]);
                        else
                            this.SetHaplotype(j + n, (int)markovP[n, 2]);
                        n++;
                    }
                    i = j - 1; //go back the start position of the uncertain region and caculate the p value (now the uncertain region is certain).
                }

            }
     
 //*** now  we should check the break point of the haplotype. (at which position the conversion P is almost zero.(which is simulated to be 1-power(0.95, sampleSize)))    
            int lastBreak = 0;  //last break point position.
            bool flag = false;  //false means that no break point. true means there is at least one break point.
            //int lastBreak = this.GetSequence().Length;
            //return true;
            for (int i = 0; i < this.GetP().Length; i++)
            {
                if (this.GetP(i) <= 0.001)//1.0 / (linkageData.GetUpperBound(0) + 1))
                {
                    if (i - lastBreak < 100)
                    {
                        lastBreak = i;
                    }
                    else
                    {
                        RareBlock newBlock = new RareBlock();
                        newBlock.InitialSequence(i - lastBreak);
                        for (int j = 0; j < newBlock.GetSequence().Length; j++)
                            newBlock.SetSequence(j, this.GetSequence(j + lastBreak));
                        newBlock.SetStart(lastBreak);
                        newBlock.SetEnd(i - 1);
                        newBlock.SetGroup(this.GetGroup());
                        pairwiseIBDList.Insert(index + 1, newBlock);
                        lastBreak = i;
                        
                    }
                    flag = true;
                }
            }
            if (flag)
            {
                pairwiseIBDList.RemoveAt(index);
                return false;
            }
            return true;
            
        }
        public double GetRandomP(SnpDataSet dataSet)
        {
            double minusLogRandonP = 0;
            for (int i = 0; i < this.GetHaplotype().Length; i++)
            { 
                if (this.GetHaplotype(i) == 2)
                    minusLogRandonP += -Math.Log10(dataSet.GetMAF(i + this.GetStart()));
                else
                    minusLogRandonP += -Math.Log10((double)1 - dataSet.GetMAF(i + this.GetStart()));
            }
            return minusLogRandonP;
        }
        public double GetRandomModelLikelihood(SnpDataSet dataSet)
        {
            double minorLikelihood = 0;
            double majorLikelihood = 0;
            for (int i = 0; i < this.GetSequence().Length; i++)
            {
                if (this.GetHaplotype(i) == 0)
                    majorLikelihood += (double)1 / (1 - dataSet.GetMAF(i + this.GetStart()));
                else
                    minorLikelihood += (double)1 / dataSet.GetMAF(i + this.GetStart());
            }
            return (double)Math.Abs(majorLikelihood / (majorLikelihood + minorLikelihood));
        }
        public double GetLinkageModelLikelihood()
        {
            double minorLikelihood = 0;
            double majorLikelihood = 0;
            for (int i = 0; i < this.GetSequence().Length; i++)
            {
                if (this.GetHaplotype(i) == 0)
                    majorLikelihood += (double)1/this.GetP(i);
                else
                    minorLikelihood += (double)1/this.GetP(i);
            }
            return (double)Math.Abs((majorLikelihood /( majorLikelihood + minorLikelihood)));
        }
    }
}
