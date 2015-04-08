using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;

namespace IBD
{
    class GroupShare
    {
        public ArrayList individualList; // the group sharing individual list
        public ArrayList sharePairList;  // pairs of individuals that share this region, used to do the extension, not for final result. // (1,2), (2,3), (1,4) => 1,2,3,4. individuallist
        public int allSharedFirstBlockID; //first block id of the all shared region
        public int allSharedLastBlockID;  //last block id of the all shared region
        public int partiallySharedFirstBlockID;
        public int partiallySharedLastBlockID;
        public ArrayList haplotypeList; // block haplotypeID list of the all shared region.
        public ArrayList allSharedSnpHaplotype; //shared haplotype for each snp in the all shared region.
        public ArrayList partiallySharedSnpHaplotype;
        public int[,] partiallySharedRegion; //start and end position of the shared region for each individual
        public ArrayList partiallySharedBlockHaplotype; //block haplotypeID for partially shared region.
        public double AValue;


        public void SetAllSharedSnpHaplotype(ArrayList individualList, BlockDictionary newBlockDictionary)
        {
            this.allSharedSnpHaplotype = new ArrayList();
            for (int blockID = this.allSharedFirstBlockID; blockID < this.allSharedLastBlockID; blockID++)
            {
                
                for (int inBlockSnpIndex = 0; inBlockSnpIndex < ((Block)(newBlockDictionary.blockList[blockID])).GetEndSnpID() - ((Block)(newBlockDictionary.blockList[blockID])).GetStartSnpID() + 1; inBlockSnpIndex++)
                    allSharedSnpHaplotype.Add(((Block)(newBlockDictionary.blockList[blockID])).Gethaplotype((int)(this.haplotypeList[blockID - this.allSharedFirstBlockID]), inBlockSnpIndex));
                if (blockID == newBlockDictionary.blockList.Count - 1 || blockID == this.allSharedLastBlockID)
                    return;
                for (int nonBlockSnpIndex = ((Block)(newBlockDictionary.blockList[blockID])).GetEndSnpID() + 1; nonBlockSnpIndex < ((Block)(newBlockDictionary.blockList[blockID + 1])).GetStartSnpID(); nonBlockSnpIndex++)
                {
                    int sharedHaplotype = 1;
                    foreach (int individual in this.individualList)
                    { 
                        if ((int)(((Individual)(individualList[individual])).genotype[nonBlockSnpIndex]) != 1)
                            sharedHaplotype = ((int)(((Individual)(individualList[individual])).genotype[nonBlockSnpIndex])) / 2;
                    }
                    allSharedSnpHaplotype.Add(sharedHaplotype);
                }
            }
            
        }
        public void SetPartiallySharedSnpHaplotype(ArrayList individualList, BlockDictionary newBlockDictionary)
        {
            this.partiallySharedSnpHaplotype = new ArrayList();
            for (int blockID = this.partiallySharedFirstBlockID; blockID <= this.partiallySharedLastBlockID; blockID++)
            {
                for (int inBlockSnpIndex = 0; inBlockSnpIndex < ((Block)(newBlockDictionary.blockList[blockID])).GetEndSnpID() - ((Block)(newBlockDictionary.blockList[blockID])).GetStartSnpID() + 1; inBlockSnpIndex++)
                    partiallySharedSnpHaplotype.Add(((Block)(newBlockDictionary.blockList[blockID])).Gethaplotype((int)(this.partiallySharedBlockHaplotype[blockID - this.partiallySharedFirstBlockID]), inBlockSnpIndex));
                if (blockID == newBlockDictionary.blockList.Count - 1 || blockID == this.partiallySharedLastBlockID)
                    return;
                for (int nonBlockSnpIndex = ((Block)(newBlockDictionary.blockList[blockID])).GetEndSnpID() + 1; nonBlockSnpIndex < ((Block)(newBlockDictionary.blockList[blockID + 1])).GetStartSnpID(); nonBlockSnpIndex++)
                {
                    int sharedHaplotype = 1;
                    for (int individual = 0; individual < this.individualList.Count; individual++)
                    {
                        if (this.partiallySharedRegion[individual, 0] > blockID || this.partiallySharedRegion[individual, 1] < blockID)
                            continue;
                        if ((int)(((Individual)(individualList[individual])).genotype[nonBlockSnpIndex]) != 1)
                            sharedHaplotype = ((int)(((Individual)(individualList[individual])).genotype[nonBlockSnpIndex])) / 2;
                    }
                    this.partiallySharedSnpHaplotype.Add(sharedHaplotype);
                }
            }
        }
        public GroupShare GetCopy()
        {
            GroupShare newGroupShare = new GroupShare();
            newGroupShare = (GroupShare)this.MemberwiseClone();
            return newGroupShare;
        }
        public void SetAValue(BlockDictionary newBlockDictionary,ArrayList newIndividualList, double[,] newLinkageData)
        {
            double[] singleP = new double[this.individualList.Count];
            for (int i = 0; i < this.individualList.Count; i++)
            {
                for (int j = this.partiallySharedRegion[i, 0]; j <= this.partiallySharedRegion[i, 1] && j <= this.partiallySharedLastBlockID; j++)
                {
                    double nextFrequency = ((Block)(newBlockDictionary.blockList[j])).GetFrequency((int)(this.partiallySharedBlockHaplotype[j - this.partiallySharedFirstBlockID]));
                    nextFrequency = 0.0 - Math.Log10(nextFrequency);
                    singleP[i] += nextFrequency;
                    if (j != this.partiallySharedRegion[i, 1] && j != this.partiallySharedLastBlockID)
                    {
                        ArrayList recombinationSnpList = new ArrayList();
                        for (int k = ((Block)(newBlockDictionary.blockList[j])).GetEndSnpID() - ((Block)(newBlockDictionary.blockList[this.partiallySharedFirstBlockID])).GetStartSnpID() + 1; k <= ((Block)(newBlockDictionary.blockList[j + 1])).GetStartSnpID() - ((Block)(newBlockDictionary.blockList[this.partiallySharedFirstBlockID])).GetStartSnpID() - 1; k++)
                            recombinationSnpList.Add(this.partiallySharedSnpHaplotype[k]);

                        nextFrequency = Math.Log10(this.GetRecombinationAValue(((Block)(newBlockDictionary.blockList[j])).GetEndSnpID() + 1, recombinationSnpList, newLinkageData));

                    }
                }
            }
            for (int i = 0; i < this.individualList.Count; i++)
            {
                this.AValue += singleP[i];
                if (i == 10)
                    break;
            }
            double sharedRegionAValue = 0;
            for (int i = allSharedFirstBlockID; i <= allSharedLastBlockID; i++)
            {
                sharedRegionAValue -= Math.Log10(((Block)(newBlockDictionary.blockList[i])).GetFrequency((int)haplotypeList[i - allSharedFirstBlockID]));
                if (i != allSharedLastBlockID)
                {
                    ArrayList recombinationSnpList = new ArrayList();
                    for (int j = ((Block)(newBlockDictionary.blockList[i])).GetEndSnpID() - ((Block)(newBlockDictionary.blockList[this.allSharedFirstBlockID])).GetStartSnpID() + 1; j <= ((Block)(newBlockDictionary.blockList[i + 1])).GetStartSnpID() - ((Block)(newBlockDictionary.blockList[this.allSharedFirstBlockID])).GetStartSnpID() - 1; j++)
                        recombinationSnpList.Add(this.allSharedSnpHaplotype[j]);
                    
                    if (sharedRegionAValue > 1000000)
                    {
                        Console.Write("here");
                    }
                }
                if (i == 10)
                    break;
            }

            double h0 = 1;
            for (int newBlockIndex = this.partiallySharedFirstBlockID; newBlockIndex <= this.partiallySharedLastBlockID; newBlockIndex++)
            {
                int shareCount = 0;
                for (int i = 0; i < this.individualList.Count; i++)
                {
                    if (this.partiallySharedRegion[i, 0] <= newBlockIndex && this.partiallySharedRegion[i, 1] >= newBlockIndex)
                        shareCount++;
                }

                h0 += Math.Log10(Math.Pow(2, shareCount)) -CombinationLog(shareCount, shareCount * 2);

            }
            h0 = h0 - this.AValue;

            int generation;
            double length;
            length = (((Block)(newBlockDictionary.blockList[this.partiallySharedLastBlockID])).GetEnd() - ((Block)(newBlockDictionary.blockList[this.partiallySharedFirstBlockID])).GetStart()) / 1000000.0;
            generation =Convert.ToInt32( 1.0 / length*100);
            double h1 = Math.Log10(0.5) * this.individualList.Count * generation;
            this.AValue = h1 - h0;
            
        }
        public double CombinationLog(int select, int all)
        {
            if (select < 1 || select >= all || all < 1)
                return 0;
            double result = 0;
            for (int i = 1; i <= all; i++)
                result += Math.Log10(i);
            for (int i = 1; i <= select; i++)
                result -= Math.Log10( i);
            for (int i = 1; i <= all - select; i++)
                result -= Math.Log10( i);
            
            return result;
        }
        public double GetRecombinationAValue(int startSnpIndex, ArrayList snpHaplotypeList, double[,] linkageData)
        {
            double recombinationP = 1;
            for (int snpIndex = startSnpIndex; snpIndex < startSnpIndex + snpHaplotypeList.Count - 1; snpIndex++)
            {
                if (snpIndex >= linkageData.GetUpperBound(0))
                    return recombinationP;
                if (((int)snpHaplotypeList[snpIndex - startSnpIndex]) == 0 && ((int)snpHaplotypeList[snpIndex + 1 - startSnpIndex]) == 0)
                    recombinationP *= linkageData[snpIndex, 2];
                else if (((int)snpHaplotypeList[snpIndex - startSnpIndex]) == 0 && ((int)snpHaplotypeList[snpIndex + 1 - startSnpIndex]) == 1)
                    recombinationP *= linkageData[snpIndex, 3];
                else if (((int)snpHaplotypeList[snpIndex - startSnpIndex]) == 1 && ((int)snpHaplotypeList[snpIndex + 1 - startSnpIndex]) == 0)
                    recombinationP *= linkageData[snpIndex, 4];
                else
                    recombinationP *= linkageData[snpIndex, 5];
            }
            return recombinationP;
        }
        public void SetPartiallySharedRegion(ArrayList[] IBDTable, BlockDictionary newBlockDictionary, ArrayList wholeIndividualList, int windowStart, int windowEnd)
        {
            this.partiallySharedRegion = new int[this.individualList.Count, 2];
            this.partiallySharedBlockHaplotype = new ArrayList();
            this.partiallySharedFirstBlockID = windowStart;
            this.partiallySharedLastBlockID = windowEnd - 1;
        //set the start point of the partially shared region.
            ArrayList lastRemainingIndividual = (ArrayList)this.individualList.Clone();
            for (int blockID = allSharedFirstBlockID - 1; blockID >= windowStart; blockID--)
            {
                int bestSharedHaplotype = -1;
                ArrayList remainningIndividual = new ArrayList();
                foreach (HaplotypeSharing newHaplotypeSharing in IBDTable[blockID]) //find the best sharing haplotype in this block.
                {
                    ArrayList intersection = this.SimpleIntersection(lastRemainingIndividual, newHaplotypeSharing.individualList);
                    
                    if (intersection.Count > remainningIndividual.Count)
                    {
                        remainningIndividual = (ArrayList)intersection.Clone();
                        bestSharedHaplotype = newHaplotypeSharing.haplotpeID;
                        if (remainningIndividual.Count == lastRemainingIndividual.Count)
                        {
                            //if (remainningIndividual.Count == this.individualList.Count)
                                //this.allSharedFirstBlockID = blockID;
                            break;
                        }
                    }
                }
                /*
                for (int recombinationSnpIndex = ((Block)(newBlockDictionary.blockList[blockID + 1])).GetStartSnpID() - 1; recombinationSnpIndex > ((Block)(newBlockDictionary.blockList[blockID])).GetEndSnpID(); recombinationSnpIndex--)
                {
                    ArrayList commonAlleleList = new ArrayList();
                    ArrayList rareAlleleList = new ArrayList();
                    for (int remainingIndex = 0; remainingIndex < remainningIndividual.Count; remainingIndex++)
                    {
                        if (((Individual)(wholeIndividualList[(int)remainningIndividual[remainingIndex]])).genotype[recombinationSnpIndex] == 0)
                            rareAlleleList.Add(remainingIndex);
                        else if (((Individual)(wholeIndividualList[(int)remainningIndividual[remainingIndex]])).genotype[recombinationSnpIndex] == 2)
                            commonAlleleList.Add(remainingIndex);
                        ArrayList removeList = new ArrayList();
                        if (rareAlleleList.Count >= commonAlleleList.Count)
                            removeList = (ArrayList)commonAlleleList.Clone();
                        else
                            removeList = (ArrayList)rareAlleleList.Clone();
                        removeList.Sort();
                        removeList.Reverse();
                        for (int removeIndex = 0; removeIndex < removeList.Count; removeIndex++)
                            remainningIndividual.RemoveAt((int)removeList[removeIndex]);
                    }
                }
                 */
                if (remainningIndividual.Count > 1)
                {
                    partiallySharedBlockHaplotype.Insert(0, bestSharedHaplotype);
                    if (remainningIndividual.Count == lastRemainingIndividual.Count)
                    {
                        continue;
                    }
                    else
                    {
                        ArrayList difference = this.GetDifferenceSet(lastRemainingIndividual, remainningIndividual);
                        foreach (int i in difference)
                        {
                            this.partiallySharedRegion[this.GetIndividualOrderInList(i), 0] = blockID + 1;
                        }
                        lastRemainingIndividual = (ArrayList)remainningIndividual.Clone();
                    }
                }
                else
                {
                    foreach (int i in lastRemainingIndividual)
                        this.partiallySharedRegion[this.GetIndividualOrderInList(i), 0] = blockID + 1;
                    this.partiallySharedFirstBlockID = blockID + 1;
                    //this.firstBlockID = blockID + 1;
                    break;
                }
            }
            for (int i = 0; i < this.individualList.Count; i++)
            {
                if (partiallySharedRegion[i, 0] == 0)
                    partiallySharedRegion[i, 0] = windowStart;
            }
            for (int i = 0; i < this.haplotypeList.Count; i++)
            {
                partiallySharedBlockHaplotype.Add(this.haplotypeList[i]);
            }
     //set the end point of the partially shared region.
            lastRemainingIndividual = (ArrayList)this.individualList.Clone();
            for (int blockID = allSharedLastBlockID + 1; blockID < windowEnd; blockID++)
            {
                int bestSharedHaplotype = -1;
                ArrayList remainningIndividual = new ArrayList();
                foreach (HaplotypeSharing newHaplotypeSharing in IBDTable[blockID]) //find the best sharing haplotype in this block.
                {
                    ArrayList intersection = this.SimpleIntersection(lastRemainingIndividual, newHaplotypeSharing.individualList);

                    if (intersection.Count > remainningIndividual.Count)
                    {
                        remainningIndividual = (ArrayList)intersection.Clone();
                        bestSharedHaplotype = newHaplotypeSharing.haplotpeID;
                        if (remainningIndividual.Count == lastRemainingIndividual.Count)
                        {
                            //if (remainningIndividual.Count == this.individualList.Count)
                                //this.allSharedLastBlockID = blockID;
                            break;
                        }
                    }
                }
                if (remainningIndividual.Count > 1)
                {
                    partiallySharedBlockHaplotype.Add(bestSharedHaplotype);
                    if (remainningIndividual.Count == lastRemainingIndividual.Count)
                    {
                        continue;
                    }
                    else
                    {
                        ArrayList difference = this.GetDifferenceSet(lastRemainingIndividual, remainningIndividual);
                        foreach (int i in difference)
                        {
                            this.partiallySharedRegion[this.GetIndividualOrderInList(i), 1] = blockID - 1;
                        }
                        lastRemainingIndividual = (ArrayList)remainningIndividual.Clone();
                    }
                }
                else
                {
                    foreach (int i in lastRemainingIndividual)
                        this.partiallySharedRegion[this.GetIndividualOrderInList(i), 1] = blockID - 1;
                    this.partiallySharedLastBlockID = blockID - 1;
                    //this.lastBlockID = blockID - 1;
                    break;
                }
            }
            for (int i = 0; i < this.individualList.Count; i++)
            {
                if (partiallySharedRegion[i, 1] == 0)
                    partiallySharedRegion[i, 1] = windowEnd;
            }
/*
            int remainIndividual = 0;
            for (int i = 0; i < this.individualList.Count; i++)
            {
                if (((Block)(newBlockDictionary.blockList[partiallySharedRegion[i, 1]])).GetEndSnpID() - ((Block)(newBlockDictionary.blockList[partiallySharedRegion[i, 0]])).GetStartSnpID() <= 10)
                {
                    this.partiallySharedRegion[i, 0] = -1;
                    this.partiallySharedRegion[i, 1] = -1;
                }
                else
                    remainIndividual++;
            }
            ArrayList newIndividualList = new ArrayList();
            int[,] newPartiallySharedRegion = new int[remainIndividual, 2];
            int individualOrder = 0;
            for (int i = 0; i < this.individualList.Count; i++)
            {
                if (this.partiallySharedRegion[i, 0] >= 0)
                {
                    newIndividualList.Add(this.individualList[i]);
                    newPartiallySharedRegion[individualOrder, 0] = this.partiallySharedRegion[i, 0];
                    newPartiallySharedRegion[individualOrder, 1] = this.partiallySharedRegion[i, 1];
                    individualOrder++;
                }
            }
            this.individualList = (ArrayList)(newIndividualList.Clone());
            this.partiallySharedRegion = (int[,])newPartiallySharedRegion.Clone();
 */ 
            if (this.individualList.Count < 2)
                return;
            while (this.allSharedFirstBlockID > windowStart)
            {
                bool foundUpperBound = false;
                for (int i = 0; i < this.individualList.Count; i++)
                {
                    if (this.partiallySharedRegion[i, 0] >= this.allSharedFirstBlockID)
                    {
                        foundUpperBound = true;
                        break;
                    }
                }
                if (foundUpperBound)
                    break;
                else
                {
                    this.allSharedFirstBlockID--;
                    this.haplotypeList.Insert(0, this.partiallySharedBlockHaplotype[this.allSharedFirstBlockID - this.partiallySharedFirstBlockID]);
                    //break;
                }
            }

            while (this.allSharedLastBlockID < this.partiallySharedLastBlockID)
            {
                bool foundLowerBound = false;
                for (int i = 0; i < this.individualList.Count; i++)
                {
                    if (this.partiallySharedRegion[i, 1] <= this.allSharedLastBlockID)
                    {
                        foundLowerBound = true;
                        break;
                    }
                }
                if (foundLowerBound)
                    break;
                else
                {
                   // if (this.allSharedLastBlockID >= windowEnd)
                   //     break;
                    this.allSharedLastBlockID++;
                    this.haplotypeList.Add(this.partiallySharedBlockHaplotype[this.allSharedLastBlockID - this.partiallySharedFirstBlockID]);
                }
            }
            if (this.partiallySharedLastBlockID == 0)
                this.partiallySharedLastBlockID = this.allSharedLastBlockID;
        }
        public int GetIndividualOrderInList(int individualID) //get the position of a individual in this individual list.
        {
            for (int i = 0; i < this.individualList.Count; i++)
                if ((int)(this.individualList[i]) == individualID)
                    return i;
            return -1;
        }
        public ArrayList SimpleIntersection(ArrayList arrayListA, ArrayList arrayListB)

        {
            ArrayList intersection = new ArrayList();
            foreach (int i in arrayListA)
            {
                foreach (int[] newPair in arrayListB)
                {
                    if (newPair[0] == i || newPair[1] == i)
                    {
                        intersection.Add(i);
                        break;
                    }
                }
            }
            return intersection;
            /*
            ArrayList newArrayListB = new ArrayList();
            int count = 0;
            foreach (int[] newPair in arrayListB)
            {
                newArrayListB.Add(newPair[0]);
                newArrayListB.Add(newPair[1]);
                count++;
                if (count == 100)
                    newArrayListB.Sort();
            }
            for (int i = 1; i < newArrayListB.Count; i++)
            {
                if ((int)(newArrayListB[i - 1]) == (int)(newArrayListB[i]))
                {
                    newArrayListB.RemoveAt(i);
                    i--;
                }
            }
            ArrayList intersection = new ArrayList();
            for (int i = 0; i < arrayListA.Count; i++)
            {
                if (newArrayListB.Contains(arrayListA[i]))
                    intersection.Add(arrayListA[i]);
            }
            return intersection;
             */ 
        }
        public ArrayList GetDifferenceSet(ArrayList arrayListA, ArrayList arrayListB) //return the A-B elements in A.
        {
            ArrayList difference = new ArrayList();
            for (int i = 0; i < arrayListA.Count; i++)
            {
                if (!arrayListB.Contains(arrayListA[i]))
                    difference.Add(arrayListA[i]);
            }
            return difference;
        }
        public bool CheckRecombinationRegion(ArrayList newIndividualList, BlockDictionary newBlockDictionary)
        {
            if (allSharedLastBlockID == newBlockDictionary.blockList.Count - 1)
                return false;
            int firstSnpID = ((Block)(newBlockDictionary.blockList[allSharedLastBlockID])).GetEndSnpID() + 1;
            int lastSnpID = ((Block)(newBlockDictionary.blockList[allSharedLastBlockID + 1])).GetStartSnpID() - 1;
            for (int i = firstSnpID + 1; i <= lastSnpID; i++)
            {
                ArrayList homoCommonList = new ArrayList();
                ArrayList homoRareList = new ArrayList();
                for (int j = 0; j < this.sharePairList.Count; j++)
                {

                    if (((Individual)(newIndividualList[((int[])(this.sharePairList[j]))[0]])).genotype[i] == 0 || ((Individual)(newIndividualList[((int[])(this.sharePairList[j]))[1]])).genotype[i] == 0)
                        homoCommonList.Add(j);
                    if (((Individual)(newIndividualList[((int[])(this.sharePairList[j]))[0]])).genotype[i] == 2 || ((Individual)(newIndividualList[((int[])(this.sharePairList[j]))[1]])).genotype[i] == 2)
                        homoRareList.Add(j);
                }
                if (homoCommonList.Count < homoRareList.Count)
                {
                    homoCommonList.Sort();
                    
                    for (int k = homoCommonList.Count - 1; k >= 0; k--)
                        this.sharePairList.RemoveAt((int)homoCommonList[k]);
                }
                else
                {
                    homoRareList.Sort();
                    
                    for (int k = homoRareList.Count - 1; k >= 0; k--)
                        this.sharePairList.RemoveAt((int)homoRareList[k]);
                }
                if (this.sharePairList.Count < 3)
                    return false;
            }
            return true;
            
        }
        private ArrayList GetListIntersection(HaplotypeSharing newGroupShare)
        {
            ArrayList intersection = new ArrayList();
            for (int i = 0; i < this.individualList.Count; i++)
            {
                if (newGroupShare.individualList.Contains(this.individualList[i]))
                    intersection.Add(this.individualList[i]);
            }
            return intersection;
        }
        private bool CheckSame(GroupShare newGroupShare)
        {
            
            if (this.allSharedLastBlockID != newGroupShare.allSharedLastBlockID)
                return false;
            if (this.allSharedFirstBlockID != newGroupShare.allSharedFirstBlockID)
                return false;
            if (this.individualList.Count != newGroupShare.individualList.Count)
                return false;
            this.individualList.Sort();
            newGroupShare.individualList.Sort();
            for (int i = 0; i < this.individualList.Count; i++)
            {
                if (((int)(this.individualList[i])) != ((int)(newGroupShare.individualList[i])))
                    return false;
            }
            /*for (int i = 0; i < this.haplotypeList.Count; i++)
            {
                if (((int)this.haplotypeList[i]) != ((int)newGroupShare.haplotypeList[i]))
                    return false;
            }
             */
            return true;
        }
        public bool CheckContain(GroupShare newGroupShare)
        {
            if (this.allSharedLastBlockID < newGroupShare.allSharedLastBlockID)
                return false;
            if (this.allSharedFirstBlockID > newGroupShare.allSharedFirstBlockID)
                return false;
            if (this.individualList.Count != newGroupShare.individualList.Count)
                return false;
            newGroupShare.individualList.Sort();
            this.individualList.Sort();
            for (int i = 0; i < this.individualList.Count; i++)
            {
                if (((int)(this.individualList[i])) != ((int)(newGroupShare.individualList[i])))
                    return false;
            }

            return true;
        }
        public bool ChekcSubGroup(GroupShare subGroup)
        {
            //this.SetIndividualList();
            //subGroup.SetIndividualList();
            this.individualList.Sort();
            subGroup.individualList.Sort();
            if (this.individualList.Count < subGroup.individualList.Count)
                return false;
            int thisIndex = 0;
            for (int subIndex = 0; subIndex < subGroup.individualList.Count; subIndex++)
            {
                while (thisIndex < this.individualList.Count)
                {
                    if ((int)(this.individualList[thisIndex]) < (int)(subGroup.individualList[subIndex]))
                        thisIndex++;
                    else if ((int)(this.individualList[thisIndex]) == (int)(subGroup.individualList[subIndex]))
                    {
                        thisIndex++;
                        break;
                    }
                    else
                        return false;
                    if (thisIndex == individualList.Count)
                        return false;
                }
            }
            return true;
        }
        public bool CheckOverlap(GroupShare newGroupShare) //combine same group of individuals that share close regions.
        {
            if (this.individualList.Count != newGroupShare.individualList.Count)
                return false;
            for (int i = 0; i < this.individualList.Count; i++)
            {
                if (((int)(this.individualList[i])) != ((int)(newGroupShare.individualList[i])))
                    return false;
            }
            if (this.allSharedLastBlockID - newGroupShare.allSharedFirstBlockID >= -2)
            {
                //if (this.lastBlockID < newGroupShare.firstBlockID)
                   // this.haplotypeList.Add(0);
                if (this.allSharedLastBlockID < newGroupShare.allSharedFirstBlockID)
                {
                    if (this.allSharedLastBlockID - newGroupShare.allSharedFirstBlockID == -2)
                        this.haplotypeList.Add(0);
                    for (int i = 0; i <= newGroupShare.allSharedLastBlockID - newGroupShare.allSharedFirstBlockID; i++)
                        this.haplotypeList.Add(newGroupShare.haplotypeList[i]);
                }
                else
                {
                    for (int i = this.allSharedLastBlockID - newGroupShare.allSharedFirstBlockID + 1; i <= newGroupShare.allSharedLastBlockID - newGroupShare.allSharedFirstBlockID; i++)
                        this.haplotypeList.Add(newGroupShare.haplotypeList[i]);
                }
                this.allSharedLastBlockID = newGroupShare.allSharedLastBlockID;
                return true;
            }
            else
                return false;
        }
        public void SetIndividualList() //generate the individual list from pair list
        {
            this.individualList = new ArrayList();
            foreach (int[] newPair in sharePairList)
            {
                if (individualList.IndexOf(newPair[0]) == -1)
                    this.individualList.Add(newPair[0]);
                if (individualList.IndexOf(newPair[1]) == -1)
                    this.individualList.Add(newPair[1]);

            }
            this.individualList.Sort();

        }

        public bool GroupShareExtention(ArrayList[] IBDTable, ArrayList newIndividualList, BlockDictionary newBlockDictionary, ArrayList newlyComingSharingList, List<GroupShare> finalResultList, ArrayList newlyAddedList)
        {
            if (this.allSharedLastBlockID >= IBDTable.GetUpperBound(0)) //the end of the block list, add to result if still some sharing left
            {
                if (this.individualList.Count > 1 && (this.allSharedLastBlockID - this.allSharedFirstBlockID >= 2)) //add to result list
                {
                    finalResultList.Add((GroupShare)this.MemberwiseClone());
                }
            }
            bool foundNext = false; //if find haplotype that can continue all the individuals in the next block, then "this" will not be added to list
            ArrayList passList = new ArrayList();
            ArrayList waitingToAddGroupList = new ArrayList();
            for (int haplotypeIndex = 0; haplotypeIndex < newlyComingSharingList.Count; haplotypeIndex++)
            {
                bool foundBefore = false;
                ArrayList pairListIntersection;
                pairListIntersection = this.GetPairListIntersection(this.sharePairList, ((GroupShare)(newlyComingSharingList[haplotypeIndex])).sharePairList);
                if (pairListIntersection.Count == this.sharePairList.Count)
                {
                    foundNext = true;   
                    this.allSharedLastBlockID += 1;
                    this.haplotypeList.Add(((GroupShare)(newlyComingSharingList[haplotypeIndex])).haplotypeList[0]);
                    if (pairListIntersection.Count == ((GroupShare)(newlyComingSharingList[haplotypeIndex])).sharePairList.Count)
                    {    
                        newlyComingSharingList.RemoveAt(haplotypeIndex);
                        haplotypeIndex--;
                    }
                    return true;
                }
                foreach (ArrayList eachList in passList) //check this intersection has already go checking next block, so do not need to do again for this subgroup.
                {
                    if (eachList.Count != pairListIntersection.Count)
                        continue;
                    else
                    {
                        for (int i = 0; i < pairListIntersection.Count; i++)
                        {
                            if (((int[])pairListIntersection[i])[0] != ((int[])(eachList[i]))[0] || ((int[])pairListIntersection[i])[1] != ((int[])(eachList[i]))[1])
                                break;
                            if (i == (pairListIntersection.Count - 1) && ((int[])pairListIntersection[i])[0] == ((int[])(eachList[i]))[0] && ((int[])pairListIntersection[i])[1] == ((int[])(eachList[i]))[1])
                                foundBefore = true;
                        }
                    }
                }
                if (foundBefore == true)  //checked before, do not need to go checking
                    break;
                if (pairListIntersection.Count > 1)
                {
                    passList.Add(pairListIntersection.Clone());
                    GroupShare newGroupShare = new GroupShare();
                    newGroupShare.sharePairList = (ArrayList)pairListIntersection.Clone();
                    newGroupShare.allSharedLastBlockID = this.allSharedLastBlockID + 1;
                    newGroupShare.allSharedFirstBlockID = this.allSharedFirstBlockID;
                    newGroupShare.haplotypeList = (ArrayList)this.haplotypeList.Clone();
                    newGroupShare.haplotypeList.Add(((GroupShare)(newlyComingSharingList[haplotypeIndex])).haplotypeList[0]);
                    newGroupShare.SetIndividualList();
                    if (newGroupShare.CheckRecombinationRegion(newIndividualList, newBlockDictionary)) // need to check the recombination region first.
                        waitingToAddGroupList.Add(newGroupShare);
                }
                if (foundNext == false && (this.allSharedLastBlockID - this.allSharedFirstBlockID >= 2))
                {
                    finalResultList.Add((GroupShare)(this.MemberwiseClone()));
                }
            }
            MyPrivateResultComparer newPrivateResultComparer = new MyPrivateResultComparer();
            waitingToAddGroupList.Sort(newPrivateResultComparer);
            for (int i = 0; i < waitingToAddGroupList.Count - 1; i++)
            {
                for (int j = i + 1; j < waitingToAddGroupList.Count; j++)
                {
                    if (((GroupShare)(waitingToAddGroupList[i])).ChekcSubGroup(((GroupShare)(waitingToAddGroupList[j]))) == true)
                    {
                        waitingToAddGroupList.RemoveAt(j);
                        j--;
                    }

                }
            }
            foreach (GroupShare newGroupShare in waitingToAddGroupList)
            {
                if (newlyAddedList.Count == 0)
                    newlyAddedList.Add(newGroupShare);
                else
                    for (int i = 0; i < newlyAddedList.Count; i++)
                    {
                        if (((GroupShare)(newlyAddedList[i])).ChekcSubGroup(newGroupShare) == true)
                            break;
                        if (newGroupShare.ChekcSubGroup((GroupShare)(newlyAddedList[i])) == true)
                        {
                            newlyAddedList[i] = newGroupShare;
                            break;
                        }
                        if (i == newlyAddedList.Count - 1)
                            newlyAddedList.Add(newGroupShare);
                    }
            }
            return false;
        }
        
        private ArrayList GetPairListIntersection(ArrayList firstList, ArrayList secondList)
        {
            MyPairComparer newMyPairComparer = new MyPairComparer();
            firstList.Sort(newMyPairComparer);
            secondList.Sort(newMyPairComparer);
            int firstIndex = 0;
            int secondIndex = 0;
            ArrayList intersection = new ArrayList();
            while (firstIndex < firstList.Count)
            {
                while (secondIndex < secondList.Count)
                {
                    if (((int[])(firstList[firstIndex]))[0] * 10000 + ((int[])(firstList[firstIndex]))[1] == ((int[])(secondList[secondIndex]))[0] * 10000 + ((int[])(secondList[secondIndex]))[1])
                    {
                        intersection.Add(firstList[firstIndex]);
                        secondIndex++;
                        firstIndex++;
                        break;
                    }
                    else if (((int[])(firstList[firstIndex]))[0] * 10000 + ((int[])(firstList[firstIndex]))[1] > ((int[])(secondList[secondIndex]))[0] * 10000 + ((int[])(secondList[secondIndex]))[1])
                    {
                        secondIndex++;
                    }
                    else
                    {
                        firstIndex++;
                        break;
                    }
                }
                if (secondIndex == secondList.Count)
                    break;
            }
            return intersection;

        }
        private class MyPairComparer : IComparer
        {
            public int Compare(Object x, Object y)
            {
                if (((int[])x)[0] > ((int[])y)[0])
                    return 1;
                else if ((((int[])x)[0] < ((int[])y)[0]))
                    return -1;
                else
                {
                    if (((int[])x)[1] > ((int[])y)[1])
                        return 1;
                    else if (((int[])x)[1] < ((int[])y)[1])
                        return -1;
                    else
                        return 0;
                }
            }
        }
        private class MyPrivateResultComparer : IComparer
        {
            public int Compare(Object x, Object y)
            {
                if (((GroupShare)y).individualList.Count != ((GroupShare)x).individualList.Count)
                    return (((GroupShare)x).individualList.Count - ((GroupShare)y).individualList.Count);
                for (int i = 0; i < ((GroupShare)x).individualList.Count; i++)
                {
                    if ((((int)((GroupShare)y).individualList[i]) - ((int)((GroupShare)x).individualList[i])) != 0)
                        return ((int)(((GroupShare)x).individualList[i]) - (int)(((GroupShare)y).individualList[i]));
                }

                if (((GroupShare)y).allSharedLastBlockID != ((GroupShare)x).allSharedLastBlockID)
                    return ((((GroupShare)x).allSharedLastBlockID - ((GroupShare)y).allSharedLastBlockID));
                if (((GroupShare)y).allSharedFirstBlockID != ((GroupShare)x).allSharedFirstBlockID)
                    return ((((GroupShare)x).allSharedFirstBlockID - ((GroupShare)y).allSharedFirstBlockID));

                return 0;

            }
        }
        
    }
}
