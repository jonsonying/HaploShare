using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;
using System.IO;

namespace IBD
{
    class Individual
    {
        private int individualID;
        public int[] genotype;
        public ArrayList[] haplotype;
        public string individualname;
        public Individual(int[] newGenotype, int newID)
        { 
            this.genotype = (int[])newGenotype.Clone();
            this.individualID = newID;
        }
        public int GetIndividualID()
        {
            return this.individualID;
        }
        public void SetIndividualID(int newId)
        {
            this.individualID = newId;
        }
        public void Haplotyping(BlockDictionary blockDictionary, MapData mapData, double[,] newLinkageData)
        {
            this.haplotype = new ArrayList[blockDictionary.blockList.Count];
            int blockID = -1;
            foreach (Block newBlock in blockDictionary.blockList)
            {
                blockID++;
                if (newBlock.GetEndSnpID() + newBlock.GetStartSnpID() == 0) continue;
                int[] blockGenotype = new int[newBlock.GetEndSnpID() - newBlock.GetStartSnpID() + 1];
                for (int snpID = newBlock.GetStartSnpID(); snpID <= newBlock.GetEndSnpID(); snpID++)
                {
                    blockGenotype[snpID - newBlock.GetStartSnpID()] = this.genotype[snpID];
                }
                ArrayList newHaplotypePairList = newBlock.GenotypeToHaplotype(blockGenotype, newLinkageData);
                this.haplotype[blockID] = newHaplotypePairList;

            }

        }
        private bool CombineAdjacentPairwiseRegion(PairwiseIBDRegion firstRegion, PairwiseIBDRegion secondRegion, Individual thatIndividual)
        {
            if (firstRegion.GetLastBlockID() + 2 < secondRegion.GetFirstBlockID())
                return false;
            else
            {
                int errorCount = 0;
                for (int i = firstRegion.GetEnd() + 1; i < secondRegion.GetStart(); i++)
                    if (Math.Abs(this.genotype[i] - thatIndividual.genotype[i]) == 2)
                        errorCount++;
                if (errorCount > 1)
                    return false;
                else
                {

                    if (firstRegion.GetLastBlockID() + 2 == secondRegion.GetFirstBlockID())
                    {
                        ArrayList errorBlockHaplotypeList = new ArrayList();
                        foreach (HaplotypePair newPair in ((ArrayList)(this.haplotype[firstRegion.GetLastBlockID() + 1])))
                        {
                            errorBlockHaplotypeList.Add(newPair.firstHaplotypeID);
                            errorBlockHaplotypeList.Add(newPair.secondHaplotypeID);
                        }
                        foreach (HaplotypePair newPair in ((ArrayList)(thatIndividual.haplotype[firstRegion.GetLastBlockID() + 1])))
                        {
                            errorBlockHaplotypeList.Add(newPair.firstHaplotypeID);
                            errorBlockHaplotypeList.Add(newPair.secondHaplotypeID);
                        }
                        errorBlockHaplotypeList.Sort();
                        for (int i = 1; i < errorBlockHaplotypeList.Count; i++)
                        {
                            if (((int)(errorBlockHaplotypeList[i])) == ((int)(errorBlockHaplotypeList[i - 1])))
                            {
                                errorBlockHaplotypeList.RemoveAt(i);
                                i--;
                            }
                        }
                        firstRegion.AddHaplotype(errorBlockHaplotypeList);
                    }
                    
                    for (int i = 0; i < secondRegion.GetHaplotypeLength(); i++)
                    {
                        firstRegion.AddHaplotype(secondRegion.GetHaplotype(i));
                    }
                    firstRegion.SetLastBlockID(secondRegion.GetLastBlockID());
                    firstRegion.SetEnd(secondRegion.GetEnd());
                    return true; // can be combined
                }
            }

        }
        public void PairwiseComparison(ArrayList pairwiseIBDRegionList, Individual thatIndividual, Block[] blockList, MapData newMapData, double lengthCut, int startBlockID, int endBlockID)
        {
            PairwiseIBDRegion newPairwiseIBDRegion = new PairwiseIBDRegion(this.GetIndividualID(), thatIndividual.GetIndividualID(), 0, 0, -1, -1);
            for (int blockID = startBlockID; blockID < endBlockID; blockID++) //foreach block, from first to end.
            {
                if (newPairwiseIBDRegion.GetLastBlockID() == -1) // last IBD region has just added or deleted, need to construc a new one.
                {
                    newPairwiseIBDRegion.SetStart(blockList[blockID].GetStartSnpID());
                    newPairwiseIBDRegion.SetEnd(blockList[blockID].GetEndSnpID());
                    newPairwiseIBDRegion.SetLastBlockID(blockID);
                    newPairwiseIBDRegion.SetFirstBlockID(blockID);
                    newPairwiseIBDRegion.ClearHaplotype();
                }
                ArrayList sharedHaplotypeID = new ArrayList();
                ArrayList thisHaplotypeList = new ArrayList();//store all the possible haplotype for this individual in this block.
                ArrayList thatHaplotypeList = new ArrayList();//store all the possible haplotype for that individual in this block.
                bool thisUnknownHaplotype = false;
                bool thatUnknownHaplotype = false;
                for (int i = 0; i < this.haplotype[blockID].Count; i++) //store all the possible haplotype for this individual in this block.
                {
                    if (((HaplotypePair)this.haplotype[blockID][i]).firstHaplotypeID != -1)
                        thisHaplotypeList.Add(((HaplotypePair)this.haplotype[blockID][i]).firstHaplotypeID);
                    else
                    {
                        thisUnknownHaplotype = true;
                        break;
                    }
                    if (((HaplotypePair)this.haplotype[blockID][i]).secondHaplotypeID != -1 && ((HaplotypePair)this.haplotype[blockID][i]).secondHaplotypeID != ((HaplotypePair)this.haplotype[blockID][i]).firstHaplotypeID)
                        thisHaplotypeList.Add(((HaplotypePair)this.haplotype[blockID][i]).secondHaplotypeID);
                    else
                    {
                        if (((HaplotypePair)this.haplotype[blockID][i]).secondHaplotypeID == -1)
                            thisUnknownHaplotype = true;
                        break;
                    }
                }
                for (int i = 0; i < thatIndividual.haplotype[blockID].Count; i++)//store all the possible haplotype for that individual in this block.
                {
                    if (((HaplotypePair)thatIndividual.haplotype[blockID][i]).firstHaplotypeID != -1)
                        thatHaplotypeList.Add(((HaplotypePair)thatIndividual.haplotype[blockID][i]).firstHaplotypeID);
                    else
                    {
                        thatUnknownHaplotype = true;
                        break;
                    }
                    if (((HaplotypePair)thatIndividual.haplotype[blockID][i]).secondHaplotypeID != -1 && ((HaplotypePair)thatIndividual.haplotype[blockID][i]).secondHaplotypeID != ((HaplotypePair)thatIndividual.haplotype[blockID][i]).firstHaplotypeID)
                        thatHaplotypeList.Add(((HaplotypePair)thatIndividual.haplotype[blockID][i]).secondHaplotypeID);
                    else
                    {
                        if (((HaplotypePair)thatIndividual.haplotype[blockID][i]).secondHaplotypeID == -1)
                            thatUnknownHaplotype = true;
                        break;
                    }
                }

                if (thisUnknownHaplotype && thatUnknownHaplotype) //found unknown haplotype in both sample 
                {
                    for (int s = 0; s < blockList[blockID].GetHaplotypeCount(); s++)
                        sharedHaplotypeID.Add(s);
                }
                else if (thisUnknownHaplotype && !thatUnknownHaplotype) // only in this sample
                { 
                    foreach (int newId in thatHaplotypeList)
                        sharedHaplotypeID.Add(newId);
                }
                else if((!thisUnknownHaplotype) && thatUnknownHaplotype) // only that sample
                {
                    foreach (int newId in thisHaplotypeList)
                        sharedHaplotypeID.Add(newId);
                }
                else                                                // normal case ,non unknown
                {
                    thisHaplotypeList.Sort();
                    thatHaplotypeList.Sort();

                    int thisHaplotypeIndex = 0;
                    int thatHaplotypeIndex = 0;
                    while (thisHaplotypeIndex < thisHaplotypeList.Count && thatHaplotypeIndex < thatHaplotypeList.Count)
                    { // store all the shared haplotype id in the sharedHaplotypeID list.
                        if ((int)thisHaplotypeList[thisHaplotypeIndex] == (int)thatHaplotypeList[thatHaplotypeIndex])
                        {
                            sharedHaplotypeID.Add(thisHaplotypeList[thisHaplotypeIndex]);
                            thisHaplotypeIndex++;
                            thatHaplotypeIndex++;
                        }
                        else if ((int)thisHaplotypeList[thisHaplotypeIndex] < (int)thatHaplotypeList[thatHaplotypeIndex])
                            thisHaplotypeIndex++;
                        else
                            thatHaplotypeIndex++;
                    }
                }
                bool shareRecombinationRegion = true;

                if (sharedHaplotypeID.Count == 0) // no shared haplotype found!
                {
                    if (newMapData.GetGeneticDistance(newPairwiseIBDRegion.GetEnd()) - newMapData.GetGeneticDistance(newPairwiseIBDRegion.GetStart()) >= lengthCut/10) // only add the pairwise IBD region whose length is >= 2cm
                    {
                        if (newPairwiseIBDRegion.GetLastBlockID() - newPairwiseIBDRegion.GetFirstBlockID() != 0)
                            if (pairwiseIBDRegionList.Count == 0 || !this.CombineAdjacentPairwiseRegion((PairwiseIBDRegion)pairwiseIBDRegionList[pairwiseIBDRegionList.Count - 1], newPairwiseIBDRegion, thatIndividual))
                                if (newMapData.GetGeneticDistance(newPairwiseIBDRegion.GetEnd()) - newMapData.GetGeneticDistance(newPairwiseIBDRegion.GetStart()) >= lengthCut)
                                    pairwiseIBDRegionList.Add(newPairwiseIBDRegion.GetCopy());
                    }
                    newPairwiseIBDRegion.SetLastBlockID(-1);
                    continue;
                }
                else
                {
                    newPairwiseIBDRegion.SetEnd(blockList[blockID].GetEndSnpID());
                    newPairwiseIBDRegion.AddHaplotype(sharedHaplotypeID);
                    newPairwiseIBDRegion.SetLastBlockID(blockID);
                    if (blockID == blockList.Length - 1)
                    {
                        if (newMapData.GetGeneticDistance(newPairwiseIBDRegion.GetEnd()) - newMapData.GetGeneticDistance(newPairwiseIBDRegion.GetStart()) >= lengthCut)
                            if (pairwiseIBDRegionList.Count == 0 || !this.CombineAdjacentPairwiseRegion((PairwiseIBDRegion)pairwiseIBDRegionList[pairwiseIBDRegionList.Count - 1], newPairwiseIBDRegion, thatIndividual))
                                if (newMapData.GetGeneticDistance(newPairwiseIBDRegion.GetEnd()) - newMapData.GetGeneticDistance(newPairwiseIBDRegion.GetStart()) >= lengthCut)
                                    pairwiseIBDRegionList.Add(newPairwiseIBDRegion.GetCopy());
                        return;
                    }
                    else
                        for (int i = blockList[blockID].GetEndSnpID() + 1; i < blockList[blockID + 1].GetStartSnpID(); i++)
                        {
                            if (Math.Abs(this.genotype[i] - thatIndividual.genotype[i]) == 2)
                            {
                                if (newMapData.GetGeneticDistance(newPairwiseIBDRegion.GetEnd()) - newMapData.GetGeneticDistance(newPairwiseIBDRegion.GetStart()) >= lengthCut) // only add the pairwise IBD region whose length is >= 2cm
                                    if (pairwiseIBDRegionList.Count == 0 || !this.CombineAdjacentPairwiseRegion((PairwiseIBDRegion)pairwiseIBDRegionList[pairwiseIBDRegionList.Count - 1], newPairwiseIBDRegion, thatIndividual))
                                        if (newMapData.GetGeneticDistance(newPairwiseIBDRegion.GetEnd()) - newMapData.GetGeneticDistance(newPairwiseIBDRegion.GetStart()) >= lengthCut)
                                            pairwiseIBDRegionList.Add(newPairwiseIBDRegion.GetCopy());
                                newPairwiseIBDRegion.SetLastBlockID(-1);
                                newPairwiseIBDRegion.SetEnd(i - 1);
                                break;
                            }
                        }
                }
            }
        }
    }
}
