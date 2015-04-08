using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;
using System.IO;

namespace IBD
{

    class IBDList
    {
        public ArrayList pairwiseRegionList; //pairwise IBD region list.
        public ArrayList groupRegionList; // group IBD region list, generated from pairwise IBD list.
        public IBDList(ArrayList newPairwiseRegionList)
        {
            this.pairwiseRegionList = newPairwiseRegionList;
        }
        /*public void findGroupIBD(ArrayList newIndividualList, BlockDictionary newBlockDictionary)
        { 
            ArrayList[] IBDTable = new ArrayList[newBlockDictionary.blockList.Count];
            for (int regionIndex = 0; regionIndex < this.pairwiseRegionList.Count; regionIndex++) //fill in the table
            {
                for (int blockID = ((PairwiseIBDRegion)this.pairwiseRegionList[regionIndex]).GetFirstBlockID(); blockID <= ((PairwiseIBDRegion)this.pairwiseRegionList[regionIndex]).GetLastBlockID(); blockID++)
                {
                    foreach (int haplotypeID in ((PairwiseIBDRegion)this.pairwiseRegionList[regionIndex]).GetHaplotype(blockID - ((PairwiseIBDRegion)this.pairwiseRegionList[regionIndex]).GetFirstBlockID()))
                    {
                        if (IBDTable[blockID] == null)
                        {
                            HaplotypeSharing newHaplotypeShare = new HaplotypeSharing(haplotypeID, ((PairwiseIBDRegion)this.pairwiseRegionList[regionIndex]).GetFirstIndividual());
                            newHaplotypeShare.individualList.Add(((PairwiseIBDRegion)this.pairwiseRegionList[regionIndex]).GetSecondIndividual());
                            IBDTable[blockID] = new ArrayList();
                            IBDTable[blockID].Add(newHaplotypeShare);
                        }
                        else
                        {
                            bool added = false;
                            for (int newHaplotypeShare = 0; newHaplotypeShare < IBDTable[blockID].Count; newHaplotypeShare++)
                            {
                                if (((HaplotypeSharing)IBDTable[blockID][newHaplotypeShare]).haplotpeID == haplotypeID)
                                {
                                    ((HaplotypeSharing)IBDTable[blockID][newHaplotypeShare]).individualList.Add(((PairwiseIBDRegion)this.pairwiseRegionList[regionIndex]).GetFirstIndividual());
                                    ((HaplotypeSharing)IBDTable[blockID][newHaplotypeShare]).individualList.Add(((PairwiseIBDRegion)this.pairwiseRegionList[regionIndex]).GetSecondIndividual());
                                    added = true;
                                }
                            }
                            if (!added)
                            {
                                HaplotypeSharing newHaplotypeShare = new HaplotypeSharing(haplotypeID, ((PairwiseIBDRegion)this.pairwiseRegionList[regionIndex]).GetFirstIndividual());
                                newHaplotypeShare.individualList.Add(((PairwiseIBDRegion)this.pairwiseRegionList[regionIndex]).GetSecondIndividual());
                                IBDTable[blockID].Add(newHaplotypeShare);
                            }
                        }
                    }
                    
                }
            }
            for (int i = 0; i <= IBDTable.GetUpperBound(0); i++)
            {
                for (int j = 0; j < IBDTable[i].Count; j++)
                {
                    ((HaplotypeSharing)IBDTable[i][j]).individualList.Sort();
                    for (int k = 0; k < ((HaplotypeSharing)IBDTable[i][j]).individualList.Count - 1; k++)
                    {
                        if ((int)(((HaplotypeSharing)IBDTable[i][j]).individualList[k]) == (int)(((HaplotypeSharing)IBDTable[i][j]).individualList[k + 1]))
                        {
                            ((HaplotypeSharing)IBDTable[i][j]).individualList.RemoveAt(k);
                            k--;
                        }
                    }
                }
            }

 
            ArrayList resultList = new ArrayList();
            for (int blockIndex = 0; blockIndex < IBDTable.GetUpperBound(0); blockIndex++)
            {
                foreach (HaplotypeSharing newHaplotypeSharing in IBDTable[blockIndex])
                {
                    if (newHaplotypeSharing.individualList.Count > 3)
                    {
                        GroupShare newGroupShare = new GroupShare();
                        newGroupShare.individualList = (ArrayList)newHaplotypeSharing.individualList.Clone();
                        newGroupShare.firstBlockID = blockIndex;
                        newGroupShare.lastBlockID = blockIndex;
                        newGroupShare.haplotypeList = new ArrayList();
                        newGroupShare.haplotypeList.Add(newHaplotypeSharing.haplotpeID);
                        newGroupShare.Extend(IBDTable, resultList,newIndividualList, newBlockDictionary);
                    }
                }
                Console.WriteLine("{0} {1}", blockIndex, resultList.Count);
            }
            for (int i = 0; i < resultList.Count; i++)
            {
                ((GroupShare)resultList[i]).individualList.Sort();
                for (int j = 0; j < ((GroupShare)resultList[i]).individualList.Count - 1; j++)
                {
                    if ((int)(((GroupShare)resultList[i]).individualList[j]) == (int)(((GroupShare)resultList[i]).individualList[j + 1]))
                    {
                        ((GroupShare)resultList[i]).individualList.RemoveAt(j + 1);
                        j--;
                    }
                }
                if (((GroupShare)resultList[i]).individualList.Count <= 3)
                {
                    resultList.RemoveAt(i);
                    i--;
                }
            }
            MyResultComparer newResultComparer = new MyResultComparer();
            resultList.Sort(newResultComparer);
            for (int i = 1; i < resultList.Count; i++)
            {
                if (((GroupShare)resultList[i]).lastBlockID != ((GroupShare)resultList[i - 1]).lastBlockID)
                    continue;
                if (((GroupShare)resultList[i]).firstBlockID < ((GroupShare)resultList[i - 1]).firstBlockID)
                    continue;
                if (((GroupShare)resultList[i]).individualList.Count != ((GroupShare)resultList[i - 1]).individualList.Count)
                    continue;
                bool same = true;
                for (int j = 0; j < ((GroupShare)resultList[i]).individualList.Count; j++)
                {
                    if ((int)(((GroupShare)resultList[i]).individualList[j]) != (int)(((GroupShare)resultList[i - 1]).individualList[j]))
                    {
                        same = false;
                        break;
                    }
                }
                if (same)
                {
                    resultList.RemoveAt(i);
                    i--;
                }
            }
            FileStream cFile = new FileStream(@"F:\20cmtest\newresult.txt", FileMode.OpenOrCreate);
            StreamWriter sw1 = new StreamWriter(cFile);
            for (int i = 0; i < resultList.Count; i++)
            {
                if (((GroupShare)resultList[i]).individualList.Count < 4)
                {
                    resultList.RemoveAt(i);
                    i--;
                    continue;
                }
                sw1.Write("{0} {1} {2}\t", ((GroupShare)resultList[i]).firstBlockID, ((GroupShare)resultList[i]).lastBlockID, ((GroupShare)resultList[i]).individualList.Count);
                for (int j = 0; j < ((GroupShare)resultList[i]).individualList.Count; j++)
                    sw1.Write("{0} ", ((GroupShare)resultList[i]).individualList[j]);
                sw1.Write("\t");
                for (int k = 0; k < ((GroupShare)resultList[i]).haplotypeList.Count; k++)
                    sw1.Write("{0}", ((GroupShare)resultList[i]).haplotypeList[k]);
                sw1.WriteLine();
            }

            foreach (GroupShare newGroupShare in resultList)
            {
                //newGroupShare.SetPartiallySharedRegion(IBDTable,newBlockDictionary);
            }
            for (int i = 1; i < resultList.Count; i++)
            {
                if (((GroupShare)resultList[i]).lastBlockID != ((GroupShare)resultList[i - 1]).lastBlockID)
                    continue;
                if (((GroupShare)resultList[i]).firstBlockID < ((GroupShare)resultList[i - 1]).firstBlockID)
                    continue;
                if (((GroupShare)resultList[i]).individualList.Count != ((GroupShare)resultList[i - 1]).individualList.Count)
                    continue;
                bool same = true;
                for (int j = 0; j < ((GroupShare)resultList[i]).individualList.Count; j++)
                {
                    if ((int)(((GroupShare)resultList[i]).individualList[j]) != (int)(((GroupShare)resultList[i - 1]).individualList[j]))
                    {
                        same = false;
                        break;
                    }
                }
                if (same)
                {
                    resultList.RemoveAt(i);
                    i--;
                }
            }

            resultList.Sort(newResultComparer);
            foreach (GroupShare newGroupShare in resultList)

            {
                newGroupShare.SetAValue(newBlockDictionary);
            }
            MyResultAValueComparer newResultAValueComparer = new MyResultAValueComparer();
            resultList.Sort(newResultAValueComparer);

            for (int i = 0; i < resultList.Count; i++)
            {
                if (((GroupShare)resultList[i]).individualList.Count < 4)
                {
                    resultList.RemoveAt(i);
                    i--;
                    continue;
                }
                sw1.Write("{0}\t", ((GroupShare)resultList[i]).AValue);
                sw1.Write("{0} {1} {2}\t", ((GroupShare)resultList[i]).firstBlockID, ((GroupShare)resultList[i]).lastBlockID, ((GroupShare)resultList[i]).individualList.Count);
                for (int j = 0; j < ((GroupShare)resultList[i]).individualList.Count; j++)
                    sw1.Write("{0} ", ((GroupShare)resultList[i]).individualList[j]);
                sw1.Write("\t");
                for (int k = 0; k < ((GroupShare)resultList[i]).haplotypeList.Count; k++)
                    sw1.Write("{0}", ((GroupShare)resultList[i]).haplotypeList[k]);
                sw1.WriteLine();
            }
            sw1.Close();
                Console.ReadKey();
        }
         */

        private class MyResultAValueComparer : IComparer
        {
            public int Compare(Object x, Object y)
            {
                if (((GroupShare)y).AValue != ((GroupShare)x).AValue)
                    return ((int)(((GroupShare)y).AValue - ((GroupShare)x).AValue)*100000);
                else
                    return 0;

            }
        }
    }
}
