using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;
using System.IO;

namespace IBD
{
    public struct HaplotypeSharing
    {
        public int haplotpeID;
        public ArrayList individualList; //list of pair of individuals: (1,2),(2,3),(4,5)...
        public HaplotypeSharing(int newHaplotypeID, ArrayList newPairList)
        {
            this.haplotpeID = newHaplotypeID;
            this.individualList = (ArrayList)(newPairList.Clone());
        }
        public HaplotypeSharing(int newHaplotypeID, int[] firstPair)
        {
            this.haplotpeID = newHaplotypeID;
            ArrayList newList = new ArrayList();
            newList.Add(firstPair);
            this.individualList = newList;
        }
    }
    class GroupShareFinder
    {
        public ArrayList pairwiseRegionList; //pairwise IBD region list.
        public ArrayList groupRegionList; // group IBD region list, generated from pairwise IBD list.
        public GroupShareFinder(ArrayList newPairwiseRegionList)
        {
            this.pairwiseRegionList = newPairwiseRegionList;
        }
        public void findGroupIBD(ArrayList newIndividualList, BlockDictionary newBlockDictionary, double[,] linkageData, List<GroupShare> result, List<GroupShare> resultList, int windowStart, int windowEnd)
        {
            ArrayList[] IBDTable = new ArrayList[newBlockDictionary.blockList.Count];
            for (int regionIndex = 0; regionIndex < pairwiseRegionList.Count; regionIndex++)//(PairwiseIBDRegion newRegion in this.pairwiseRegionList)
            {
                PairwiseIBDRegion newRegion = (PairwiseIBDRegion)(pairwiseRegionList[regionIndex]);
                for (int blockID = newRegion.GetFirstBlockID(); blockID <= newRegion.GetLastBlockID(); blockID++)
                {
                    foreach (int haplotypeID in newRegion.GetHaplotype(blockID - newRegion.GetFirstBlockID()))
                    {
                        if (IBDTable[blockID] == null)
                        {
                            int[] newPair = new int[2];
                            newPair[0] = newRegion.GetFirstIndividual();
                            newPair[1] = newRegion.GetSecondIndividual();
                            HaplotypeSharing newHaplotypeShare = new HaplotypeSharing(haplotypeID,newPair);
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
                                    int[] newPair = new int[2];
                                    newPair[0] = newRegion.GetFirstIndividual();
                                    newPair[1] = newRegion.GetSecondIndividual();
                                    ((HaplotypeSharing)IBDTable[blockID][newHaplotypeShare]).individualList.Add(newPair);
                                    added = true;
                                }
                            }
                            if (!added)
                            {
                                int[] newPair = new int[2];
                                newPair[0] = newRegion.GetFirstIndividual();
                                newPair[1] = newRegion.GetSecondIndividual();
                                HaplotypeSharing newHaplotypeShare = new HaplotypeSharing(haplotypeID, newPair);
                                IBDTable[blockID].Add(newHaplotypeShare);
                            }
                        }
                    }
                }
                pairwiseRegionList[regionIndex] = null;
            }

            pairwiseRegionList = null;
            GC.Collect();
            MyResultComparer newResultComparer = new MyResultComparer();

            if (IBDTable[IBDTable.GetUpperBound(0)] == null)
                IBDTable[IBDTable.GetUpperBound(0)] = new ArrayList();
            for (int blockIndex = windowStart; blockIndex < windowEnd; blockIndex++)
            {
                ArrayList newlyComingSharingList = new ArrayList();
                if (IBDTable[blockIndex] == null)
                {
                    foreach (GroupShare newGroupShare in resultList)
                        result.Add(newGroupShare);
                    resultList = new List<GroupShare>();
                    IBDTable[blockIndex] = new ArrayList();
                    continue;
                }
                foreach (HaplotypeSharing newHaplotypeSharing in IBDTable[blockIndex])
                {
                    if (newHaplotypeSharing.individualList.Count > 2 && newHaplotypeSharing.individualList.Count < 0.2*newIndividualList.Count) //at least 2 pairs sharing this haplotype.
                    {
                        GroupShare newGroupShare = new GroupShare();
                        newGroupShare.sharePairList = (ArrayList)newHaplotypeSharing.individualList.Clone();
                        newGroupShare.allSharedFirstBlockID = blockIndex;
                        newGroupShare.allSharedLastBlockID = blockIndex;
                        newGroupShare.haplotypeList = new ArrayList();
                        newGroupShare.haplotypeList.Add(newHaplotypeSharing.haplotpeID);
                        newGroupShare.SetIndividualList();
                        newlyComingSharingList.Add(newGroupShare);

                    }
                }
                ArrayList newlyAddedList = new ArrayList();


                //for the next block, if no haplotype contains more than 1 pair, stop the extension and record all the region waiting to be extended.
                if (newlyComingSharingList.Count == 0)
                {
                    foreach (GroupShare newGroup in resultList)
                        result.Add(newGroup);
                } 
                
                for (int i = 0; i < resultList.Count; i++)
                {
                    GroupShare newGroupshare1 = (GroupShare)(resultList[i]);
                    if (newGroupshare1.GroupShareExtention(IBDTable, newIndividualList, newBlockDictionary, newlyComingSharingList,result,newlyAddedList)== false)
                    {
                        resultList.RemoveAt(i);
                        i--;
                    }
                }
                foreach (GroupShare newGroupShare in newlyAddedList)
                {
                    bool found = false;
                    foreach (GroupShare checkingGroup in resultList)
                    {
                        if (checkingGroup.ChekcSubGroup(newGroupShare) == true)
                        {
                            found = true;
                            break;
                        }
                    }
                    if (found == false)
                        resultList.Add(newGroupShare);
                }
                foreach (GroupShare newGroupShare in newlyComingSharingList)
                {
                    bool found = false;
                    foreach (GroupShare checkingGroup in resultList)
                    {
                        
                        if (checkingGroup.ChekcSubGroup(newGroupShare) == true)
                        {
                            found = true;
                            break;
                        }
                    }
                    if (found == false)
                        resultList.Add(newGroupShare);
                }
                
                resultList.Sort(newResultComparer);

           }


            foreach (GroupShare newShare in result)
            {
                newShare.SetIndividualList();
            }


            result.Sort(newResultComparer);

            for (int i = 0; i < result.Count - 1; i++)
            {
                if (((GroupShare)(result[i])).CheckContain((GroupShare)(result[i + 1])) == true)
                {
                    result.RemoveAt(i + 1);
                    i--;
                }
                else if (((GroupShare)(result[i + 1])).CheckContain(((GroupShare)(result[i]))) == true)
                {
                    result.RemoveAt(i);
                    i--;
                }
            }


            for (int i = 0; i < result.Count - 1; i++)
            {
                if (((GroupShare)(result[i])).CheckOverlap((GroupShare)(result[i + 1])))
                {
                    result.RemoveAt(i + 1);
                    i--;
                }
            }

            foreach (GroupShare newShare in result)
            {
                newShare.SetIndividualList();
            }

            for (int i = 0; i < result.Count; i++)
            {
                ((GroupShare)(result[i])).SetPartiallySharedRegion(IBDTable, newBlockDictionary, newIndividualList, windowStart, windowEnd);
                Console.SetCursorPosition(16, Console.CursorTop);
                Console.Write("{0} / {1}", i +1, result.Count);

            }

            for (int i = 0; i < result.Count; i++)
            {
                ((GroupShare)(result[i])).SetAllSharedSnpHaplotype(newIndividualList,newBlockDictionary);
            }
            for (int i = 0; i < result.Count; i++)
            {
                ((GroupShare)(result[i])).SetPartiallySharedSnpHaplotype(newIndividualList, newBlockDictionary);
            }

            for (int i = 0; i < result.Count; i++)
            {
                ((GroupShare)(result[i])).SetAValue(newBlockDictionary,newIndividualList,linkageData);
            }
            Console.WriteLine("\t\t\t\t\tFinished.");
            result.Sort(newResultComparer);

            for (int i = 0; i < result.Count - 1; i++)
            {
                if (((GroupShare)(result[i])).CheckContain((GroupShare)(result[i + 1])))
                {
                    result.RemoveAt(i + 1);
                    i--;
                }
            }
            for (int i = 0; i < result.Count - 1; i++)
            {
                if (((GroupShare)(result[i])).CheckOverlap((GroupShare)(result[i + 1])))
                {
                    result.RemoveAt(i + 1);
                    i--;
                }
            }

 
        }

        private class MyResultAValueComparer : IComparer
        {
            public int Compare(Object x, Object y)
            {
                if (((GroupShare)y).AValue != ((GroupShare)x).AValue)
                    return ((int)((((GroupShare)y).AValue - ((GroupShare)x).AValue) * 100000));
                else
                    return 0;

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
