using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;
using System.IO;

namespace IBD
{
    class ResultList
    {
        public List<GroupShare> result;//case result
        public List<List<GroupShare>> controlResult;// list[2] stands for data for 3 individuals, same below.
        public List<List<double>> controlResultParameters;
        public double[,] pTable;

        public void SetpTable()
        { 
            pTable = new double[801,2];
            char[] splitArray = new char[] { ' ' };
            FileStream aFile = new FileStream("distribution.txt", FileMode.Open);
            StreamReader sr = new StreamReader(aFile);
            for (int i = 0; i < 801; i++)
            {
                string strLine = sr.ReadLine();
                string[] strArray = strLine.Split(splitArray);
                pTable[i, 0] = Convert.ToDouble(strArray[0]);
                pTable[i, 1] = Convert.ToDouble(strArray[1]);
            }

        }
        public void SetControlResult (int totalCase, double maxIndRate) //set up the control result list from permutation
        {
            controlResult = new List<List<GroupShare>>(Convert.ToInt32(Math.Round(totalCase * maxIndRate)));
            for (int i = 0; i < Convert.ToInt32(Math.Round(totalCase * maxIndRate)); i++)
            {
                controlResult.Add(new List<GroupShare>());
            }
            
        }
        public void AddNewControlList(List<GroupShare> newList, int totalCase, double maxIndRate) //add the permutation result to the list for further analysis
        {
            List<GroupShare> currentList = new List<GroupShare>(Convert.ToInt32(Math.Round(totalCase * maxIndRate)));
            for (int i = 0; i < Convert.ToInt32(Math.Round(totalCase * maxIndRate)); i++)
            {
                currentList.Add(null);
            }
            foreach (GroupShare newGroupShare in newList)
            {
                if (newGroupShare.individualList.Count > totalCase * maxIndRate)
                    break;
                if (currentList[newGroupShare.individualList.Count - 1] == null || currentList[newGroupShare.individualList.Count - 1].AValue < newGroupShare.AValue)
                {
                    bool flag = true;
                    foreach (GroupShare preGroupShare in this.controlResult[newGroupShare.individualList.Count - 1])
                    {
                        if (this.CheckDuplication(preGroupShare, newGroupShare))
                        {
                            flag = false;
                            break;
                        }
                            
                    }
                    if (flag)
                        currentList[newGroupShare.individualList.Count - 1] = newGroupShare;
                }

            }
            for (int i = 2; i < totalCase * maxIndRate; i++)
            {
                if (currentList[i] == null)
                    currentList[i] = currentList[i - 1];
                this.controlResult[i].Add(currentList[i]);
            }
            
        }
        public void CalculateControlResultParameters() //calculate the parameters for each control result group.
        {
            this.controlResultParameters = new List<List<double>>();
            foreach (List<GroupShare> newList in this.controlResult)
            {
                double average = 0;
                double std = 0;
                double count = 0;
                foreach (GroupShare newGroup in newList)
                {
                    average += newGroup.AValue;
                    count++;

                }
                average = average / count;
                foreach (GroupShare newGroup in newList)
                {
                    std += Math.Pow((newGroup.AValue - average), 2);
                    count++;

                }
                std = Math.Sqrt(std);
                List<double> newData = new List<double>();
                newData.Add(average);
                newData.Add(std);
                this.controlResultParameters.Add(newData);
            }
        }
        public double GetPValue(GroupShare newGroupShare) //calculate the p value from A value for all results from case.
        { 
            int indSize = newGroupShare.individualList.Count;
            double zScore = (newGroupShare.AValue - this.controlResultParameters[indSize - 1][0]) / this.controlResultParameters[indSize - 1][1];

            if (zScore < -4)
                return 1.0;
            if (zScore > 3.9)
                return 0.0001;
            zScore = (zScore + 4.0) * 100.0;
            return this.pTable[Convert.ToInt32(zScore), 0];

        }
        public bool CheckDuplication (GroupShare groupA, GroupShare groupB) //true for duplicated region.
        {
            if (groupA.individualList.Count != groupB.individualList.Count)
                return false;
            for (int i = 0; i < groupA.individualList.Count; i++)
            {
                if (groupA.individualList[i] != groupB.individualList[i])
                    return false;
            }
            if ((groupA.allSharedFirstBlockID >= groupB.allSharedFirstBlockID && groupA.allSharedFirstBlockID <= groupB.allSharedLastBlockID) || (groupA.allSharedLastBlockID >= groupB.allSharedFirstBlockID && groupA.allSharedLastBlockID <= groupB.allSharedLastBlockID))
            //check the region is overlaped or not
            {
                return true;
            }
            else
                return false;
        }
        public ResultList()
        {
            this.result = new List<GroupShare>();
            
        }
        public void Combine (List<GroupShare> newList)
        {
            foreach (GroupShare newGroupShare in newList)
                this.result.Add(newGroupShare.GetCopy());
            newList = null;
        }
        public void ClearReplication()
        {
            MyResultComparer newResultComparer = new MyResultComparer();
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

        }

        public void calculateAllPvalue(List<GroupShare> caseResultList, int totalCase, double maxIndRate, Commands newCommands)
        { 
            double pValue = 0;
            foreach (GroupShare newGroupShare in caseResultList)
            {
                if (newGroupShare.individualList.Count > totalCase * maxIndRate)
                {
                    newGroupShare.AValue = 1;
                    continue;
                }
                pValue = this.GetPValue(newGroupShare);
                //Console.WriteLine("\t\tAvalue: {0}\tind:{1}\taverage:{2}\tstd:{3}\tpvalue:{4}", newGroupShare.AValue, newGroupShare.individualList.Count, this.controlResultParameters[newGroupShare.individualList.Count - 1][0], this.controlResultParameters[newGroupShare.individualList.Count - 1][1], pValue);
                newGroupShare.AValue = pValue;
            }

        }
        public void PrintResult(BlockDictionary newBlockDictionary, int totalCase, double maxIndRate, Commands newCommands)
        {
            MyResultAValueComparer newComparer = new MyResultAValueComparer();
            result.Sort(newComparer);

            FileStream cFile = new FileStream(newCommands.outputName, FileMode.Append);
            StreamWriter sw1 = new StreamWriter(cFile);
            if (result.Count == 0)
                return;
            
            int count = 0;
            for (int i = 0; i < result.Count; i++)
            {
                
                if (i>200&&(((GroupShare)result[i]).AValue > 0.05 || ((GroupShare)result[i]).individualList.Count > totalCase*maxIndRate))
                {
                    //limit++;
                    continue;
                }
                Block firstBlock = ((Block)(newBlockDictionary.blockList[((GroupShare)result[i]).allSharedFirstBlockID]));
                Block lastBlock = ((Block)(newBlockDictionary.blockList[((GroupShare)result[i]).allSharedLastBlockID]));
                double a = 1.0 / ((GroupShare)result[i]).AValue;
                sw1.Write("{0}\t", ((GroupShare)result[i]).AValue);
                sw1.Write("{0}\t{1}\t{2:00000.0}\t{3}:{4}\t{5}\t", newCommands.chromosome, ((GroupShare)result[i]).individualList.Count, (lastBlock.GetEnd() - firstBlock.GetStart()) / 1000, firstBlock.GetStart(), lastBlock.GetEnd(), lastBlock.GetEndSnpID() - firstBlock.GetStartSnpID() + 1);
                for (int j = 0; j < ((GroupShare)result[i]).individualList.Count; j++)
                {
                    if (j < ((GroupShare)result[i]).individualList.Count - 1)
                        sw1.Write("{0},", ((GroupShare)result[i]).individualList[j]);
                    else
                        sw1.Write("{0}\t", ((GroupShare)result[i]).individualList[j]);
                }
                sw1.Write("{0}:{1}\t", firstBlock.GetStartSnpID(), lastBlock.GetEndSnpID());

                for (int j = 0; j < ((GroupShare)result[i]).individualList.Count; j++)
                {
                    if (j < ((GroupShare)result[i]).individualList.Count - 1)
                        sw1.Write("{0}:{1},", (((Block)(newBlockDictionary.blockList[((GroupShare)result[i]).partiallySharedRegion[j, 0]])).GetStartSnpID()), (((Block)(newBlockDictionary.blockList[((GroupShare)result[i]).partiallySharedRegion[j, 1]])).GetEndSnpID()));
                    else
                        sw1.Write("{0}:{1}\t", (((Block)(newBlockDictionary.blockList[((GroupShare)result[i]).partiallySharedRegion[j, 0]])).GetStartSnpID()), (((Block)(newBlockDictionary.blockList[((GroupShare)result[i]).partiallySharedRegion[j, 1]])).GetEndSnpID()));
                }

                for (int blockID = ((GroupShare)result[i]).allSharedFirstBlockID; blockID <= ((GroupShare)result[i]).allSharedLastBlockID; blockID++)
                {
                    Block currentBlock = ((Block)(newBlockDictionary.blockList[blockID]));
                    sw1.Write("|");
                    for (int snpID = 0; snpID <= currentBlock.GetEndSnpID() - currentBlock.GetStartSnpID(); snpID ++)
                        sw1.Write("{0}",currentBlock.Gethaplotype((int)(((GroupShare)result[i]).haplotypeList[blockID - ((GroupShare)result[i]).allSharedFirstBlockID]),snpID));
                    if (blockID == ((GroupShare)result[i]).allSharedLastBlockID)
                        break;
                    sw1.Write("|");
                    for (int snpID = currentBlock.GetEndSnpID()+1; snpID < ((Block)(newBlockDictionary.blockList[blockID + 1])).GetStartSnpID(); snpID++)
                        sw1.Write("-");
                }

                sw1.Write("\t");
                for (int blockID = ((GroupShare)result[i]).allSharedFirstBlockID; blockID <= ((GroupShare)result[i]).allSharedLastBlockID; blockID++)
                { 
                    sw1.Write("{0},",((Block)(newBlockDictionary.blockList[blockID])).GetFrequency((int)(((GroupShare)result[i]).haplotypeList[blockID - ((GroupShare)result[i]).allSharedFirstBlockID])).ToString("0.0000"));
                }
                sw1.WriteLine();
                if (i > 4)
                    count++;

            }
            Console.Write("{0} regions found in cases, {1} regions found to be with significant p value (0.05).\nFinished writing significant results in outputfile.",this.result.Count,count);
            sw1.Close();
        }

        private class MyResultAValueComparer : IComparer<GroupShare>
        {
            public int Compare(GroupShare x, GroupShare y)
            {
                if (((GroupShare)y).AValue != ((GroupShare)x).AValue)
                    return -((int)((((GroupShare)y).AValue - ((GroupShare)x).AValue) * 100000));
                else
                    return 0;

            }
        }

    }


}
