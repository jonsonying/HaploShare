using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;
using System.IO;

namespace IBD
{
    class Main:Object
    {
        public void Evalation(GeneticDistance newGeneticDistance, MapData newMapData, BlockDictionary newBlockDictionary, Commands newCommands)
        {
                ResultList finalResultList = new ResultList();
                finalResultList.SetpTable();
                CaseControlEvaluation(newGeneticDistance, newMapData, newBlockDictionary, newCommands, finalResultList);
        }
        public void CaseControlEvaluation(GeneticDistance newGeneticDistance, MapData newMapData, BlockDictionary newBlockDictionary, Commands newCommands, ResultList finalResultList)
        {
            SnpDataSet caseData = new SnpDataSet();  //Creat Case data set
            SnpDataSet controlData = new SnpDataSet();//Creat Control data set

            Console.Write("\tReading ped data...");
            Inputdata inputData = new Inputdata();
            inputData.data(caseData, controlData, newCommands);  //start getting the data and put into case/comtrol set 
            ArrayList controlIndividualList = new ArrayList();
            ArrayList caseIndividualList = new ArrayList();
            Console.Write("finished.\n");
            inputData = null;
            GC.Collect();
            Console.Write("Analyzing LD data...");
            double[,] controlLinkageData = new double[controlData.GetSnpCount(), 7];
            //caseData.Linkage(caseLinkageData);
            controlData.SetIndividualList(controlIndividualList);
            controlData = null;
            GC.Collect();
            Console.Write("finished.\n");
            finalResultList.SetControlResult(caseData.individual.Length, 0.2);
            
            Console.WriteLine("Haplotyping for all the control samples...");


            for (int indIndex = 0; indIndex < controlIndividualList.Count; indIndex++)
                ((Individual)(controlIndividualList[indIndex])).SetIndividualID(indIndex);


            int ind = 0;
            foreach (Individual newIndividual in controlIndividualList) //haplotyping for each individual.
            {
                newIndividual.Haplotyping(newBlockDictionary, newMapData, controlLinkageData);
                ind++;
                Console.SetCursorPosition(8, Console.CursorTop);
                Console.Write("{0} / {1}", ind, controlIndividualList.Count);
            }
            Console.Write("\t\tFinished.\n");


            
            for (int n = 0; n < newCommands.permutation; n++)  //permutation
            {
                Console.Write("Iteration: {0} / {1} \n", n+1, newCommands.permutation);
                
                List<GroupShare> resultList = new List<GroupShare>(); //unfinished result from last window
                List<int> selectControlList = new List<int>();


                selectControlList = this.Select(controlIndividualList.Count, caseData.individual.Length, new Random());
                //selectControlList = this.Select(controlIndividualList.Count, 30, new Random());

                for (int windowStart = 0; windowStart <= newBlockDictionary.blockList.Count - 1; windowStart += 30000)
                {
                    List<GroupShare> newResultList = new List<GroupShare>();
                    int windowEnd = Math.Min(windowStart + 30000, newBlockDictionary.blockList.Count - 1);
                    ArrayList wholePairwiseIBDRegionList = new ArrayList();
                    Console.WriteLine("\t\tSearching for pairwise sharing regions");
                    for (int i = 0; i < selectControlList.Count; i++)
                    {
                        for (int j = i + 1; j < selectControlList.Count; j++)
                        {

                            ((Individual)controlIndividualList[selectControlList[i]]).PairwiseComparison(wholePairwiseIBDRegionList, (Individual)controlIndividualList[selectControlList[j]], (Block[])newBlockDictionary.blockList.ToArray(typeof(Block)), newMapData, newCommands.cut, windowStart, windowEnd);
                        }
                        Console.SetCursorPosition(16, Console.CursorTop);
                        Console.Write("{0} / {1}\tTotal regions: {2}", i + 1, selectControlList.Count, wholePairwiseIBDRegionList.Count);
                    }
                    Console.WriteLine("\t\tFinished.");
                    MyIBDComparer newComparer = new MyIBDComparer();
                    wholePairwiseIBDRegionList.Sort(newComparer);
                    GroupShareFinder newIBDList = new GroupShareFinder(wholePairwiseIBDRegionList);
                    Console.WriteLine("\t\tSrart searching for group sharing regions: ");

                    newIBDList.findGroupIBD(controlIndividualList, newBlockDictionary, controlLinkageData, newResultList, resultList, windowStart, windowEnd);
                    finalResultList.AddNewControlList(newResultList,caseData.individual.Length, 0.2);
                    

                }

            }
            Console.WriteLine("Finished permutation in controls");

            finalResultList.CalculateControlResultParameters();
            Console.WriteLine("Finished puermutation result calculation");



            caseData.SetIndividualList(caseIndividualList);
            caseData = null;
            GC.Collect();
            
           

            Console.WriteLine("Haplotyping for all case samples...");

            for (int indIndex = 0; indIndex < caseIndividualList.Count; indIndex++)
                ((Individual)(caseIndividualList[indIndex])).SetIndividualID(indIndex);


            ind = 0;
            foreach (Individual newIndividual in caseIndividualList) //haplotyping for each individual.
            {
                newIndividual.Haplotyping(newBlockDictionary, newMapData, controlLinkageData);
                ind++;
                Console.SetCursorPosition(8, Console.CursorTop);
                Console.Write("{0} / {1}", ind, caseIndividualList.Count);
            }
            Console.Write("\t\tFinished.\n");



            for (int n = 0; n < 1; n++)  //case evaluation
            {
                Console.Write("\tTotally {0} blocks will be analyzed.\n", newBlockDictionary.blockList.Count);

                List<GroupShare> resultList = new List<GroupShare>(); //unfinished result from last window


                for (int windowStart = 0; windowStart <= newBlockDictionary.blockList.Count - 1; windowStart += 30000)
                {
                    List<GroupShare> newResultList = new List<GroupShare>();
                    int windowEnd = Math.Min(windowStart + 30000, newBlockDictionary.blockList.Count - 1);
                    Console.WriteLine("\tAnalyzing block {0} to {1}", windowStart, windowEnd - 1);
                    ArrayList wholePairwiseIBDRegionList = new ArrayList();
                    Console.WriteLine("\t\tSearching for pairwise sharing regions");
                    for (int i = 0; i < caseIndividualList.Count; i++)
                    {
                        for (int j = i + 1; j < caseIndividualList.Count; j++)
                        {
                            ((Individual)caseIndividualList[i]).PairwiseComparison(wholePairwiseIBDRegionList, (Individual)caseIndividualList[j], (Block[])newBlockDictionary.blockList.ToArray(typeof(Block)), newMapData, newCommands.cut, windowStart, windowEnd);
                        }
                        Console.SetCursorPosition(16, Console.CursorTop);
                        Console.Write("{0} / {1}\tTotal regions: {2}", i + 1, caseIndividualList.Count, wholePairwiseIBDRegionList.Count);
                    }
                    Console.WriteLine("\n\t\tFinished searching pairwise sharing regions");

                    MyIBDComparer newComparer = new MyIBDComparer();
                    wholePairwiseIBDRegionList.Sort(newComparer);
                    GroupShareFinder newIBDList = new GroupShareFinder(wholePairwiseIBDRegionList);
                    Console.WriteLine("\t\tSrart searching for group sharing regions:");

                    newIBDList.findGroupIBD(caseIndividualList, newBlockDictionary, controlLinkageData, newResultList, resultList, windowStart, windowEnd);
                    finalResultList.calculateAllPvalue(newResultList, caseIndividualList.Count, 0.2, newCommands);
                    finalResultList.Combine(newResultList);
                    finalResultList.ClearReplication();
                    
                }

            }
            finalResultList.PrintResult(newBlockDictionary, caseIndividualList.Count, 0.2, newCommands);

            
        }
        
        public List<int> Select(int total, int select, Random rd)
        {


            List<double[]> randomList = new List<double[]>();
            for (int i = 0; i < total; i++)
            {
                double[] newRandom = { i, rd.NextDouble() };
                randomList.Add(newRandom);

            }
            RandomComparer newComparer = new RandomComparer();
            randomList.Sort(newComparer);
            //return randomList;
            List<int> returnList = new List<int>();
            for (int i = 0; i < select; i++)
            {
                returnList.Add((int)(randomList[i][0]));
            }
            returnList.Sort();

            return returnList;
        }
        
        private class MyIBDComparer : IComparer
        {
            public int Compare(Object x, Object y)
            {                                                                       //the order of the blocks
                if (((PairwiseIBDRegion)x).GetStart() - ((PairwiseIBDRegion)y).GetStart() != 0)      //start position smaller first
                    return (((PairwiseIBDRegion)x).GetStart() - ((PairwiseIBDRegion)y).GetStart());
                else if (((PairwiseIBDRegion)x).GetEnd() - ((PairwiseIBDRegion)y).GetEnd() != 0)     //end position smaller first
                    return (((PairwiseIBDRegion)x).GetEnd() - ((PairwiseIBDRegion)y).GetEnd());
                else if (((PairwiseIBDRegion)x).GetFirstIndividual() - ((PairwiseIBDRegion)y).GetFirstIndividual() != 0)
                    return (((PairwiseIBDRegion)x).GetFirstIndividual() - ((PairwiseIBDRegion)y).GetFirstIndividual());
                else if (((PairwiseIBDRegion)x).GetSecondIndividual() - ((PairwiseIBDRegion)y).GetSecondIndividual() != 0)
                    return (((PairwiseIBDRegion)x).GetSecondIndividual() - ((PairwiseIBDRegion)y).GetSecondIndividual());

                return 0;                                                                               //all the same
            }
        }
        private class RandomComparer : IComparer<double[]>
        {
            public int Compare(double[] x, double[] y)
            {
                if (((double[])(x))[1] > ((double[])(y))[1])
                    return 1;
                else if (((double[])(x))[1] < ((double[])(y))[1])
                    return -1;
                return 0;
            }
        }

    }
}
