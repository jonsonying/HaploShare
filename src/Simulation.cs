using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;
using System.IO;
using System.Threading;


namespace IBD
{
    class Simulation
    {
        private int ancenstor;
        private int[] offsprings;
        private int generationLength;
        private int mutationPosition;
        public int[,] inheritPosition; // [offspingID, Start/End]
        private int chromosomeLength;
        public Simulation(int offspringLength, int generation, int chrLength) //constructor
        {
            
            this.offsprings = new int[offspringLength];
            this.generationLength = generation;
            this.SetMutationPosition(chrLength);
            this.chromosomeLength = chrLength;
        }
        private void SetOffspring(int index, int id)
        {
            this.offsprings[index] = id;
        }
        private void SetMutationPosition(int chrLength)
        { 
            this.mutationPosition = this.GetRandomInt(chrLength - 1)/2 + chrLength/4;
        }
        private double GetRandomDouble() //this max can't be reached
        {
            Random randomseed = new Random();
            return randomseed.NextDouble();
        }
        private int GetRandomInt(int max)//this max can be reached
        {
            Random randomseed = new Random();
            return randomseed.Next(max + 1);
        }
        public int GetOffsping(int index)
        {
            return this.offsprings[index];
        }
        public int GetAncenstor()
        {
            return this.ancenstor;
        }
        public void PerformeSimulation(int populationSize, GeneticDistance geneticDistance, MapData mapData, ArrayList individualList, BlockDictionary blockDictionary, double[,] caseLinkageData) //random choose offspringID, simulat all offspring
        {
            inheritPosition = new int[this.offsprings.Length, 2];
            Random randomSeed = new Random();
            this.ancenstor = randomSeed.Next(populationSize);

            for (int i = 0; i < this.offsprings.Length; i++) //sumulate every offspring id.
            {
                int newOffspring = new int();
                newOffspring = randomSeed.Next(populationSize);
                bool found = false;
                for (int j = 0; j < i; j++)     //simulate offspring id, random choose from population
                {
                    if (this.offsprings[j] == newOffspring || newOffspring == this.ancenstor)  //avoid same id.
                    {
                        i--;
                        found = true;
                        break;
                    }
                }
                if (found)
                    continue;
                else
                {
                    this.SetOffspring(i, newOffspring);
                    this.SingleSimulation(i, geneticDistance, mapData);
                    if (this.inheritPosition[i, 0] == 0 && this.inheritPosition[i, 1] == 0)
                        i--;
                }
            }
            this.SetSimulatedData(individualList, blockDictionary, mapData, caseLinkageData);
//Console.ReadKey();
        }
        private void SingleSimulation(int offspringIndex, GeneticDistance geneticDistance, MapData mapData) // simulate one offspring inherit region
        {
            ArrayList recombinationPosition = new ArrayList();
            Random rd = new Random();
            for (int generation = 0; generation < this.generationLength; generation++)
            {
                int databaseSnpIndex = 0;
                
                double lastCM = 0;
                for (int sampleSnpIndex = 0; sampleSnpIndex < mapData.physicalDistance.Length; sampleSnpIndex++)
                { 
                    while (databaseSnpIndex < 300000)
                    {
                        int samplePosition = mapData.GetphysicalDistance(sampleSnpIndex);
                        int databasePosition = geneticDistance.GetPhysicalDistance(0, databaseSnpIndex);
                        
                        if (mapData.GetphysicalDistance(sampleSnpIndex) > geneticDistance.GetPhysicalDistance(0, databaseSnpIndex))
                        {
                            databaseSnpIndex++;
                        }
                        else if (mapData.GetphysicalDistance(sampleSnpIndex) == geneticDistance.GetPhysicalDistance(0, databaseSnpIndex))
                        {
                            if (sampleSnpIndex == 0)
                            {
                                lastCM = geneticDistance.GetGeneticDistance(0, databaseSnpIndex);
                                break;
                            }
                            double random = rd.NextDouble();
                            double limit = geneticDistance.GetGeneticDistance(0, databaseSnpIndex) - lastCM;
                            limit = limit / 10;
                            limit = Math.Exp(-limit) * limit;
                            if (random <= limit)
                            {
                                random = rd.NextDouble();
                                //if (random*3 <= 2)
                                    recombinationPosition.Add(sampleSnpIndex - 1);
                            }
                            lastCM = geneticDistance.GetGeneticDistance(0, databaseSnpIndex);
                            databaseSnpIndex++;

                            break;
                        }
                        else //this snp is not found in the database, set the prob = 0, check next.
                        {
                            break;
                        }
                    }
                }
            }
            recombinationPosition.Sort(); //all are the c event
            double randomCrossOver = rd.NextDouble();
            int crossoverPosition  = (int)(randomCrossOver*5);
            for (int i = crossoverPosition; i < recombinationPosition.Count - 5; i+=5)
            {
                if ((int)recombinationPosition[i] < this.mutationPosition && (int)recombinationPosition[i + 5] >= this.mutationPosition)
                {
                    this.inheritPosition[offspringIndex, 0] = (int)recombinationPosition[i];
                    
                    this.inheritPosition[offspringIndex, 1] = (int)recombinationPosition[i+5];
                    break;
                }
                
            }
        }
        private void SetSimulatedData(ArrayList individualList, BlockDictionary blockDictionary, MapData mapData, double[,] caseLinkageData) // haplotype the ancestor's genotype data, substitude one haplotype for each simulated offspring.
        {
            int[] ancenstorHaplotype1 = new int[((Individual)individualList[this.ancenstor]).genotype.Length];
            
           int[] ancenstorHaplotype2 = new int[((Individual)individualList[this.ancenstor]).genotype.Length];
          FileStream cFile = new FileStream(@"simulation.txt", FileMode.Create);
          StreamWriter sw1 = new StreamWriter(cFile);
            for (int simulationIndex = 0; simulationIndex < this.offsprings.Length + 1; simulationIndex++) //haplotyping for the ancenstor and the offspings.
            {
                int haplotypingIndex;
                if (simulationIndex == 0)
                    haplotypingIndex = this.ancenstor;
                else
                    haplotypingIndex = this.offsprings[simulationIndex - 1];
                ((Individual)individualList[haplotypingIndex]).haplotype = new ArrayList[blockDictionary.blockList.Count];
                int blockID = -1;
                foreach (Block newBlock in blockDictionary.blockList)
                {
                    blockID++;
                    if (newBlock.GetEndSnpID() + newBlock.GetStartSnpID() == 0) continue;
                    int[] blockGenotype = new int[newBlock.GetEndSnpID() - newBlock.GetStartSnpID() + 1];
                    for (int snpID = newBlock.GetStartSnpID(); snpID <= newBlock.GetEndSnpID(); snpID++)
                    {
                        blockGenotype[snpID - newBlock.GetStartSnpID()] = ((Individual)individualList[haplotypingIndex]).genotype[snpID];
                    }
                    ArrayList newHaplotypePairList = newBlock.GenotypeToHaplotype(blockGenotype, caseLinkageData);
                    ((Individual)individualList[haplotypingIndex]).haplotype[blockID] = newHaplotypePairList;

                }
                //Haplotyping for the ancenstor.
                int[] haplotype1 = new int[((Individual)individualList[haplotypingIndex]).genotype.Length];
                int[] haplotype2 = new int[((Individual)individualList[haplotypingIndex]).genotype.Length];

                blockID = 0;
                Random randomseed = new Random();
                for (int snpIndex = 0; snpIndex < ((Individual)individualList[haplotypingIndex]).genotype.Length; snpIndex++)
                {
                    if (blockID >= blockDictionary.blockList.Count || mapData.physicalDistance[snpIndex] < ((Block)blockDictionary.blockList[blockID]).GetStart()) // this snp is not in the current block.
                    {
                        if (snpIndex == 0) // start of the simulated region.
                        {
                            haplotype1[snpIndex] = ((Individual)individualList[haplotypingIndex]).genotype[snpIndex] / 2;
                            haplotype1[snpIndex] = (((Individual)individualList[haplotypingIndex]).genotype[snpIndex] + 1) / 2;
                        }
                        else
                        {
                            if (((Individual)individualList[haplotypingIndex]).genotype[snpIndex] == 0) //two major
                            {
                                haplotype1[snpIndex] = 0;
                                haplotype2[snpIndex] = 0;
                            }
                            else if (((Individual)individualList[haplotypingIndex]).genotype[snpIndex] == 2) //two minor
                            {
                                haplotype1[snpIndex] = 1;
                                haplotype2[snpIndex] = 1;
                            }
                            else // major + minor combination at that snp.
                            {
                                if (haplotype1[snpIndex - 1] == haplotype2[snpIndex - 1]) //if the last snp is not the major+minor combination
                                {
                                    if (randomseed.NextDouble() <= 0.5)
                                    {
                                        haplotype1[snpIndex] = 0;
                                        haplotype2[snpIndex] = 1;
                                    }
                                    else
                                    {
                                        haplotype1[snpIndex] = 1;
                                        haplotype2[snpIndex] = 0;
                                    }
                                }
                                else
                                {
                                    if (randomseed.NextDouble() <= caseLinkageData[snpIndex - 1, 6])
                                    {
                                        haplotype1[snpIndex] = haplotype1[snpIndex - 1];
                                        haplotype2[snpIndex] = haplotype2[snpIndex - 1];

                                    }
                                    else
                                    {
                                        haplotype1[snpIndex] = haplotype2[snpIndex - 1];
                                        haplotype2[snpIndex] = haplotype1[snpIndex - 1];
                                    }
                                }
                            }
                        }
                    }
                    else //the snp is in or after the current block
                    {
                        if (mapData.physicalDistance[snpIndex] > ((Block)blockDictionary.blockList[blockID]).GetEnd()) //after the current block.
                        {
                            blockID++;
                            snpIndex--;
                        }
                        else //in the current block
                        {
                            // the two haplotype in the block is with the same start allele or the last snp is not major+minor.
                            if (((HaplotypePair)((Individual)individualList[haplotypingIndex]).haplotype[blockID][0]).secondHaplotypeID == -1)
                            {
                                for (snpIndex = ((Block)blockDictionary.blockList[blockID]).GetStartSnpID(); snpIndex <= ((Block)blockDictionary.blockList[blockID]).GetEndSnpID(); snpIndex++)
                                {
                                    haplotype1[snpIndex] = ((Block)blockDictionary.blockList[blockID]).Gethaplotype(((HaplotypePair)((Individual)individualList[haplotypingIndex]).haplotype[blockID][0]).firstHaplotypeID, snpIndex - ((Block)blockDictionary.blockList[blockID]).GetStartSnpID());
                                    haplotype2[snpIndex] = ((Individual)individualList[haplotypingIndex]).genotype[snpIndex] - haplotype1[snpIndex];
                                }
                                snpIndex--;
                            }
                            else
                            {
                                if (snpIndex == 0 || ((Individual)individualList[haplotypingIndex]).genotype[snpIndex] != 1 || ((Individual)individualList[haplotypingIndex]).genotype[snpIndex - 1] != 1)// the two haplotype in the block is with the same start allele or the last snp is not major+minor.
                                {
                                    if (randomseed.NextDouble() <= 0.5)
                                        for (snpIndex = ((Block)blockDictionary.blockList[blockID]).GetStartSnpID(); snpIndex <= ((Block)blockDictionary.blockList[blockID]).GetEndSnpID(); snpIndex++)
                                        {
                                            haplotype1[snpIndex] = ((Block)blockDictionary.blockList[blockID]).Gethaplotype(((HaplotypePair)((Individual)individualList[haplotypingIndex]).haplotype[blockID][0]).firstHaplotypeID, snpIndex - ((Block)blockDictionary.blockList[blockID]).GetStartSnpID());
                                            haplotype2[snpIndex] = ((Individual)individualList[haplotypingIndex]).genotype[snpIndex] - haplotype1[snpIndex];
                                        }
                                    else
                                        for (snpIndex = ((Block)blockDictionary.blockList[blockID]).GetStartSnpID(); snpIndex <= ((Block)blockDictionary.blockList[blockID]).GetEndSnpID(); snpIndex++)
                                        {

                                            haplotype2[snpIndex] = ((Block)blockDictionary.blockList[blockID]).Gethaplotype(((HaplotypePair)((Individual)individualList[haplotypingIndex]).haplotype[blockID][0]).firstHaplotypeID, snpIndex - ((Block)blockDictionary.blockList[blockID]).GetStartSnpID());
                                            haplotype1[snpIndex] = ((Individual)individualList[haplotypingIndex]).genotype[snpIndex] - haplotype2[snpIndex];
                                        }
                                    snpIndex--;
                                }
                                else // 1 + 1 major+minor before major + minor in block.
                                {
                                    if (randomseed.NextDouble() <= caseLinkageData[snpIndex - 1, 6]) //AaBb -> AB+ab
                                    {
                                        if (haplotype1[snpIndex - 1] == ((Block)blockDictionary.blockList[blockID]).Gethaplotype(((HaplotypePair)((Individual)individualList[haplotypingIndex]).haplotype[blockID][0]).firstHaplotypeID, 0))
                                        {
                                            for (snpIndex = ((Block)blockDictionary.blockList[blockID]).GetStartSnpID(); snpIndex <= ((Block)blockDictionary.blockList[blockID]).GetEndSnpID(); snpIndex++)
                                            {
                                                haplotype1[snpIndex] = ((Block)blockDictionary.blockList[blockID]).Gethaplotype(((HaplotypePair)((Individual)individualList[haplotypingIndex]).haplotype[blockID][0]).firstHaplotypeID, snpIndex - ((Block)blockDictionary.blockList[blockID]).GetStartSnpID());
                                                haplotype2[snpIndex] = ((Individual)individualList[haplotypingIndex]).genotype[snpIndex] - haplotype1[snpIndex];
                                            }
                                            snpIndex--;
                                        }
                                        else
                                        {
                                            for (snpIndex = ((Block)blockDictionary.blockList[blockID]).GetStartSnpID(); snpIndex <= ((Block)blockDictionary.blockList[blockID]).GetEndSnpID(); snpIndex++)
                                            {

                                                haplotype2[snpIndex] = ((Block)blockDictionary.blockList[blockID]).Gethaplotype(((HaplotypePair)((Individual)individualList[haplotypingIndex]).haplotype[blockID][0]).firstHaplotypeID, snpIndex - ((Block)blockDictionary.blockList[blockID]).GetStartSnpID());
                                                haplotype1[snpIndex] = ((Individual)individualList[haplotypingIndex]).genotype[snpIndex] - haplotype2[snpIndex];
                                            }
                                            snpIndex--;
                                        }
                                    }
                                    else  //AaBb -> Ab+aB
                                    {
                                        if (haplotype1[snpIndex - 1] == ((Block)blockDictionary.blockList[blockID]).Gethaplotype(((HaplotypePair)((Individual)individualList[haplotypingIndex]).haplotype[blockID][0]).firstHaplotypeID, 0))
                                        {
                                            for (snpIndex = ((Block)blockDictionary.blockList[blockID]).GetStartSnpID(); snpIndex <= ((Block)blockDictionary.blockList[blockID]).GetEndSnpID(); snpIndex++)
                                            {
                                                haplotype2[snpIndex] = ((Block)blockDictionary.blockList[blockID]).Gethaplotype(((HaplotypePair)((Individual)individualList[haplotypingIndex]).haplotype[blockID][0]).firstHaplotypeID, snpIndex - ((Block)blockDictionary.blockList[blockID]).GetStartSnpID());
                                                haplotype1[snpIndex] = ((Individual)individualList[haplotypingIndex]).genotype[snpIndex] - haplotype2[snpIndex];
                                            }
                                            snpIndex--;
                                        }
                                        else
                                        {
                                            for (snpIndex = ((Block)blockDictionary.blockList[blockID]).GetStartSnpID(); snpIndex <= ((Block)blockDictionary.blockList[blockID]).GetEndSnpID(); snpIndex++)
                                            {

                                                haplotype1[snpIndex] = ((Block)blockDictionary.blockList[blockID]).Gethaplotype(((HaplotypePair)((Individual)individualList[haplotypingIndex]).haplotype[blockID][0]).firstHaplotypeID, snpIndex - ((Block)blockDictionary.blockList[blockID]).GetStartSnpID());
                                                haplotype2[snpIndex] = ((Individual)individualList[haplotypingIndex]).genotype[snpIndex] - haplotype1[snpIndex];
                                            }
                                            snpIndex--;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if (simulationIndex == 0)
                {
                    ancenstorHaplotype1 = (int[])haplotype1.Clone();
                    ancenstorHaplotype2 = (int[])haplotype2.Clone();
                }
                else
                {
                    for (int m = this.inheritPosition[simulationIndex - 1, 0]; m <= this.inheritPosition[simulationIndex - 1, 1]; m++)
                    {
                        haplotype2[m] = ancenstorHaplotype1[m];
                        ((Individual)individualList[haplotypingIndex]).genotype[m] = haplotype1[m] + haplotype2[m];
                    }
                }

                if (simulationIndex == 0)
                    sw1.Write("{0}\t{1} to {2}\t", haplotypingIndex, 0, 3569);
                
                else
                    sw1.Write("{0}\t{1} to {2}\t", haplotypingIndex, this.inheritPosition[simulationIndex - 1, 0], this.inheritPosition[simulationIndex - 1, 1]);
                
                for (int i = 0; i < ((Individual)individualList[haplotypingIndex]).genotype.Length; i++)
                {
                    if (simulationIndex == 0 || i >= this.inheritPosition[simulationIndex - 1, 0] && i <= this.inheritPosition[simulationIndex - 1, 1])
                        sw1.Write("{0}", ((Individual)individualList[haplotypingIndex]).genotype[i]);
                    else
                        sw1.Write("-");
                }
                sw1.WriteLine();
                sw1.WriteLine();
          
                
            }
            sw1.Close();
        }
    }
}
