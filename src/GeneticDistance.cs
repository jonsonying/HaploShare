using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;
using System.IO;

namespace IBD
{
    class GeneticDistance
    {
        private double[,,] geneticDistance; //[chr, snpIndex, 0.physicalPosition, 1.recombinationRate, 2.geneticDistance, 3.recombination probability]
        public double GetGeneticDistance(int chr, int snpIndex)
        {
            return this.geneticDistance[chr, snpIndex, 2];
        }
        public double GetRecombinationRate(int chr, int snpIndex)
        {
            return this.geneticDistance[chr, snpIndex, 1];
        }
        public double GetRecombinationProbability(int chr, int snpIndex)
        {
            return this.geneticDistance[chr, snpIndex, 3];
        }
        public int GetPhysicalDistance(int chr, int snpIndex)
        {
            return (int)this.geneticDistance[chr, snpIndex,0];
        }
        public void SetGeneticDistance(int chr, int snpIndex, double distance)
        {
            this.geneticDistance[chr, snpIndex, 2] = distance;
        }
        public void SetRecombinationRate(int chr, int snpIndex, double rate)
        {
            this.geneticDistance[chr, snpIndex, 1] = rate;
        }
        public void SetPhysicalDistance(int chr, int snpIndex, double position)
        {
            this.geneticDistance[chr, snpIndex, 0] = position;
        }
        public void SetRecombinationProbability(int chr, int snpIndex, double probability)
        {
            this.geneticDistance[chr, snpIndex, 3] = probability;
        }
        public void InitialGeniticDistance(int chr, int length, Commands newCommands)
        { 
            this.geneticDistance = new double[1, length, 4]; //should be 23 chromosome
            string strLine;  //store every line of data
            string[] strLineArray; //store every line of data into array
            char[] splitArray = new char[] { ' ' }; //split by a space
            ArrayList collectionArray;
            try
            {
                collectionArray = new ArrayList();
                FileStream aFile = new FileStream("GeneticDistance\\"+"genetic_map_chr"+newCommands.chromosome+"_b36.txt", FileMode.Open);
                StreamReader sr = new StreamReader(aFile);
                strLine = sr.ReadLine();
                while (strLine != null)
                {
                    collectionArray.Add(strLine);
                    strLine = sr.ReadLine();
                }
                sr.Close();
                collectionArray.RemoveAt(0);
                for (int snpID = 0; snpID < collectionArray.Count; snpID++)
                {
                    String R = (String)collectionArray[snpID];
                    strLineArray = R.Split(splitArray);
                    this.SetPhysicalDistance(chr, snpID, Convert.ToDouble(strLineArray[0]));
                    this.SetRecombinationRate(chr, snpID, Convert.ToDouble(strLineArray[1]));
                    this.SetGeneticDistance(chr, snpID, Convert.ToDouble(strLineArray[2]));
                    if (snpID == 0)
                    {
                        //this.SetRecombinationProbability(chr, snpID, 0);
                        continue;
                    }
                    else
                    {
                        double probability = (4 + 1) * 2*(this.GetGeneticDistance(chr, snpID) - this.GetGeneticDistance(chr, snpID - 1)) / (double)100;
                        probability = Math.Exp(-probability) * probability;
                        //Console.WriteLine("{0}  {1}", this.GetGeneticDistance(chr, snpID + 1), this.GetGeneticDistance(chr, snpID));
                        this.SetRecombinationProbability(chr, snpID - 1, probability);
                    }
                }
            }
            catch (IOException e)
            {
                Console.WriteLine("An IO exception has been thrown!");
                Console.WriteLine(e.ToString());
                Console.ReadLine();
                return;
            }
        }

    }
}
