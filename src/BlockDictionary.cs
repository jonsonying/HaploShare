using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;
using System.IO;

namespace IBD
{

    class BlockDictionary
    {
        private struct Haplotype
        {
            public string hapString;
            public double frequency;
            public Haplotype(string newHapString, double newFrequency)
            {
                this.hapString = newHapString;
                this.frequency = newFrequency;
            }
        }
        public ArrayList blockList;
        public int[,] alleleMap; // for mapping from atgc to 01; 0 for major and 1 for minor allele. exp: [1,]: A 1 G 0; 1 for A and 0 for G.
        public int GetLength()
        {
            return this.blockList.Count;
        }
        public BlockDictionary(Commands newCommends)
        {
            this.blockList = new ArrayList();
            FileStream dictionary = new FileStream(@"newhap_freq_chr" + newCommends.chromosome + ".txt", FileMode.Open);
            StreamReader sr = new StreamReader(dictionary);
            string dictionaryLine = sr.ReadLine();
            char[] splitArray = new char[] { '\t' };
            char[] singleSplitArray = new char[] {' '};
            Block newBlock = new Block();
            ArrayList haplotypeList = new ArrayList();
            while (dictionaryLine != null)
            {
                if (dictionaryLine == "#")
                {
                    if (haplotypeList.Count != 0)
                    {
                        HaplotypeFrequencyCompare haplotypeFrequencyCompare = new HaplotypeFrequencyCompare();
                        haplotypeList.Sort(haplotypeFrequencyCompare);
                        ArrayList newhaplotypeStringList = new ArrayList();
                        double[] frequency = new double[haplotypeList.Count];
                        int j = 0;
                        foreach (Haplotype newHaplotype in haplotypeList)
                        {
                            //string[] tmp = ((String)newHaplotype.hapString).Split(splitArray);
                            char[] tmp2 = new char[((String)newHaplotype.hapString).Length];
                            for (int i = 0; i < ((String)newHaplotype.hapString).Length; i++)
                                //tmp2[i] = Convert.ToChar(tmp[i]);
                                tmp2[i] = newHaplotype.hapString[i];
                            newhaplotypeStringList.Add(tmp2);
                            frequency[j] = newHaplotype.frequency;
                            j++;
                        }
                        newBlock.SetHaplotype(newhaplotypeStringList);
                        newBlock.SetFrequency(frequency);
                        blockList.Add(newBlock.GetCopy());
                        haplotypeList.Clear();
                    }
                    newBlock = new Block();
                    dictionaryLine = sr.ReadLine();
                    if (dictionaryLine != null)
                    {
                        string[] tmp3 = ((String)dictionaryLine).Split(splitArray);
                        newBlock.SetStart(Convert.ToInt32(tmp3[0]));
                        newBlock.SetEnd(Convert.ToInt32(tmp3[1]));
                    }
                }
                else
                { 
                    string[] tmp;
                    tmp = ((String)dictionaryLine).Split(splitArray);
                    Haplotype newHaplotype = new Haplotype(tmp[0], Convert.ToDouble(tmp[1]));
                    if (Convert.ToDouble(tmp[1]) >= 0.001)
                        haplotypeList.Add(newHaplotype);
                }
                dictionaryLine = sr.ReadLine();
            }
            sr.Close();
        }
        public void Localization(MapData localMap)
        {   
            int mapIndex = 0;
            int blockIndex = 0;
            ((Block)this.blockList[blockIndex]).SetStartSnpID(-1);
            while (mapIndex < localMap.physicalDistance.Length && blockIndex < this.blockList.Count)
            {
                if (localMap.physicalDistance[mapIndex] > ((Block)this.blockList[blockIndex]).GetEnd())
                {
                    blockIndex++;
                }
                else if (localMap.physicalDistance[mapIndex] < ((Block)this.blockList[blockIndex]).GetStart())
                {
                    mapIndex++;
                }
                else
                {
                    if ((((Block)this.blockList[blockIndex]).GetStartSnpID() == 0 && blockIndex != 0) || (((Block)this.blockList[blockIndex]).GetStartSnpID() == -1 && blockIndex == 0))
                        ((Block)this.blockList[blockIndex]).SetStartSnpID(mapIndex);
                    else if (((Block)this.blockList[blockIndex]).GetEndSnpID() < mapIndex)
                        ((Block)this.blockList[blockIndex]).SetEndSnpID(mapIndex);
                    mapIndex++;
                }
            }
        }
        private class HaplotypeFrequencyCompare : IComparer
        {
            public int Compare(Object x, Object y)
            {                                                                       //the order of the haplotypes
                return (int)(((Haplotype)y).frequency * 10000 - ((Haplotype)x).frequency * 10000);
            }
        }
    }
}
