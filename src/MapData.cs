using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;
using System.IO;


namespace IBD
{
    class MapData
    {
        public int[] physicalDistance;
        public double[] geneticDistance;
        public void SetphysicalDistance(int snpID, int position)
        {
            this.physicalDistance[snpID] = position;
        }
        public int GetphysicalDistance(int snpID)
        {
            return this.physicalDistance[snpID];
        }
        public double GetGeneticDistance(int snpID)
        {
            return this.geneticDistance[snpID];
        }
        public void ReadGeneticDistance(GeneticDistance geneticDistance)
        {
            int databaseIndex = 0;
            for (int mapIndex = 0; mapIndex < this.physicalDistance.Length; mapIndex++)
            {
                while (geneticDistance.GetGeneticDistance(0,databaseIndex) != 0 || databaseIndex == 0)
                {
                    if (geneticDistance.GetPhysicalDistance(0, databaseIndex) == this.physicalDistance[mapIndex])
                    {
                        this.geneticDistance[mapIndex] = geneticDistance.GetGeneticDistance(0, databaseIndex);
                        if (mapIndex != 0 && this.geneticDistance[mapIndex - 1] == 0)
                        { 
                            int i = 1;
                            while (mapIndex > i && this.geneticDistance[mapIndex - i] == 0)
                                i++;
                            double geneticRate = (this.geneticDistance[mapIndex] - this.geneticDistance[mapIndex - i]) / (this.physicalDistance[mapIndex] - this.physicalDistance[mapIndex - i]);
                            for (int j = i - 1; j > 0; j--)
                                this.geneticDistance[mapIndex - j] = this.geneticDistance[mapIndex - i] + geneticRate * (this.physicalDistance[mapIndex - j] - this.physicalDistance[mapIndex - i]);
                        }
                          
                        databaseIndex++;
                        break;
                    }
                    else if (geneticDistance.GetPhysicalDistance(0, databaseIndex) <= this.physicalDistance[mapIndex])
                        databaseIndex++;
                    else
                        break;
                }
            }
        }
        public void ReadMapFile(Commands newCommands)
        {
            
            ArrayList collectionArray = new ArrayList();
            string strLine;  //store every line of data
            string[] strLineArray; //store every line of data into array
            char[] splitArray = new char[] { '\t'}; //split by a tab
            FileStream aFile = new FileStream(newCommands.mapName, FileMode.Open);
            StreamReader sr = new StreamReader(aFile);
            strLine = sr.ReadLine();
            while (strLine != null)
            {
                collectionArray.Add(strLine);
                strLine = sr.ReadLine();
            }
            sr.Close();
            physicalDistance = new int[collectionArray.Count];
            for (int snpID = 0; snpID < collectionArray.Count; snpID++)
            {
                String R = (String)collectionArray[snpID];
                strLineArray = R.Split(splitArray);
                this.SetphysicalDistance(snpID, Convert.ToInt32(strLineArray[3]));
            }
            this.geneticDistance = new double[this.physicalDistance.Length];
        }
    }
}
