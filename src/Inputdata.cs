using System;
using System.Collections;
using System.Linq;
using System.Text;
using System.IO;

namespace IBD
{
    class Inputdata
    {
        public String[,] finalData;
        public int[,] rareallel;
        public int[,] caseData;
        public int[,] controlData;
        public string[,] alleleMap;
        public int populationSize;
        int caseSize;
        int controlSize;
        int Tem; 

        public void data(SnpDataSet caseDataSet, SnpDataSet controlDataSet, Commands newCommands)
        {   
            string strLine;
            string[] strArray;
            char[] splitArray = new char[] { ' '};
            string[] individualName;
            ArrayList collectionArray;
            try
            {
                populationSize = 0;
                collectionArray = new ArrayList();
                FileStream aFile = new FileStream(newCommands.pedName, FileMode.OpenOrCreate);
                StreamReader sr = new StreamReader(aFile);
                strLine = sr.ReadLine();
                strArray = strLine.Split(splitArray);
                Tem = strArray.Length;
                while (strLine != null)
                {
                    populationSize++;
                    strArray = strLine.Split(splitArray);
                    if (strArray[5] == "2")
                        caseSize++;
                    else
                        controlSize++;
                    collectionArray.Add(strLine);
                    strLine = sr.ReadLine();
                }
                sr.Close();

                SetAlleleMap(newCommands);

                individualName = new string[populationSize];
            
                finalData = new String[populationSize, Tem - 5];
                for (int u = 0; u < collectionArray.Count; u++)
                {
                    String R = (String)collectionArray[u];
                    strArray = R.Split(splitArray);
                    individualName[u] = strArray[1];
                    for (int s = 0; s < Tem - 5; s++)
                    {
         
                        finalData[u, s] = strArray[s + 5];
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

            
            int[,] TransTab = new int[populationSize, Tem - 5];
            rareallel = new int[populationSize,((Tem - 6)/ 2+1)];
            for (int x = 0; x < populationSize; x++)
            {
                rareallel[x, 0] = (int)finalData[x, 0][0]-48;

            }

            for (int m = 1; m < Tem - 6; m = m + 2)
            {

                for (int x = 0; x < populationSize; x++)
                { 
                    if (finalData[x,m][0] == Convert.ToChar(this.alleleMap[(m - 1)/2,0]))
                        TransTab[x, m] = 1;
                    else
                        TransTab[x, m] = 0;
                    if (finalData[x, m + 1][0] == Convert.ToChar(this.alleleMap[(m - 1)/2, 0]))
                        TransTab[x, m + 1] = 1;
                    else
                        TransTab[x, m + 1] = 0;
                    if (finalData[x, m][0] == Convert.ToChar("0"))
                    {
                        TransTab[x, m] = 1;
                        TransTab[x, m + 1] = 0;
                    }
                    rareallel[x, (m + 1) / 2] = TransTab[x, m] + TransTab[x, m + 1];
                }
            }

            controlData = new int[controlSize, (Tem - 6) / 2];
            caseData = new int[caseSize, (Tem - 6) / 2];
            ArrayList caseIndividualName = new ArrayList();
            ArrayList controlIndividualName = new ArrayList();
            ArrayList caseIndividualId = new ArrayList();
            ArrayList controlIndividualId = new ArrayList();
            int casePtr = 0;
            int controlPtr = 0;
            for (int u = 0; u < populationSize; u++)
            {
                    if (rareallel[u, 0] == 1)
                    {
                        controlIndividualName.Add(individualName[u]);
                        controlIndividualId.Add(u);
                        String S = (String)collectionArray[u];
                        for (int m = 1; m <= (Tem - 6) / 2; m++)
                            controlData[controlPtr, m - 1] = rareallel[u, m];
                        controlPtr++;
                    }
                    else if (rareallel[u, 0] == 2)
                    {
                        caseIndividualName.Add(individualName[u]);
                        caseIndividualId.Add(u);
                        String S = (String)collectionArray[u];
                        for (int m = 1; m <= (Tem - 6) / 2; m++)

                            caseData[casePtr, m - 1] = rareallel[u, m];
                        casePtr++;
                    }
            }
            caseDataSet.SetSnp(caseData);
            caseDataSet.SetIndividualName((string[])caseIndividualName.ToArray(typeof(string)));
            caseDataSet.SetIndividual((int[])caseIndividualId.ToArray(typeof(int)));
            caseDataSet.SetSnpCount((Tem - 6) / 2);
            controlDataSet.SetSnp(controlData);
            controlDataSet.SetIndividualName((string[])controlIndividualName.ToArray(typeof(string)));
            controlDataSet.SetIndividual((int[])controlIndividualId.ToArray(typeof(int)));
            controlDataSet.SetSnpCount((Tem - 6) / 2);


        }
        private void SetAlleleMap(Commands newCommands)
        {
            FileStream aFile = new FileStream("chr" + newCommands.chromosome + ".frq", FileMode.OpenOrCreate);
            StreamReader sr = new StreamReader(aFile);
            this.alleleMap = new string[(this.Tem - 6) / 2, 3];
            string newLine = sr.ReadLine();
            newLine = sr.ReadLine();
            string[] strLineArray;
            int count = 0;
            while (newLine != null)
            {
                strLineArray = System.Text.RegularExpressions.Regex.Split(newLine, @"\s+");

                this.alleleMap[count, 0] = strLineArray[3];
                this.alleleMap[count, 1] = strLineArray[5];
                this.alleleMap[count, 2] = strLineArray[4];

                    count++;
                newLine = sr.ReadLine();
            }
        }
    }

}
