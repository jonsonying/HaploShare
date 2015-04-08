using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;


namespace IBD
{
    class Window
    {
        public int start;
        public int end;
        public int id;

        public void inWindowExtend (ArrayList group, int snpPosition, SnpDataSet dataSet, ArrayList shareWindowListRef)
        {
            if (snpPosition == this.end)                //when it reached the end of the window
            {
                ShareWindow newShareWindow = new ShareWindow();             //construct a shareWindow object
                newShareWindow.windowStart = this;
                newShareWindow.windowLast = this;
                newShareWindow.group = group;
                group.Sort();
                shareWindowListRef.Add(newShareWindow);  // Add newShareWindow to shareWindow List;
                return;
                
            }

            ArrayList group0 = new ArrayList();         //store the individuals that the data of this snp is 0
            ArrayList group1 = new ArrayList();         //store the individuals that the data of this snp is 1
            ArrayList group2 = new ArrayList();         //store the individuals that the data of this snp is 2
            
            foreach (int individual in group)
            { 
               if (dataSet.GetSnp(individual,snpPosition) == 0)  // the data is 0
               {
                   group0.Add(individual);
               }
               else if (dataSet.GetSnp(individual, snpPosition) == 1) //the data is 1
               {
                   group1.Add(individual);
               }
               else                                       //the data is 2                    
               {
                   group2.Add(individual);
               }
            }

            if (group2.Count + group1.Count >= 3)           // construct a new group if group2+group1>=3
            {
                ArrayList group12 = new ArrayList();        
                foreach (int ind in group2)                 //add group2 and group1 into group21
                {
                    group12.Add(ind);
                }
                foreach (int ind in group1)                 
                {
                    group12.Add(ind);
                }
                //group12.ToArray();
                inWindowExtend(group12, snpPosition+1,dataSet, shareWindowListRef);     //use group21 to do the next Extend
            }

            if (group0.Count + group1.Count >= 3)            // construct a new group if group0+group1>=3, similar to last.
            {
                ArrayList group01 = new ArrayList();
                foreach (int ind in group0)
                {
                    group01.Add(ind);
                }
                foreach (int ind in group1)
                {
                    group01.Add(ind);
                }
                //group01.ToArray();
                inWindowExtend(group01, snpPosition+1, dataSet, shareWindowListRef);  //do the next Extend
            }
        }
    }
}
