using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;
using System.IO;

namespace IBD
{
    class ShareWindow   //store the individuals that can share the whole window.
    {
        public Window windowStart;
        public Window windowLast;
        public ArrayList group;
        private void Remove(ArrayList windowInSmaePosition, ArrayList removeWindowList)
        {
            foreach (ShareWindow sW1 in windowInSmaePosition)
            {
                foreach (ShareWindow sW2 in windowInSmaePosition)
                {
                    if (sW1 == sW2) { continue; }
                    int i = 0;
                    foreach (int ind1 in sW1.group)
                    {
                        foreach (int ind2 in sW2.group)
                        {
                            if (ind1 == ind2) i++;
                            
                        }
                    }
                    if (i == sW1.group.Count)
                    {
                        //ShareWindowInSmaePosition.Remove(sW1);
                        removeWindowList.Add(sW1);
                        break;
                    }
                }
            }
        }
        public void RemoveOverlap(ArrayList shareWindowList)
        {   
            int i=-1;
            ArrayList ShareWindowInSamePosition = new ArrayList();
            ArrayList removeWindowList = new ArrayList();
            
            foreach(ShareWindow shareWindow1 in shareWindowList)
            {
                if (i < shareWindow1.windowStart.id)
                {
                    if (i >= 0)
                    {
                        this.Remove(ShareWindowInSamePosition,removeWindowList);
                    }
                    ShareWindowInSamePosition.Clear();
                    ShareWindowInSamePosition.Add(shareWindow1);
                    i++;
                }
                else 
                { 
                    ShareWindowInSamePosition.Add(shareWindow1);
                }
            }
            foreach (ShareWindow shareWindow2 in removeWindowList)
            {
                shareWindowList.Remove(shareWindow2);
            }
        }
        public void Combine(ArrayList shareWindowList, ArrayList combinedShareWindowList)   //combine a sharewindow with next window(id+1)
        {
            
            foreach (ShareWindow shareWindowItem in shareWindowList)
            {
                if (shareWindowItem.windowStart.id == this.windowLast.id + 1)
                {
                    //Console.Write("{0}\t{1}\t{2}\r\n", shareWindowItem.windowStart.id,this.group.Count,shareWindowItem.group.Count);
                    //Console.ReadKey();
                    ArrayList newGroup = new ArrayList();
                    foreach (int individualA in this.group)  //find the share individuals between two groups
                    {
                        foreach (int individualB in shareWindowItem.group)
                        {
                            if (individualA == individualB)
                            {
                                newGroup.Add(individualA);
                            }
                        }
                    }
                    if (newGroup.Count >= 5)   //gourp size should be >=3 for combined window
                    {
                        ShareWindow newShareWindow = new ShareWindow();
                        newShareWindow.group = newGroup;
                        newShareWindow.windowStart = this.windowStart;
                        newShareWindow.windowLast = shareWindowItem.windowLast;
                        newShareWindow.Combine(shareWindowList, combinedShareWindowList);
                        if ((this.group.Count - newGroup.Count > 1) && (newShareWindow.windowLast.id-newShareWindow.windowStart.id)>3)
                        {
                            combinedShareWindowList.Add(newShareWindow);
                            Console.Write("{0}\t{1}\t{2}\r\n", newShareWindow.group.Count, newShareWindow.windowStart.id, newShareWindow.windowLast.id);
                        }
                    }
                    else
                    {
                        if (this.windowLast.id - this.windowStart.id > 0)
                        {
                            combinedShareWindowList.Add(this);
                        }
                    }
                }
                else if (shareWindowItem.windowStart.id - this.windowLast.id > 3)
                {
                    break;
                }
            }
        }
    }
}
