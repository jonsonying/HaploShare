using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace IBD
{
    class MyResultComparer : IComparer<GroupShare>
    {
        public int Compare(GroupShare x, GroupShare y)
        {
            if (((GroupShare)y).individualList.Count != ((GroupShare)x).individualList.Count)
                return (((GroupShare)x).individualList.Count - ((GroupShare)y).individualList.Count);
            for (int i = 0; i < ((GroupShare)x).individualList.Count; i++)
            {
                if ((((int)((GroupShare)y).individualList[i]) - ((int)((GroupShare)x).individualList[i])) != 0)
                    return ((int)(((GroupShare)x).individualList[i]) - (int)(((GroupShare)y).individualList[i]));
            }
            //if (((GroupShare)y).lastBlockID != ((GroupShare)x).lastBlockID)
            //  return ((((GroupShare)x).lastBlockID - ((GroupShare)y).lastBlockID));
            if (((GroupShare)y).allSharedLastBlockID != ((GroupShare)x).allSharedLastBlockID)
                return ((((GroupShare)x).allSharedLastBlockID - ((GroupShare)y).allSharedLastBlockID));
            if (((GroupShare)y).allSharedFirstBlockID != ((GroupShare)x).allSharedFirstBlockID)
                return ((((GroupShare)x).allSharedFirstBlockID - ((GroupShare)y).allSharedFirstBlockID));

            return 0;

        }
    }
}
