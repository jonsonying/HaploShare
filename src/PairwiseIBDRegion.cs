using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;

namespace IBD
{
    class PairwiseIBDRegion
    {
        private int firstIndividualID;
        private int secondIndividualID;
        private int start;
        private int end;
        private int firstBlockID;
        private int lastBlockID;
        private ArrayList haplotype;
        
        public int GetStart()
        {
            return this.start;
        }
        public int GetEnd()
        {
            return this.end;
        }
        public ArrayList GetHaplotype(int index)
        {
            return (ArrayList)this.haplotype[index];
        }
        public void SetStart(int newStart)
        {
            this.start = newStart;
        }
        public void SetEnd(int newEnd)
        {
            this.end = newEnd;
        }
        public int GetFirstIndividual()
        {
            return this.firstIndividualID;
        }
        public int GetSecondIndividual()
        {
            return this.secondIndividualID;
        }
        public void SetLastBlockID(int newLastBlockID)
        {
            this.lastBlockID = newLastBlockID;
        }
        public int GetLastBlockID()
        {
            return this.lastBlockID;
        }
        public int GetFirstBlockID()
        {
            return this.firstBlockID;
        }
        public void SetFirstBlockID(int newID)
        {
            this.firstBlockID = newID;
        }
        public void AddHaplotype(ArrayList newHaplotypeList)
        {
            this.haplotype.Add(newHaplotypeList);
        }
        public void ClearHaplotype()
        {
            this.haplotype.Clear();
        }
        public int GetHaplotypeLength()
        {
            return this.haplotype.Count;
        }
        public PairwiseIBDRegion GetCopy()
        {
            PairwiseIBDRegion newPairwiseRegion = new PairwiseIBDRegion(this.firstIndividualID, this.secondIndividualID, this.start, this.end, this.lastBlockID, this.firstBlockID);
            newPairwiseRegion.haplotype = (ArrayList)(this.haplotype).Clone();
            return newPairwiseRegion;
        }
        public PairwiseIBDRegion(int firstID, int secondID, int newStart, int newEnd, int newLastBlock, int newFirstBlock)
        {
            this.firstIndividualID = firstID;
            this.secondIndividualID = secondID;
            this.start = newStart;
            this.end = newEnd;
            this.lastBlockID = newLastBlock;
            this.firstBlockID = newFirstBlock;
            this.haplotype = new ArrayList();
        }
    }
}
