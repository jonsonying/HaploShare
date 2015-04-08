using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Text.RegularExpressions;
using System.Collections;
using System.Diagnostics;

namespace IBD
{

    class Program
    {

        static void Main(string[] args)
        {
            Commands newCommends = new Commands();
            if (!newCommends.ReadCommand(args))
                return;
            Stopwatch stw = new Stopwatch();
            stw.Start();

            if (newCommends.chromosome == 0)
            {
                Regex r = new Regex(@"[0-9]+");
                FileStream cFile = new FileStream(newCommends.outputName, FileMode.Create);
                cFile.Close();
                for (int i = 1; i < 23; i++)
                {
                    //Console.WriteLine(i);
                    newCommends.chromosome = i;

                    newCommends.pedName = r.Replace(newCommends.pedName, i.ToString());
                    newCommends.mapName = r.Replace(newCommends.mapName, i.ToString());
                    //newCommends.
                    Console.WriteLine("Start reading file data:");
                    Console.Write("\tReading genetic distance data...");
                    GeneticDistance geneticdistance = new GeneticDistance(); //store all the snp genetic distance and recombination data.
                    geneticdistance.InitialGeniticDistance(0, 300000, newCommends);
                    Console.Write("finished.\n");
                    Console.Write("\tReading map data...");

                    MapData mapdata = new MapData();                       //store all the snp physical position in the map file.
                    mapdata.ReadMapFile(newCommends);
                    mapdata.ReadGeneticDistance(geneticdistance);
                    Console.Write("finished.\n");
                    Console.Write("\tReading haplotype frequency data...");

                    BlockDictionary newDictionary = new BlockDictionary(newCommends); //store all the block haplotype information in the Dictionary, the object type is Block.
                    newDictionary.Localization(mapdata);
                    Console.Write("finished.\n");
                    Main newEvaluation = new Main();
                    newEvaluation.Evalation(geneticdistance, mapdata, newDictionary, newCommends);

                }

                return;
            }
            else
            {
                FileStream cFile = new FileStream(newCommends.outputName, FileMode.Create);
                cFile.Close();
                Console.WriteLine("Start reading file data:");
                Console.Write("\tReading genetic distance data...");
                GeneticDistance geneticdistance = new GeneticDistance(); //store all the snp genetic distance and recombination data.
                geneticdistance.InitialGeniticDistance(0, 300000, newCommends);
                Console.Write("finished.\n");
                Console.Write("\tReading map data...");

                MapData mapdata = new MapData();                       //store all the snp physical position in the map file.
                mapdata.ReadMapFile(newCommends);
                mapdata.ReadGeneticDistance(geneticdistance);
                Console.Write("finished.\n");
                Console.Write("\tReading haplotype frequency data...");

                BlockDictionary newDictionary = new BlockDictionary(newCommends); //store all the block haplotype information in the Dictionary, the object type is Block.
                newDictionary.Localization(mapdata);
                Console.Write("finished.\n");
                Main newEvaluation = new Main();
                newEvaluation.Evalation(geneticdistance, mapdata, newDictionary, newCommends);
                return;
            }
            stw.Stop();

        }


    }
}
