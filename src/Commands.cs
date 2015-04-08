using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.IO;
using System.Collections;

namespace IBD
{
    class Commands
    {
        public string filename;
        public string pedName;
        public string mapName;
        public string outputName;
        public int permutation;
        public double cut;
        public int chromosome;
        public void HelpInfo()
        {
            Console.Write(
                "Commands:\n"+
                "\t-h\thelp\n"+
                "\t-f\t.ped and .map filename (should be the same)\n" +
                "\t-c\tchromosome 1-23, type 'all' for analyze all chromosomes\n" +
                "\t-o\toutput file name\n");
        }
        public bool ReadCommand(string[] args)
        {
            try
            {
                Regex fileCommand = new Regex(@"^-file");
                Regex chrCommand = new Regex(@"^-chr");
                Regex helpCommand = new Regex(@"^-h");
                Regex outCommand = new Regex(@"^-out");
                Regex numberCommand = new Regex(@"[0-9]+");
                Regex permutationCommand = new Regex(@"-perm");
                Regex cutoffCommand = new Regex(@"-cut");
                Regex commandMatch = new Regex(@"-");
                if (args.Length == 0 || helpCommand.IsMatch(args[0]))
                {
                    this.HelpInfo();
                    return false;
                }
                if (args.Length < 2)
                {
                    Console.WriteLine("\tPlease input your file name.");
                    this.HelpInfo();
                    return false;
                }
                else if (fileCommand.IsMatch(args[0]) || !commandMatch.IsMatch(args[1]))
                {
                    this.filename = args[1];
                    this.pedName = args[1] + ".ped";
                    this.mapName = args[1] + ".map";
                    //Console.WriteLine("ped file is {0}", this.pedName);
                    try
                    {
                        FileStream pFile = new FileStream(this.pedName, FileMode.Open);
                        pFile.Close();
                        FileStream mFile = new FileStream(this.mapName, FileMode.Open);
                        mFile.Close();
                    }
                    catch (IOException e)
                    {
                        Console.WriteLine("\tCan not open file {0} or {1}.", this.pedName, this.mapName);
                        this.HelpInfo();
                        //Console.WriteLine(e.ToString());
                        return false;
                    }
                    Console.WriteLine("\tOpen file {0} and {1} successfully.", this.pedName, this.mapName);
                    if (args.Length < 4)
                    {
                        this.HelpInfo();
                        return false;
                    }
                    if (chrCommand.IsMatch(args[2]) || !commandMatch.IsMatch(args[3]))
                    {
                        if (args[3] == "all")
                        {
                            this.chromosome = 0;
                            Console.WriteLine("\tAll chromosomes will be analyzed.");

                        }
                        else if (numberCommand.IsMatch(args[3]) && Convert.ToInt32(args[3]) >= 0 && Convert.ToInt32(args[3]) < 24)
                        {
                            this.chromosome = Convert.ToInt32(args[3]);
                            Console.WriteLine("\tChromosome {0} will be analyzed.", this.chromosome);
                            try
                            {
                                FileStream gFIle = new FileStream("GeneticDistance\\" + "genetic_map_chr" + this.chromosome + "_b36.txt", FileMode.Open);
                                gFIle.Close();
                            }
                            catch (IOException e)
                            {
                                Console.WriteLine("Can not open file genetic_map_chr{0}_b36.txt.", this.chromosome);
                                Console.WriteLine("Please check it is in the geneticDistance folder.");
                                Console.WriteLine(e.ToString());
                                return false;
                            }
                            finally
                            {
                                Console.WriteLine("\tOpen file genetic_map_chr{0}_b36.txt successfully.", this.chromosome);
                            }
                        }
                        else
                        {
                            Console.WriteLine("\tPlease specify which chromosome to analyse.");
                            this.HelpInfo();
                            return false;
                        }


                        if (permutationCommand.IsMatch(args[4]) || !commandMatch.IsMatch(args[5]))
                        {
                            this.permutation = Convert.ToInt32(args[5]);
                            
                            
                            Console.WriteLine("\tThe number of permutation is set to be {0} times.", this.permutation);
                            
                        }
                        else
                        {
                            Console.WriteLine("\tPlease specify the number of permutation in controls.");
                            return false;
                        }

                        if (permutationCommand.IsMatch(args[6]) || !commandMatch.IsMatch(args[7]))
                        {
                            this.cut = Convert.ToDouble(args[7]);


                            Console.WriteLine("\tThe cutoff for the length of pairwise IBD.is set to be {0}cM.", this.cut);

                        }
                        else
                        {
                            Console.WriteLine("\tPlease specify the cutoff for length of pairwise IBD (cm).");
                            return false;
                        }

                        if (args.Length < 10)
                        {
                            Console.WriteLine("\tThe outfile name will be defult.");
                            return true;
                        }
                        if (outCommand.IsMatch(args[8]) || !commandMatch.IsMatch(args[9]))
                        {
                            this.outputName = args[9] + ".txt";
                            Console.WriteLine("\tThe output file name is {0}", this.outputName);
                            return true;
                        }
                        else
                        {
                            Console.WriteLine("\tThe outfile name will be defult.");
                            return true;
                        }

                    }
                    else
                    {
                        Console.WriteLine("\tPlease specify the chromosome to analyse:");
                        this.HelpInfo();
                    }

                }
                else
                {
                    Console.WriteLine("\tPlease specify the file name.");
                    this.HelpInfo();
                }
                return true;
            }
            catch (Exception e){
                this.HelpInfo();
                return false;
            }
        }
    }
}
