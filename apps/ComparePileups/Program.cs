using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;

namespace ComparePileups
{
    class Program
    {
        class BAMAndSAMPair
        {
            public string bamFilename;
            public string nameSortedSamFilename;

            public FileStream nameSortedSamFile;

            List<string> readNamesMappingLocus = new List<string>();
        } // BAMAndSAMPair
        static void Main(string[] args)
        {
            if (args.Count() != 5)
            {
                Console.WriteLine("usage: ComparePileups locus alignment0.bam alignment0.name-sorted.sam alignment1.bam alignment1.name-sorted.sam");
                return;
            } // If we don't have the right number of args

            var locus = args[0];
            var inputFiles = new BAMAndSAMPair[2];
            inputFiles[0].bamFilename = args[1];
            inputFiles[0].nameSortedSamFilename = args[2];
            inputFiles[1].bamFilename = args[3];
            inputFiles[1].nameSortedSamFilename = args[4];

            for (int i = 0; i < 2; i++)
            {
                inputFiles[i].nameSortedSamFile = File.OpenRead(inputFiles[i].nameSortedSamFilename);
                if (inputFiles[i].nameSortedSamFile == null)
                {
                    Console.WriteLine("Unable to open " + inputFiles[i].nameSortedSamFilename);
                    return;
                }
            } // foreach input file





        } // Main
    } // Program
} // Namespace
