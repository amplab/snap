using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;
using ASELib;

namespace MakeRepetitiveRegionMap
{
    class Program
    {
        static ASETools.Configuration configuration;
        static StreamWriter outputFile;
        static long totalBasesInDuplicateRegions;

        static void Main(string[] args)
        {
            configuration = ASETools.Configuration.loadFromFile(args);

            if (null == configuration)
            {
                Console.WriteLine("Giving up because we were unable to load configuration.");
                return;
            }

            if (configuration.commandLineArgs.Count() != 0)
            {
                Console.WriteLine("usage: MakeRepetitiveRegionMap");
                return;
            }

            var timer = new Stopwatch();
            timer.Start();

            outputFile = ASETools.CreateStreamWriterWithRetry(configuration.redundantChromosomeRegionFilename);
            outputFile.WriteLine("Chromosome\tBegin\tEnd");

            var queue = new List<int>();
            for (int i = 1; i <= ASETools.nHumanAutosomes; i++)
            {
                queue.Add(i);
            }

            var threading = new ASETools.WorkerThreadHelper<int, int>(queue, HandleOneChromosome, null, null, 1);

            threading.run();

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine("Took " + ASETools.ElapsedTimeInSeconds(timer) + " to find " + totalBasesInDuplicateRegions + " bases in duplicate regions");
        } // Main

        static void HandleOneChromosome(int chromosomeNumber, int unused)
        {
            var inputFilename = configuration.chromosomeMapsDirectory + "chr" + chromosomeNumber + ASETools.allLociAlignedExtension;
            var inputFile = ASETools.CreateStreamReaderWithRetry(inputFilename);

            if (null == inputFile)
            {
                Console.WriteLine("Unable to open input file " + inputFilename);
                return;
            }

            string inputLine;

            var duplicateLoci = new HashSet<int>();
            int largestDuplicateLocus = 0;

            while (null != (inputLine = inputFile.ReadLine()))
            {
                if (inputLine.Count() == 0 || inputLine[0] == '@')
                {
                    //
                    // Ignore the header lines.
                    //
                    continue;
                }

                var samLine = new ASETools.SAMLine(inputLine);

                if (samLine.isUnmapped() || !samLine.isSecondaryAlignment())
                {
                    //
                    // We only care about secondary alignments here.  Since the input reads are just the reference genome, there should
                    // be a primary alignment to the generated location (or to a location with the same bases) for every read.
                    //
                    continue;
                }

                //
                // The read name is of the form ChromosomeNumber.generatedLocus.  If we get this far, it's a duplicate locus, so add it
                // to the naughty list.
                //
                if (!samLine.qname.Contains('.'))
                {
                    Console.WriteLine("Can't parse read name " + samLine.qname + ": no dot.");
                    continue;
                }

                int locus = Convert.ToInt32(samLine.qname.Substring(samLine.qname.IndexOf('.') + 1));

                duplicateLoci.Add(locus);
                largestDuplicateLocus = Math.Max(locus, largestDuplicateLocus);
            } // for each line in the SAM file

            lock (outputFile)
            {
                bool inDuplicateRegion = false;
                int startOfDuplicateRegion = 0;
                for (int locus = 1; locus <= largestDuplicateLocus; locus++)
                {
                    if (duplicateLoci.Contains(locus))
                    {
                        if (!inDuplicateRegion)
                        {
                            inDuplicateRegion = true;
                            startOfDuplicateRegion = locus;
                        }
                    } else if (inDuplicateRegion)
                    {
                        outputFile.WriteLine(chromosomeNumber + "\t" + startOfDuplicateRegion + "\t" + locus);
                        totalBasesInDuplicateRegions += locus - startOfDuplicateRegion;
                        inDuplicateRegion = false;
                    }
                }// For every locus in this chromosome

                if (inDuplicateRegion)
                {
                    //
                    // It's the largest duplicate locus, so unless there are no duplicate regions in this chromosome, we
                    // should always come here.
                    //
                    outputFile.WriteLine(chromosomeNumber + "\t" + startOfDuplicateRegion + "\t" + largestDuplicateLocus);
                }
            } // lock outputFile

            inputFile.Close();
        }
    }
}
