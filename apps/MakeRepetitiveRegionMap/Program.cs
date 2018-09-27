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
        static long totalBasesInDuplicateRegions = 0;
        static long nBasesExcludedOnlyBecauseOfTranscriptome = 0;
        static long nBasesInTranscriptome = 0;
        static long nBasesExcludedFromTranscriptome = 0;

        static Dictionary<string, HashSet<int>> duplicateRegionsFromTranscriptome = new Dictionary<string,HashSet<int>>();

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

            for (int i = 1; i <= ASETools.nHumanAutosomes; i++)
            {
                duplicateRegionsFromTranscriptome.Add("chr" + i, new HashSet<int>());
            }

            var queue = new List<int>();
            for (int i = 1; i <= ASETools.nHumanAutosomes; i++)
            {
                queue.Add(i);
            }

            Console.WriteLine("Processing transcriptome (one dot/chromosome):");
            ASETools.PrintNumberBar(ASETools.nHumanAutosomes);
            var threading = new ASETools.WorkerThreadHelper<int, int>(queue, (a, b) => HandleOneChromosome(true, a, b), null, null, 1);
            threading.run();
            Console.WriteLine(ASETools.ElapsedTimeInSeconds(timer));
            Console.WriteLine(nBasesInTranscriptome + " bases in transcriptome, of which " + nBasesExcludedFromTranscriptome + " (" + ((int)(((double)nBasesExcludedFromTranscriptome / nBasesInTranscriptome) * 100 + 0.5)) + "%) are duplicates.");

            Console.WriteLine("Processing autosomes (one dot/chromosome):");
            ASETools.PrintNumberBar(ASETools.nHumanAutosomes);

            outputFile = ASETools.CreateStreamWriterWithRetry(configuration.redundantChromosomeRegionFilename);
            outputFile.WriteLine("Chromosome\tBegin\tEnd");

            threading = new ASETools.WorkerThreadHelper<int, int>(queue, (a,b) => HandleOneChromosome(false, a, b), null, null, 1);

            threading.run();

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine("Took " + ASETools.ElapsedTimeInSeconds(timer) + " to find " + totalBasesInDuplicateRegions + " bases in duplicate regions, of which " + nBasesExcludedOnlyBecauseOfTranscriptome + " are there only because of the transcriptome.");
        } // Main


        class TranscriptomeContig
        {
            public readonly int start;
            public readonly int end;
            public int length() { return end - start; }
            public readonly string chromosome;

            TranscriptomeContig(string chromosome_, int start_, int end_)
            {
                chromosome = chromosome_;
                start = start_;
                end = end_;

                if (end < start)
                {
                    throw new Exception("TranscriptomeContig: end(" + end + ") < start(" + start + ")");
                }
            }

            public static List<TranscriptomeContig> parse (string contigName)
            {
                var contigTokens = contigName.Split('_');
                if (contigTokens.Count() < 4 || contigTokens[0] != "transcript" || contigTokens[1].Substring(0,3) != "chr" || contigTokens.Count() % 2 != 0)
                {
                    throw new Exception("Malformed transcriptome contig name: " + contigName);
                }

                var contigs = new List<TranscriptomeContig>();
                for (int i = 2; i < contigTokens.Count(); i += 2)
                {
                    contigs.Add(new TranscriptomeContig(contigTokens[1], Convert.ToInt32(contigTokens[i]), Convert.ToInt32(contigTokens[i + 1])));
                }

                return contigs;
            }
        }

        static void HandleOneChromosome(bool handlingTranscriptome, int chromosomeNumber, int unused)
        {
            string chromosome = ((!handlingTranscriptome) ? "chr": "transcriptome_") + chromosomeNumber;
            var inputFilename = configuration.chromosomeMapsDirectory + chromosome + ASETools.allLociAlignedExtension;
            var inputFile = ASETools.CreateStreamReaderWithRetry(inputFilename);
            long nDuplicateLocations = 0;
            long totalLocations = 0;

            if (null == inputFile)
            {
                Console.WriteLine("Unable to open input file " + inputFilename);
                throw new Exception("Unable to open aligned reads from " + inputFilename);
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
                    totalLocations++;
                    continue;
                }

                if (handlingTranscriptome)
                {
                    //
                    // The same region of the genome is represented more than once in our transcriptome map.  So, even secondary mappings aren't
                    // necessarily indicative of a duplication if they're really just the same place in the reference.  Check that here.
                    //

                    var qnameTokens = samLine.qname.Split('.');
                    if (qnameTokens.Count() < 4 || qnameTokens[0] != "trans" || qnameTokens[1].Substring(0,3) != "chr")
                    {
                        throw new Exception("qname in transcriptome read has too few fields or otherwise is malformed: " + samLine.qname);
                    }

                    var initialChromosome = qnameTokens[1];
                    var initialLocus = Convert.ToInt32(qnameTokens[2]);

                    var transcriptomeContigs = TranscriptomeContig.parse(samLine.rname);
                    var mappedChromosome = transcriptomeContigs[0].chromosome;

                    int whichTransctiptomeContig = 0;
                    int mappingOffsetUnaccountedFor = samLine.pos - 1; // -1 because pos is 1-based.
                    while (mappingOffsetUnaccountedFor >= transcriptomeContigs[whichTransctiptomeContig].length())
                    {
                        mappingOffsetUnaccountedFor -= transcriptomeContigs[whichTransctiptomeContig].length();
                        whichTransctiptomeContig++;
                    }

                    int mappedLocus = transcriptomeContigs[whichTransctiptomeContig].start + mappingOffsetUnaccountedFor;

                    if (mappedChromosome == initialChromosome && mappedLocus == initialLocus)
                    {
                        //
                        // Mapped to a different copy of the same place.
                        //
                        continue;
                    }

                    if (!duplicateRegionsFromTranscriptome[initialChromosome].Contains(initialLocus))
                    {
                        duplicateRegionsFromTranscriptome[initialChromosome].Add(initialLocus);
                        nDuplicateLocations++;
                    }
                }
                else
                {
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
                }
            } // for each line in the SAM file

            if (handlingTranscriptome)
            {
                lock (configuration)
                {
                    nBasesInTranscriptome += totalLocations;
                    nBasesExcludedFromTranscriptome += nDuplicateLocations;
                }
            }
            else lock (outputFile)
            {
                bool inDuplicateRegion = false;
                int startOfDuplicateRegion = 0;
                for (int locus = 1; locus <= largestDuplicateLocus; locus++)
                {
                    bool referenceDuplicate;
                    if ((referenceDuplicate = duplicateLoci.Contains(locus)) || duplicateRegionsFromTranscriptome[chromosome].Contains(locus))
                    {
                        if (!referenceDuplicate)
                        {
                            nBasesExcludedOnlyBecauseOfTranscriptome++;
                        }

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
                    totalBasesInDuplicateRegions += largestDuplicateLocus - startOfDuplicateRegion;
                }
            } // lock outputFile

            inputFile.Close();
        }
    }
}
