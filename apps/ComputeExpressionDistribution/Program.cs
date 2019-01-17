using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;
using ASELib;
using System.Threading;

namespace ComputeExpressionDistribution
{
    class Program
    {
        static ASETools.CommonData commonData;

        class ReadDepth
        {
            public readonly string contig;
            public readonly int location;
            public readonly double fractionOfHQMappedBases;

            public ReadDepth(string contig_, int location_, double fractionOfHQMappedBases_)
            {
                contig = contig_;
                location = location_;
                fractionOfHQMappedBases = fractionOfHQMappedBases_;
            }
        }

        class CaseAndResultQueue
        {
            public readonly ASETools.Case case_;
            public readonly ASETools.ThrottledParallelQueue<ReadDepth> queue = new ASETools.ThrottledParallelQueue<ReadDepth>(1024,1);

            public CaseAndResultQueue(ASETools.Case case__)
            {
                case_ = case__;
            }
        }

        static string chromosome; // This is in non-chr form.
        static int basesInChromosome;

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);
            if (null == commonData)
            {
                return;
            }

            if (commonData.configuration.commandLineArgs.Count() != 2 || !commonData.cases.Any(x => x.Value.disease() == commonData.configuration.commandLineArgs[1]))
            {
                Console.WriteLine("usage: ComputeExpressionDistribution chromosome disease");
                return;
            }

            var disease = commonData.configuration.commandLineArgs[1];
            chromosome = ASETools.chromosomeNameToNonChrForm(commonData.configuration.commandLineArgs[0]);
            basesInChromosome = ASETools.chromosomeSizesByName[chromosome].size;

            if (commonData.cases.Select(x => x.Value).Any(x => x.disease() == disease && x.tumor_rna_allcount_filename == "" || x.tumor_rna_mapped_base_count_filename == ""))
            {
                Console.WriteLine("Some cases don't have a tumor RNA allcount or mapped base count file.");
                return;
            }

            var casesToProcess = commonData.cases.Select(x => x.Value).Where(x => x.disease() == disease).Select(x => new CaseAndResultQueue(x)).ToList();
            //
            // Choose an output directory randomly.
            //
            var random = new Random();
            var whichServer = random.Next() % commonData.configuration.dataDirectories.Count();
            var outputFilename = commonData.configuration.dataDirectories[whichServer] + @"..\" + ASETools.ExpressionDistrbutionByChromosomeDirectory + ASETools.Expression_distribution_filename_base + chromosome + "_" + disease;

            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename, 16*1024*1024); // The 16M parameter sets the output buffer size to avoid the horrid fragmentation that happens with the default 64K

            outputFile.Write("Chromosome\tLocus\tmin\tmax");
            for (int i = 1; i <= 9; i++)
            {
                outputFile.Write("\t" + i + "0th %ile");
            }
            outputFile.WriteLine();

            var processLociThread = new Thread(() => ProcessLoci(outputFile, casesToProcess));
            processLociThread.Start();

            var threading = new ASETools.WorkerThreadHelper<CaseAndResultQueue, int>(casesToProcess, HandleOneCase, null, null);

            // One thread per case.  We read all of the values for a given locus of all of the threads at one time, write them to the file and
            // repeat, so we don't have enormous memory requirements.

            threading.start(casesToProcess.Count());
            processLociThread.Join();

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine("Run time " + ASETools.ElapsedTimeInSeconds(commonData.timer));
        } // Main

        static void HandleOneCase(CaseAndResultQueue caseAndResultQueue, int state)
        {
            var allcountReader = new ASETools.AllcountReader(caseAndResultQueue.case_.tumor_rna_allcount_filename);

            long nMappedHQNuclearReads;
            int numContigs;

            if (!allcountReader.openFile(out nMappedHQNuclearReads, out numContigs))
            {
                throw new Exception("Unable to open allcount file " + caseAndResultQueue.case_.tumor_rna_allcount_filename);
            }

            var mappedBaseCount = ASETools.MappedBaseCount.readFromFile(caseAndResultQueue.case_.tumor_rna_mapped_base_count_filename);

            if (!allcountReader.ReadAllcountFile((contig, locus, readCount) => processBase(caseAndResultQueue, mappedBaseCount.mappedBaseCount, contig, locus, readCount), chromosome))  // Only reads the chromosome we want
            { 
                throw new Exception("allcount reader failed for " + caseAndResultQueue.case_.tumor_rna_allcount_filename);
            }

            caseAndResultQueue.queue.TerminateWriter();
        }

        static void processBase(CaseAndResultQueue caseAndResultQueue, long mappedBaseCount, string contig, int locus, int readCount)
        {
                caseAndResultQueue.queue.Enqueue(new ReadDepth(contig, locus, ((double)readCount / mappedBaseCount)));
        }

        class QueueAndNextReadDepth : IComparable<QueueAndNextReadDepth>
        {
            public ASETools.ThrottledParallelQueue<ReadDepth> queue;
            public ReadDepth nextElement;
            public bool done;

            public QueueAndNextReadDepth(ASETools.ThrottledParallelQueue<ReadDepth> queue_)
            {
                queue = queue_;
                done = !queue.Dequeue(out nextElement);
            }

            public int CompareTo(QueueAndNextReadDepth peer)
            {
                if (done)
                {
                    if (peer.done) return 0;
                    return 1;
                }

                if (peer.done)
                {
                    return -1;
                }

                if (nextElement.contig != peer.nextElement.contig)
                {
                    return nextElement.contig.CompareTo(peer.nextElement.contig);
                }

                return nextElement.location.CompareTo(peer.nextElement.location);
            }

            public static QueueAndNextReadDepth Min(List<QueueAndNextReadDepth> queue)
            {
                var smallestSoFar = queue[0];

                int n = queue.Count();

                for (int i = 1; i < n; i++)
                {
                    if (smallestSoFar.CompareTo(queue[i]) > 0)
                    {
                        smallestSoFar = queue[i];
                    }
                }

                return smallestSoFar;
            }

            public void advance()
            {
                done = !queue.Dequeue(out nextElement);
            }
        } // QueueAndNextReadDepth


        // 
        // The consumer thread.
        //
        static void ProcessLoci(StreamWriter outputFile, List<CaseAndResultQueue> casesToProcess)
        {
            var inputQueues = casesToProcess.Select(_ => new QueueAndNextReadDepth(_.queue)).ToList();

            int nCases = casesToProcess.Count();
            int minCasesForValue = (nCases * 3) / 10;    // Need at least 30% to make a call at all.

            int nBasesPerDot = 1;
            long nBasesProcessed = 0;
            long lastDotAt = 0;
            long lastChromosomeBoundaryAt = 0;
            int lastBaseProcessed = 0;
            bool anyProcessed = false;

            while (true)
            {
                var minElement = QueueAndNextReadDepth.Min(inputQueues);

                if (minElement.done)
                {
                    //
                    // An element that's done always compares greater to any other element, so if the
                    // min is done, there's no more data to process.
                    //
                    Console.WriteLine();
                    return;
                }

                if (!anyProcessed)
                {
                    Console.WriteLine("First base processed at " + ASETools.ElapsedTimeInSeconds(commonData.timer));
                    ASETools.PrintMessageAndNumberBar("Processing chromosome " + chromosome, "bases", basesInChromosome, out nBasesPerDot);

                    anyProcessed = true;
                }

                if (minElement.nextElement.location < lastBaseProcessed)
                {
                    //
                    // Chromosome change.
                    //
                    lastChromosomeBoundaryAt = nBasesProcessed;
                    lastBaseProcessed = 0;
                }

                nBasesProcessed += minElement.nextElement.location - lastBaseProcessed;
                while (lastDotAt + nBasesPerDot < nBasesProcessed)
                {
                    Console.Write(".");
                    lastDotAt += nBasesPerDot;
                }

                lastBaseProcessed = minElement.nextElement.location;

                var queuesAtThisLocus = inputQueues.Where(_ => _.CompareTo(minElement) == 0).ToList();
                int nAtThisLocus = queuesAtThisLocus.Count();
                if (nAtThisLocus >= minCasesForValue)
                {
                    var depthsAtThisLocus = queuesAtThisLocus.Select(_ => _.nextElement.fractionOfHQMappedBases).ToList();
                    depthsAtThisLocus.Sort();

                    var output = minElement.nextElement.contig + "\t" + minElement.nextElement.location + "\t" + depthsAtThisLocus[0] + "\t" + depthsAtThisLocus[nAtThisLocus - 1];
                    for (int i = 1; i <= 9; i++)
                    {
                        output += "\t" + depthsAtThisLocus[nAtThisLocus * i / 10];
                    }
                    outputFile.WriteLine(output);
                }

                queuesAtThisLocus.ForEach(_ => _.advance()); // Move up everyone we used (or ignored, if there weren't enough)
            } // while true
        } // ProcessLoci
    }
}
