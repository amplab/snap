using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;
using System.Threading;
using System.Diagnostics;

namespace CountReadsByContig
{
    class Program
    {
        static ASETools.ThrottledParallelQueue<List<string>> workQueue = new ASETools.ThrottledParallelQueue<List<string>>(1000, 1);
        static Dictionary<bool, Dictionary<string, long>> globalMappedCounts = new Dictionary<bool, Dictionary<string, long>>();  // High quality->contig->count
        static long globalNUnmapped = 0;
        static long globalNSecondary = 0;
        static long globalNReads = 0;

        static void WorkerThread()
        {
            List<string> lines;
            var mappedCounts = new Dictionary<bool, Dictionary<string, long>>();
            long nUnmapped = 0;
            long nSecondary = 0;

            foreach (var hq in ASETools.BothBools)
            {
                mappedCounts.Add(hq, new Dictionary<string, long>());
            }

            while (workQueue.Dequeue(out lines))
            {
                foreach (var nextLine in lines)
                {
                    if (nextLine.StartsWith("@"))
                    {
                        continue;
                    }

                    var samLine = new ASETools.SAMLine(nextLine);

                    if (samLine.isUnmapped())
                    {
                        nUnmapped++;
                        continue;
                    }

                    if (samLine.isSecondaryAlignment())
                    {
                        nSecondary++;
                        continue;
                    }

                    if (!mappedCounts[true].ContainsKey(samLine.rname))
                    {
                        foreach (var hq in ASETools.BothBools)
                        {
                            mappedCounts[hq].Add(samLine.rname, 0);
                        }
                    }

                    mappedCounts[false][samLine.rname]++;
                    if (samLine.mapq >= 10)
                    {
                        mappedCounts[true][samLine.rname]++;
                    }
                } // each line

                lock(globalMappedCounts)
                {
                    if (globalNReads / 1000000 != (globalNReads + lines.Count()) / 1000000)
                    {
                        Console.Write(".");
                    }
                    globalNReads += lines.Count();
                }
            } // each batch of lines from the work queue

            lock (globalMappedCounts)
            {
                globalNUnmapped += nUnmapped;
                globalNSecondary += nSecondary;

                foreach (var hq in ASETools.BothBools) {
                    foreach (var contigEntry in mappedCounts[hq])
                    {
                        if (!globalMappedCounts[hq].ContainsKey(contigEntry.Key))
                        {
                            globalMappedCounts[true].Add(contigEntry.Key, 0);
                            globalMappedCounts[false].Add(contigEntry.Key, 0);
                        }
                        globalMappedCounts[hq][contigEntry.Key] += contigEntry.Value;
                    } // foreach contig
                }// hq
            } // lock
        } // WorkerThread

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            if (args.Count() != 1 || !(args[0].ToLower().EndsWith(".sam") || args[0].ToLower().EndsWith(".bam") || args[0].ToLower().EndsWith(".sam.gz")))
            {
                Console.WriteLine("usage: CountReadsByContig samOrBamFile");
                return;
            }

            foreach (var hq in ASETools.BothBools)
            {
                globalMappedCounts.Add(hq, new Dictionary<string, long>());
            }

            string inputFilename = args[0].ToLower();

            StreamReader inputFile = null;

            if (inputFilename.EndsWith(".sam"))
            {
                inputFile = ASETools.CreateStreamReaderWithRetry(args[0]);
            } else if (inputFilename.EndsWith(".sam.gz"))
            {
                inputFile = ASETools.CreateCompressedStreamReaderWithRetry(args[0]);
            } else if (inputFilename.EndsWith(".bam"))
            {
                var startInfo = new ProcessStartInfo(@"c:\bolosky\bin\samtools.exe", "view -h " + args[0]);
                startInfo.RedirectStandardOutput = true;
                startInfo.UseShellExecute = false;

                Process process =  Process.Start(startInfo);

                inputFile = process.StandardOutput;
            }

            if (inputFile == null)
            {
                Console.WriteLine("Unable to open input file");
                return;
            }

            string nextLine;

            Console.Write("Processing reads (one dot/million): ");

            var workerThreads = new List<Thread>();
            for (int i = 0; i < Environment.ProcessorCount; i++)
            {
                workerThreads.Add(new Thread(WorkerThread));
            }

            workerThreads.ForEach(_ => _.Start());

            int batchSize = 1000;
            var nextBatch = new List<string>();
            while (null != (nextLine = inputFile.ReadLine()))
            {
                nextBatch.Add(nextLine);
                if (nextBatch.Count() >= batchSize)
                {
                    workQueue.Enqueue(nextBatch);
                    nextBatch = new List<string>();
                }
            }

            workQueue.Enqueue(nextBatch);   // This is OK even if it's empty.
            workQueue.TerminateWriter();

            workerThreads.ForEach(_ => _.Join());

            Console.WriteLine();
            Console.WriteLine("Contig\tmapped count (MAPQ >= 10)\tpercentage (MAPQ => 10)\tmapped count (all)\tpercentage (all)");
            Console.WriteLine("unmapped\t" + ASETools.NumberWithCommas(globalNUnmapped) + "\t" + (double)globalNUnmapped * 100 / globalNReads + @"%");
            Console.WriteLine("secondary\t" + ASETools.NumberWithCommas(globalNSecondary) + "\t" + (double)globalNSecondary * 100 / globalNReads + @"%");

            foreach (var mappingEntry in globalMappedCounts[true])
            {
                Console.Write(mappingEntry.Key);

                foreach (var hq in ASETools.BothBools) {
                    Console.Write("\t" + ASETools.NumberWithCommas(globalMappedCounts[hq][mappingEntry.Key]) + "\t" + (double)globalMappedCounts[hq][mappingEntry.Key] * 100 / globalNReads + @"%");
                }
                Console.WriteLine();
            }

            Console.WriteLine();
            Console.WriteLine(ASETools.NumberWithCommas(globalNReads) + " total reads processed in " + ASETools.ElapsedTimeInSeconds(timer));
        }
    }
}
