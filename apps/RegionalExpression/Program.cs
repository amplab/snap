using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.IO.Compression;
using System.Diagnostics;
using System.Threading;
using ASELib;


namespace RegionalExpression
{
    class Program
    {
        class Region
        {
            public Region(int regionSize_, ASETools.ExpressionFile expression_, Dictionary<string, int> highestOffsetForEachContig_, long nHighQualityMappedNuclearReads_, StreamWriter outputFile_) {
                regionSize = regionSize_;
                expression = expression_;
                highestOffsetForEachContig = highestOffsetForEachContig_;
                nHighQualityMappedNuclearReads = nHighQualityMappedNuclearReads_;
                outputFile = outputFile_;
            }

            ASETools.ExpressionFile expression;
            Dictionary<string, int> highestOffsetForEachContig;
            int regionSize;
            long nHighQualityMappedNuclearReads;
            StreamWriter outputFile;
            string currentContig = "";

            int baseOffset = 0;
            int lastBaseSeen = 0;
            long nBasesExpressed = 0;
            long nBasesExpressedWithBaselineExpression = 0;
            long nBasesExpressedWithoutBaselineExpression = 0;
            long totalReadsMappedToBasesWithoutBaselineExpression = 0;
            long totalReadsMappedToBasesWithBaselineExpression = 0;
            long nBasesWithBaselineButNoLocalExpression = 0;
            double totalZForBasesWithBaselineExpression = 0;
            double totalMuForBasesWithBaselineExpression = 0;   // Instead of standard deviations above/below the mean, just means.  This is 0-based (0 expression => 0 mu)
            double minZForBasesWithBaselineExpression = 1000000000;
            double maxZForBasesWithBaselineExpression = -10000000000;

            void ZeroRegion(int regionBaseOffset)
            {
                baseOffset = regionBaseOffset;
                lastBaseSeen = baseOffset - 1;

                closeRegion();
            }

            public void processBase(string contig, int offset, long mappedReadCount)
            {
                if (contig != currentContig) {
                    //
                    // Zero out the rest of the contig.
                    //
                    if (currentContig != "")
                    {
                        closeRegion();

                        for (int i = baseOffset + regionSize; i <= highestOffsetForEachContig[currentContig]; i+= regionSize) {
                            ZeroRegion(i);
                        }
                    }

                    baseOffset = offset - offset % regionSize;
                    lastBaseSeen = baseOffset - 1;
                    currentContig = contig;

                } else if (baseOffset + regionSize < offset) {
                    //
                    // Finish this region, and zero any with no expression.
                    //
                    closeRegion();

                    int baseOffsetForNextRegionWithExpression = offset - offset % regionSize;
                    for (int i = baseOffset + regionSize; i < baseOffsetForNextRegionWithExpression; i += regionSize) {
                        ZeroRegion(i);
                    }

                    baseOffset = offset - offset % regionSize;
                    lastBaseSeen = baseOffset - 1;
                }
 
 
                for (int i = lastBaseSeen + 1; i < offset; i++)
                {
                    processSkippedBase(i);
                }

                nBasesExpressed++;
                ASETools.MeanAndStdDev meanAndStdDev;
                if (expression.getValue(contig, offset, out meanAndStdDev))
               {
                    nBasesExpressedWithBaselineExpression++;
                    double z = (((double)mappedReadCount / (double)nHighQualityMappedNuclearReads) - meanAndStdDev.mean) / meanAndStdDev.stddev;

                    totalZForBasesWithBaselineExpression += z;
                    minZForBasesWithBaselineExpression = Math.Min(z, minZForBasesWithBaselineExpression);
                    maxZForBasesWithBaselineExpression = Math.Max(z, maxZForBasesWithBaselineExpression);
                    totalReadsMappedToBasesWithBaselineExpression += mappedReadCount;

                    totalMuForBasesWithBaselineExpression += ((double)mappedReadCount / (double)nHighQualityMappedNuclearReads) / meanAndStdDev.mean;
                }
                else
                {
                    nBasesExpressedWithoutBaselineExpression++;
                    totalReadsMappedToBasesWithoutBaselineExpression += mappedReadCount;
                }

                lastBaseSeen = offset;
            }

            public void closeRegion()
            {
                for (int i = lastBaseSeen + 1; i < baseOffset + regionSize; i++)
                {
                    processSkippedBase(i);
                }

                outputFile.WriteLine(currentContig + "\t" + baseOffset + "\t" + nBasesExpressed + "\t" + nBasesExpressedWithBaselineExpression + "\t" + nBasesExpressedWithoutBaselineExpression + "\t" + totalReadsMappedToBasesWithBaselineExpression + "\t" +
                    totalReadsMappedToBasesWithoutBaselineExpression + "\t" + nBasesWithBaselineButNoLocalExpression + "\t" + totalZForBasesWithBaselineExpression + "\t" + minZForBasesWithBaselineExpression + "\t" + maxZForBasesWithBaselineExpression + "\t" +
                    totalZForBasesWithBaselineExpression / (double)regionSize + "\t" + totalMuForBasesWithBaselineExpression / regionSize);

                nBasesExpressed = 0;
                nBasesExpressedWithBaselineExpression = 0;
                nBasesExpressedWithoutBaselineExpression = 0;
                totalReadsMappedToBasesWithoutBaselineExpression = 0;
                totalReadsMappedToBasesWithBaselineExpression = 0;
                totalZForBasesWithBaselineExpression = 0;
                totalMuForBasesWithBaselineExpression = 0;
                nBasesWithBaselineButNoLocalExpression = 0;
                minZForBasesWithBaselineExpression = 1000000000;
                maxZForBasesWithBaselineExpression = -10000000000;
            }

            void processSkippedBase(int offset)
            {
                ASETools.MeanAndStdDev meanAndStdDev;
                if (expression.getValue(currentContig, offset, out meanAndStdDev))
                {
                    nBasesWithBaselineButNoLocalExpression++;
                    totalZForBasesWithBaselineExpression -= meanAndStdDev.mean / meanAndStdDev.stddev; // It's one mean below the mean: ie. 0 expression
                    // No need to update totalMu, since 0 expression adds 0 there.
                }
            }

            static public void printHeader(StreamWriter outputFile)
            {
                outputFile.WriteLine("Contig\tContig Offset\tn Bases Expressed\tn Bases Expressed With Baseline Expression\tn Bases Expressed Without Baseline Expression\tTotal Reads Mapped To Bases With Baseline Expression\t" +
                    "Total Reads Mapped To Bases Without Baseline Expression\tCount of bases with baseline expression but not in this sample\tTotal z For Bases With Baseline Expression\tMin z For Bases With Baseline Expression\tMax z For Bases With BaselineExpression\t" +
                    "Mean z for Bases With Baseline Expression\tMean mu for Bases with Baseline Expression");
            }
        }

        static int regionSize;

        class OneRun
        {
            public string allcountFilename;
            public string tumor_rna_file_id;
            public string case_id;
        }

        static void ProcessRuns(List<OneRun> runs, ASETools.ExpressionFile expression)
        {
            while (true) {

                string allcountFilename;
                string tumor_rna_file_id;
                string case_id;
                lock (runs)
                {
                    if (runs.Count() == 0)
                    {
                        //
                        // No more work, we're done.
                        //
                        return;
                    }

                    allcountFilename = runs[0].allcountFilename;
                    tumor_rna_file_id = runs[0].tumor_rna_file_id;
                    case_id = runs[0].case_id;

                    runs.RemoveAt(0);
                }

                //
                // Run through the allcount file
                //
                var allcountTimer = new Stopwatch();
                allcountTimer.Start();

                var allcountReader = new ASETools.AllcountReader(allcountFilename);
                long mappedHQNuclearReads;
                int numContigs;

                bool fileCorrupt = !allcountReader.openFile(out mappedHQNuclearReads, out numContigs);

                if (fileCorrupt)
                {
                    Console.WriteLine("Ignoring corrupt allcount file " + allcountFilename);
                    continue;
                }

                //
                // Now process the body of the file.
                //


                int indexOfLastSlash = allcountFilename.LastIndexOf('\\');
                if (-1 == indexOfLastSlash) {
                    Console.WriteLine("Couldn't find a backslash in allcount pathname, which is supposed to be absolute: " + allcountFilename);
                    continue;
                }

                string directory = ASETools.GetDirectoryFromPathname(allcountFilename);
                var outputFilename = directory + @"\" + tumor_rna_file_id + ".regional_expression.txt";
                var outputFile = new StreamWriter(outputFilename);

                outputFile.WriteLine("RegionalExpression v3.0\t" + tumor_rna_file_id + "\t" + allcountFilename + "\t" + regionSize);
                outputFile.WriteLine("NumContigs: " + numContigs);
                Region.printHeader(outputFile);

                Region region = new Region(regionSize, expression, expression.higestOffsetForEachContig, mappedHQNuclearReads, outputFile);

                ASETools.AllcountReader.ProcessBase processBase = (x, y, z) => region.processBase(x, y, z);

                fileCorrupt = !allcountReader.ReadAllcountFile(processBase);

                region.closeRegion();

                if (fileCorrupt)
                {
                    outputFile.Close();
                    File.Delete(outputFilename);
                }
                else
                {
                    outputFile.WriteLine("**done**");
                    outputFile.Close();

                    allcountTimer.Stop();
                    lock (runs)
                    {
                        Console.WriteLine("Processed " + case_id + " in " + (allcountTimer.ElapsedMilliseconds + 500) / 1000 + "s, " + runs.Count() + " remain" +  ((runs.Count() == 1) ? "s" : "") + " in the queue.");
                    }
                }
            }
        }

        static void Main(string[] args)
        {
            var configuration = ASETools.ASEConfirguation.loadFromFile(args);

            if (configuration.commandLineArgs.Count() < 3 || configuration.commandLineArgs[2] == "-f" && configuration.commandLineArgs.Count() != 4)
            {
                Console.WriteLine("usage: RegionalExpression expression_disease_file regionSize <one or more caseIDs with the same disease | -f inputFilename> ");
                return;
            }

            try
            {
                regionSize = Convert.ToInt32(configuration.commandLineArgs[1]);
            }
            catch (FormatException)
            {
                Console.WriteLine("Couldn't parse region size from command line");
                return;
            }

            if (regionSize < 1)
            {
                Console.WriteLine("Bogus regionSize from command line");
                return;
            }

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

            ASETools.ExpressionFile expression;
            

            Stopwatch timer;
            try
            {
                Console.Write("Loading expression file " + configuration.commandLineArgs[0] + "...");
                timer = new Stopwatch();
                timer.Start();

                expression = new ASETools.ExpressionFile();

                expression.LoadFromFile(configuration.commandLineArgs[0]);   // This can take forever!

                timer.Stop();
                Console.WriteLine((timer.ElapsedMilliseconds + 500) / 1000 + "s");
            }
            catch (OutOfMemoryException)
            {
                Console.WriteLine("Out of memory exception loading the expression file (I'm really not sure why this happens when there's plenty of memory).  Running in the visual studio debugger seems to help, though.");
                return;
            }

            var runs = new List<OneRun>();

            List<string> case_ids;

            if (configuration.commandLineArgs[2] == "-f")
            {
                case_ids = File.ReadAllLines(configuration.commandLineArgs[3]).ToList();
            }
            else
            {
                case_ids = new List<string>();
                for (int i = 2; i < configuration.commandLineArgs.Count(); i++)
                {
                    case_ids.Add(configuration.commandLineArgs[i]);
                }
            }
           
            foreach (var case_id in case_ids)  // for each person we're processing
            {
                if (!cases.ContainsKey(case_id))
                {
                    Console.WriteLine("Couldn't find case_id " + case_id + ", are you sure it's a correct case ID? Ignoring");
                    continue;
                }

                var case_ = cases[case_id];

                string allcountFilename = case_.tumor_rna_allcount_filename;
                if (allcountFilename == "")
                {
                    Console.WriteLine("Case " + case_id + " doesn't appear to have a tumor RNA allcount file");
                    continue;
                }

                var run = new OneRun();
                run.allcountFilename = allcountFilename;
                run.tumor_rna_file_id = case_.tumor_rna_file_id;
                run.case_id = case_id;

                runs.Add(run);
            }

            //
            // Process the runs in parallel
            //
            timer.Reset();
            timer.Start();

            var threads = new List<Thread>();
            for (int i = 0; i < Environment.ProcessorCount; i++) {
                threads.Add(new Thread(() => ProcessRuns(runs, expression)));
            }

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            timer.Stop();
            Console.WriteLine("Processed " + (configuration.commandLineArgs.Count() - 2) + " cases in " + ASETools.ElapsedTimeInSeconds(timer));
         }
    }
}
