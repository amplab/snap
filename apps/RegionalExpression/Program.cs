using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.IO.Compression;
using System.Diagnostics;
using System.Threading;
using ExpressionLib;


namespace RegionalExpression
{
    class Program
    {

        class Contig
        {
            public string name = "";
            public long length = -1;
        }

        class Region
        {
            public Region(int regionSize_, Dictionary<string, Dictionary<int, ExpressionTools.MeanAndStdDev>> expression_, Dictionary<string, int> highestOffsetForEachContig_, long nHighQualityMappedNuclearReads_, StreamWriter outputFile_) {
                regionSize = regionSize_;
                expression = expression_;
                highestOffsetForEachContig = highestOffsetForEachContig_;
                nHighQualityMappedNuclearReads = nHighQualityMappedNuclearReads_;
                outputFile = outputFile_;
            }

            Dictionary<string, Dictionary<int, ExpressionTools.MeanAndStdDev>> expression;
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
                if (expression.ContainsKey(contig) && expression[contig].ContainsKey(offset))
                {
                    nBasesExpressedWithBaselineExpression++;
                    double z = (((double)mappedReadCount / (double)nHighQualityMappedNuclearReads) - expression[contig][offset].mean) / expression[contig][offset].stddev;

                    totalZForBasesWithBaselineExpression += z;
                    minZForBasesWithBaselineExpression = Math.Min(z, minZForBasesWithBaselineExpression);
                    maxZForBasesWithBaselineExpression = Math.Max(z, maxZForBasesWithBaselineExpression);
                    totalReadsMappedToBasesWithBaselineExpression += mappedReadCount;

                    totalMuForBasesWithBaselineExpression += ((double)mappedReadCount / (double)nHighQualityMappedNuclearReads) / (double)expression[contig][offset].mean;
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
                if (expression[currentContig].ContainsKey(offset))
                {
                    nBasesWithBaselineButNoLocalExpression++;
                    totalZForBasesWithBaselineExpression -= expression[currentContig][offset].mean / expression[currentContig][offset].stddev; // It's one mean below the mean: ie. 0 expression
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

        static Dictionary<string, Dictionary<int, ExpressionTools.MeanAndStdDev>> expression = null;
        static Dictionary<string, int> highestOffsetForEachContig = new Dictionary<string, int>();

        static int regionSize;

        class OneRun
        {
            public string allcountFilename;
            public string analysis_id;
            public string participantId;
        }

        static void ProcessRuns(List<OneRun> runs)
        {
            while (true) {

                string allcountFilename;
                string analysis_id;
                string participantId;
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
                    analysis_id = runs[0].analysis_id;
                    participantId = runs[0].participantId;

                    runs.RemoveAt(0);
                }

                //
                // Run through the allcount file
                //
                var allcountTimer = new Stopwatch();
                allcountTimer.Start();
                var allcountReader = new StreamReader(new GZipStream(new StreamReader(allcountFilename).BaseStream, CompressionMode.Decompress));

                //
                // The format is two header lines like:
                //      CountReadsCovering v1.1 \\msr-genomics-0\d$\sequence\indices\grch37-24 -a \\msr-genomics-0\e$\tcga\hnsc\tumor\rna\0415c9cc-48e6-46e7-9076-06b6b56bb4be\PCAWG.e7697e88-96df-4c96-ad8d-6c54df95d29b.STAR.v1.bam -
                //      83303682 mapped high quality nuclear reads, 9891899 low quality reads, 0 mitochondrial reads, 100740272 total reads
                // followed by a blank line followed by
                //      NumContigs: 93
                //      ContigName\tLength
                //      <ContigName\tLength> x numContigs
                // Then for each contig:
                //      >ContigName
                // And within each contig one of three line types:
                //      offset\tcount
                // Which is an offset in the contig (in HEX) followed by the count of high quality reads that map there (also in hex) or:
                //      xnumber
                // which is the number of consecutive loci with the same number of HQ reads mapped as the previous locus (x is literal, count is in hex) or:
                //      count
                // which is the new count of reads mapped for the next locus (in hex)
                //
                // And the final line in the file is
                //      **done**
                // 
                //

                var line = allcountReader.ReadLine();

                const string headerBeginning = "CountReadsCovering v";

                if (null == line || line.Count() < headerBeginning.Count() + 1 || line.Substring(0,headerBeginning.Count()) != headerBeginning)
                {
                    Console.WriteLine("Empty or corrupt allcount file " + allcountFilename);
                    continue;
                }

                if (line[headerBeginning.Count()] != '1')
                {
                    Console.WriteLine("Unsupported major version of allcount file " + allcountFilename + ".  Header line: " + line);
                    continue;
                }

                line = allcountReader.ReadLine();
                if (null == line)
                {
                    Console.WriteLine("Corrupt or tuncated allcount file " + allcountFilename);
                    continue;
                }

                var fields = line.Split(' ');
                if (fields.Count() != 16)
                {
                    Console.WriteLine("Corrupt or tuncated allcount file " + allcountFilename + ".  Second line has " + fields.Count() + " fields: " + line);
                    continue;
                }

                long mappedHQNuclearReads = 0;
                try
                {
                    mappedHQNuclearReads = Convert.ToInt64(fields[0]);
                }
                catch (FormatException)
                {
                    Console.WriteLine("Format exception parsing mapped HQ read count for file " + allcountFilename + " from line: " + line);
                    continue;
                }

                line = allcountReader.ReadLine();   // The blank line

                line = allcountReader.ReadLine();
                if (null == line) {
                    Console.WriteLine("Allcount file truncated before contig count: " + allcountFilename);
                    continue;
                }

                const string numContigsLineBeginning = "NumContigs: ";
                if (line.Count() < numContigsLineBeginning.Count() + 1 || line.Substring(0, numContigsLineBeginning.Count()) != numContigsLineBeginning) {
                    Console.WriteLine("Malformed NumContigs line in " + allcountFilename + ": " + line);
                    continue;
                }

                int numContigs = 0;
                try {
                    numContigs = Convert.ToInt32(line.Substring(numContigsLineBeginning.Count()));
                } catch (FormatException) {
                    Console.WriteLine("Couldn't parse NumContigs line in file " + allcountFilename + ": " + line);
                    continue;
                }
                
                if (numContigs < 1) {
                    Console.WriteLine("Invalid numContigs in " + allcountFilename + ": " + line);
                    continue;
                }

                line = allcountReader.ReadLine();   // The header line for the contigs.

                var contigs = new Contig[numContigs];
                bool fileCorrupt = false;

                int whichContig;

                for (whichContig = 0; whichContig < numContigs; whichContig++)
                {
                    contigs[whichContig] = new Contig();

                    line = allcountReader.ReadLine();
                    if (null == line) {
                        Console.WriteLine("File truncated in contig list " + allcountFilename);
                        fileCorrupt = true;
                        break;
                    }

                    fields = line.Split('\t');
                    if (fields.Count() != 2) {
                        Console.WriteLine("Incorrect contig line format in file " + allcountFilename + ": " + line);
                        fileCorrupt = true;
                        break;
                    }

                    contigs[whichContig].name = fields[0].ToLower();
                    try {
                        contigs[whichContig].length = Convert.ToInt64(fields[1]);
                    } catch (FormatException) {
                        Console.WriteLine("Incorrect contig line format in file " + allcountFilename + ": " + line);
                        fileCorrupt = true;
                        break;
                    }
                } // for all expected contigs

                if (fileCorrupt) {
                    continue;
                }

                //
                // Now process the body of the file.
                //

                whichContig = -1;

                int indexOfLastSlash = allcountFilename.LastIndexOf('\\');
                if (-1 == indexOfLastSlash) {
                    Console.WriteLine("Couldn't find a backslash in allcount pathname, which is supposed to be absolute: " + allcountFilename);
                    continue;
                }

                string directory = allcountFilename.Substring(0, indexOfLastSlash + 1);  // Includes trailing backslash
                var outputFilename = directory + analysis_id + ".regional_expression.txt";
                var outputFile = new StreamWriter(outputFilename);

                outputFile.WriteLine("RegionalExpression v3.0\t" + analysis_id + "\t" + allcountFilename + "\t" + regionSize);
                outputFile.WriteLine("NumContigs: " + numContigs);
                Region.printHeader(outputFile);

                Region region = new Region(regionSize, expression, highestOffsetForEachContig, mappedHQNuclearReads, outputFile);

                int currentOffset = -1;
                int currentMappedReadCount = -1;

                bool sawDone = false;
                string strippedContigName = "";

                while (null != (line = allcountReader.ReadLine()))
                {
                    if (sawDone)
                    {
                        Console.WriteLine("File " + allcountFilename + " continues after **done** line: " + line);
                        fileCorrupt = true;
                        break;
                    }

                    if ("**done**" == line)
                    {
                        sawDone = true;
                        continue;
                    }

                    if (line.Count() == 0)
                    {
                        Console.WriteLine("Unexpected blank line in " + allcountFilename);
                        fileCorrupt = true;
                        break;
                    }

                    if (line[0] == '>')
                    {
                        whichContig++;
                        if (whichContig >= numContigs)
                        {
                            Console.WriteLine("Saw too many contigs in " + allcountFilename + ": " + line);
                            fileCorrupt = true;
                            break;
                        }

                        if (line.Substring(1).ToLower() != contigs[whichContig].name)
                        {
                            Console.WriteLine("Unexpected contig in " + allcountFilename + ".  Expected " + contigs[whichContig].name + ", got ", line.Substring(1));
                            fileCorrupt = true;
                            break;
                        }

                        if (line.Count() > 4 && line.Substring(1, 3).ToLower() == "chr")
                        {
                            strippedContigName = line.Substring(4).ToLower();
                        }
                        else
                        {
                            strippedContigName = line.Substring(1).ToLower();
                        }


                        currentOffset = -1;
                        currentMappedReadCount = -1;
                        continue;
                    }

                    if (-1 == whichContig)
                    {
                        Console.WriteLine("Expected contig line after list of contigs, got " + line);
                        fileCorrupt = true;
                        break;
                    }

                    fields = line.Split('\t');
                    if (fields.Count() == 1)
                    {
                        //
                        // Either xRepeatCount or newCount
                        //
                        if (line[0] == 'x')
                        {
                            if (currentMappedReadCount < 0 || currentOffset < 0) {
                                Console.WriteLine("Got unexpected x line " + line);
                                fileCorrupt = true;
                                break;
                            }

                            int repeatCount;
                            try
                            {
                                repeatCount = Convert.ToInt32(line.Substring(1), 16); // 16 means the string is in hex
                                if (repeatCount <= 0)
                                {
                                    Console.WriteLine("Bogus repeat count " + line);
                                    fileCorrupt = true;
                                    break;
                                }
                            }
                            catch (FormatException)
                            {
                                Console.WriteLine("Format exception processing x line " + line);
                                fileCorrupt = true;
                                break;
                            }
                            for (; repeatCount > 0;repeatCount-- )
                            {
                                currentOffset++;
                                region.processBase(strippedContigName, currentOffset, currentMappedReadCount);
                            }
                        }
                        else
                        {
                            //
                            // A new mapped read count line
                            //
                            if (currentOffset <= 0)
                            {
                                Console.WriteLine("Got unexpected mapped read count line " + line);
                                fileCorrupt = true;
                                break;
                            }
                            try
                            {
                                currentMappedReadCount = Convert.ToInt32(line, 16); // 16 means the string is in hex
                                if (currentMappedReadCount <= 0)
                                {
                                    Console.WriteLine("Bogus current count " + line);
                                    fileCorrupt = true;
                                    break;
                                }
                            }
                            catch (FormatException)
                            {
                                Console.WriteLine("Format exception processing count line " + line);
                                fileCorrupt = true;
                                break;
                            }

                            currentOffset++;
                            region.processBase(strippedContigName, currentOffset, currentMappedReadCount);
                        }

                        continue;
                    }

                    if (fields.Count() != 2)
                    {
                        Console.WriteLine("Saw too many fields in line " + line);
                        fileCorrupt = true;
                        continue;
                    }

                    //
                    // An offset + count line
                    //
                    try
                    {
                        currentOffset = Convert.ToInt32(fields[0], 16); // 16 means the string is in hex
                        currentMappedReadCount = Convert.ToInt32(fields[1], 16); // 16 means the string is in hex

                        if (currentOffset <= 0 || currentMappedReadCount <= 0)
                        {
                            Console.WriteLine("Bogus offset + count line " + line);
                            fileCorrupt = true;
                            break;
                        }
                    }
                    catch (FormatException)
                    {
                        Console.WriteLine("Unable to parse offset + count line " + line);
                        fileCorrupt = true;
                        break;
                    }

                    region.processBase(strippedContigName, currentOffset, currentMappedReadCount);
                }

                if (!sawDone)
                {
                    Console.WriteLine("Truncated allcount file " + allcountFilename);
                    fileCorrupt = true;
                }

                region.closeRegion();
                outputFile.WriteLine("**done**");
                outputFile.Close();
                if (fileCorrupt)
                {
                    File.Delete(outputFilename);
                }
                else
                {
                    allcountTimer.Stop();
                    lock (runs)
                    {
                        Console.WriteLine("Processed " + participantId + " in " + (allcountTimer.ElapsedMilliseconds + 500) / 1000 + "s, " + runs.Count() + " remain" +  ((runs.Count() == 1) ? "s" : "") + " in the queue.");
                    }
                }
            }
        }

        static void Main(string[] args)
        {


            if (args.Count() < 3) {
                Console.WriteLine("usage: RegionalExpression expression_disease_file regionSize <one or more participantIDs with the same disease>");
                return;
            }

            Stopwatch timer;
            try
            {
                Console.Write("Loading expression file " + args[0] + "...");
                timer = new Stopwatch();
                timer.Start();

                expression = null;  // Let the garbage collector get rid of the previous one while we're loading the next

                expression = ExpressionTools.LoadExpressionFile(args[0]);   // This can take forever!

                //
                // Now build the list of highest offsets for each contig.
                //
                foreach (var expressionEntry in expression)
                {
                    int highestSoFar = 0;
                    foreach (var regionEntry in expressionEntry.Value)
                    {
                        highestSoFar = Math.Max(highestSoFar, regionEntry.Key);
                    }

                    highestOffsetForEachContig.Add(expressionEntry.Key, highestSoFar);
                }

                timer.Stop();
                Console.WriteLine((timer.ElapsedMilliseconds + 500) / 1000 + "s");
            }
            catch (OutOfMemoryException)
            {
                Console.WriteLine("Out of memory exception loading the expresson file (I'm really not sure why this happens when there's plenty of memory).  Running the binary (& ExpressionLib.dll) locally rather than from a share seems to help, though.");
                return;
            }

            try {
                regionSize = Convert.ToInt32(args[1]);
            } catch(FormatException) {
                Console.WriteLine("Couldn't parse region size from command line");
                return;
            }

            if (regionSize < 1)
            {
                Console.WriteLine("Bogus regionSize from command line");
                return;
            }


            string[] experiments = null;

            bool threw = false;
            do
            {
                try
                {
                    threw = false;
                    experiments = File.ReadAllLines(@"\\gcr\scratch\b99\bolosky\experiments.txt");
                }
                catch (IOException)
                {
                    Console.WriteLine("Error opening experiments.txt.  It's probably open, most likely in Excel.  Sleeping 10s and retrying");
                    threw = true;
                    Thread.Sleep(10000);
                }
            } while (threw);

            var experimentsByParticipantId = new Dictionary<string, string[]>();
            foreach (var experimentLine in experiments)
            {
                var fields = experimentLine.Split('\t');
                experimentsByParticipantId.Add(fields[2], fields);
            }

            var runs = new List<OneRun>();

            for (int i = 2; i < args.Count(); i++)  // for each person we're processing
            {
                string participantId = args[i];
                if (!experimentsByParticipantId.ContainsKey(participantId))
                {
                    Console.WriteLine("Couldn't find participant " + args[i] + ", are you sure it's a correct participant ID? Ignoring");
                    continue;
                }

                string allcountFilename = experimentsByParticipantId[participantId][16];
                if (allcountFilename == "")
                {
                    Console.WriteLine("Participant " + participantId + " doesn't appear to have an allcount file");
                    continue;
                }

                var run = new OneRun();
                run.allcountFilename = allcountFilename;
                run.analysis_id = experimentsByParticipantId[participantId][3];
                run.participantId = participantId;

                runs.Add(run);
            }

            //
            // Process the runs in parallel
            //
            int totalNumberOfExperiments = experiments.Count();
            timer.Reset();
            timer.Start();

            var threads = new List<Thread>();
            for (int i = 0; i < Environment.ProcessorCount; i++) {
                threads.Add(new Thread(() => ProcessRuns(runs)));
            }

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            timer.Stop();
            Console.WriteLine("Processed " + (args.Count() - 2) + " experiments in " + (timer.ElapsedMilliseconds + 500) / 1000 + " seconds");
         }
    }
}
