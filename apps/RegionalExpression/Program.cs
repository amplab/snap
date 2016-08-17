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
            public Region(int regionSize_, Dictionary<string, Dictionary<int, ExpressionTools.MeanAndStdDev>> expression_, long nHighQualityMappedNuclearReads_, StreamWriter outputFile_) {
                regionSize = regionSize_;
                expression = expression_;
                nHighQualityMappedNuclearReads = nHighQualityMappedNuclearReads_;
                outputFile = outputFile_;
            }

            Dictionary<string, Dictionary<int, ExpressionTools.MeanAndStdDev>> expression;
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
            double minZForBasesWithBaselineExpression = 1000000000;
            double maxZForBasesWithBaselineExpression = -10000000000;


            public void processBase(string contig, int offset, long mappedReadCount)
            {
                if (contig != currentContig || baseOffset + regionSize < offset)
                {
                    if (currentContig != "")
                    {
                        closeRegion();
                    }

                    currentContig = contig;
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
                    totalZForBasesWithBaselineExpression / (double)regionSize);

                nBasesExpressed = 0;
                nBasesExpressedWithBaselineExpression = 0;
                nBasesExpressedWithoutBaselineExpression = 0;
                totalReadsMappedToBasesWithoutBaselineExpression = 0;
                totalReadsMappedToBasesWithBaselineExpression = 0;
                totalZForBasesWithBaselineExpression = 0;
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
                }
            }

            static public void printHeader(StreamWriter outputFile)
            {
                outputFile.WriteLine("Contig\tContig Offset\tn Bases Expressed\tn Bases Expressed With Baseline Expression\tn Bases Expressed Without Baseline Expression\tTotal Reads Mapped To Bases With Baseline Expression\t" +
                    "Total Reads Mapped To Bases Without Baseline Expression\tCount of bases with baseline expression but not in this sample\tTotal z For Bases With Baseline Expression\tMin z For Bases With Baseline Expression\tMax z For Bases With BaselineExpression\t" +
                    "Mean z for Bases With Baseline Expression");
            }
        }

        static Dictionary<string, Dictionary<int, ExpressionTools.MeanAndStdDev>> expression = null;
        static Dictionary<string, string> sampleToParticipantIDMap = null;
        static Dictionary<string, ExpressionTools.Participant> participants = null;
        static List<ExpressionTools.Experiment> experiments = null;
        static Dictionary<string, ExpressionTools.Experiment> participantToExperimentMapping = null;
        static Dictionary<string, ExpressionTools.TCGARecord> tcgaRecords = null;
        static Dictionary<string, ExpressionTools.Sample> allSamples = null;
        static int regionSize;

        static void ProcessPeople(List<ExpressionTools.Experiment> experiments) {
            while (true) {

                ExpressionTools.Experiment experiment;
                lock(experiments) {
                    if (experiments.Count() == 0) {
                        //
                        // No more work, we're done.
                        //
                        return;
                    }
                    experiment = experiments[0];
                    experiments.RemoveAt(0);
                }


                //
                // Run through the allcount file
                //
                var allcountTimer = new Stopwatch();
                allcountTimer.Start();
                var allcountReader = new StreamReader(new GZipStream(new StreamReader(experiment.TumorRNAAnalysis.allcountFileName).BaseStream, CompressionMode.Decompress));

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
                    Console.WriteLine("Empty or corrupt allcount file " + experiment.TumorRNAAnalysis.allcountFileName);
                    continue;
                }

                if (line[headerBeginning.Count()] != '1')
                {
                    Console.WriteLine("Unsupported major version of allcount file " + experiment.TumorRNAAnalysis.allcountFileName + ".  Header line: " + line);
                    continue;
                }

                line = allcountReader.ReadLine();
                if (null == line)
                {
                    Console.WriteLine("Corrupt or tuncated allcount file " + experiment.TumorRNAAnalysis.allcountFileName);
                    continue;
                }

                var fields = line.Split(' ');
                if (fields.Count() != 16)
                {
                    Console.WriteLine("Corrupt or tuncated allcount file " + experiment.TumorRNAAnalysis.allcountFileName + ".  Second line has " + fields.Count() + " fields: " + line);
                    continue;
                }

                long mappedHQNuclearReads = 0;
                try
                {
                    mappedHQNuclearReads = Convert.ToInt64(fields[0]);
                }
                catch (FormatException)
                {
                    Console.WriteLine("Format exception parsing mapped HQ read count for file " + experiment.TumorRNAAnalysis.allcountFileName + " from line: " + line);
                    continue;
                }

                line = allcountReader.ReadLine();   // The blank line

                line = allcountReader.ReadLine();
                if (null == line) {
                    Console.WriteLine("Allcount file truncated before contig count: " + experiment.TumorRNAAnalysis.allcountFileName);
                    continue;
                }

                const string numContigsLineBeginning = "NumContigs: ";
                if (line.Count() < numContigsLineBeginning.Count() + 1 || line.Substring(0, numContigsLineBeginning.Count()) != numContigsLineBeginning) {
                    Console.WriteLine("Malformed NumContigs line in " + experiment.TumorRNAAnalysis.allcountFileName + ": " + line);
                    continue;
                }

                int numContigs = 0;
                try {
                    numContigs = Convert.ToInt32(line.Substring(numContigsLineBeginning.Count()));
                } catch (FormatException) {
                    Console.WriteLine("Couldn't parse NumContigs line in file " + experiment.TumorRNAAnalysis.allcountFileName + ": " + line);
                    continue;
                }
                
                if (numContigs < 1) {
                    Console.WriteLine("Invalid numContigs in " + experiment.TumorRNAAnalysis.allcountFileName + ": " + line);
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
                        Console.WriteLine("File truncated in contig list " + experiment.TumorRNAAnalysis.allcountFileName);
                        fileCorrupt = true;
                        break;
                    }

                    fields = line.Split('\t');
                    if (fields.Count() != 2) {
                        Console.WriteLine("Incorrect contig line format in file " + experiment.TumorRNAAnalysis.allcountFileName + ": " + line);
                        fileCorrupt = true;
                        break;
                    }

                    contigs[whichContig].name = fields[0].ToLower();
                    try {
                        contigs[whichContig].length = Convert.ToInt64(fields[1]);
                    } catch (FormatException) {
                        Console.WriteLine("Incorrect contig line format in file " + experiment.TumorRNAAnalysis.allcountFileName + ": " + line);
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

                int indexOfLastSlash = experiment.TumorRNAAnalysis.allcountFileName.LastIndexOf('\\');
                if (-1 == indexOfLastSlash) {
                    Console.WriteLine("Couldn't find a backslash in allcount pathname, which is supposed to be absolute: " + experiment.TumorRNAAnalysis.allcountFileName);
                    continue;
                }

                string directory = experiment.TumorRNAAnalysis.allcountFileName.Substring(0,indexOfLastSlash + 1);  // Includes trailing backslash
                var outputFilename = directory + experiment.TumorRNAAnalysis.analysis_id + ".regional_expression.txt";
                var outputFile = new StreamWriter(outputFilename);

                outputFile.WriteLine("RegionalExpression v1.0\t" + experiment.TumorRNAAnalysis.analysis_id + "\t" + experiment.TumorRNAAnalysis.allcountFileName + "\t" + regionSize);
                outputFile.WriteLine("NumContigs: " + numContigs);
                Region.printHeader(outputFile);

                Region region = new Region(regionSize, expression, mappedHQNuclearReads, outputFile);

                int currentOffset = -1;
                int currentMappedReadCount = -1;

                bool sawDone = false;
                string strippedContigName = "";

                while (null != (line = allcountReader.ReadLine()))
                {
                    if (sawDone)
                    {
                        Console.WriteLine("File " + experiment.TumorRNAAnalysis.allcountFileName + " continues after **done** line: " + line);
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
                        Console.WriteLine("Unexpected blank line in " + experiment.TumorRNAAnalysis.allcountFileName);
                        fileCorrupt = true;
                        break;
                    }

                    if (line[0] == '>')
                    {
                        whichContig++;
                        if (whichContig >= numContigs)
                        {
                            Console.WriteLine("Saw too many contigs in " + experiment.TumorRNAAnalysis.allcountFileName + ": " + line);
                            fileCorrupt = true;
                            break;
                        }

                        if (line.Substring(1).ToLower() != contigs[whichContig].name)
                        {
                            Console.WriteLine("Unexpected contig in " + experiment.TumorRNAAnalysis.allcountFileName + ".  Expected " + contigs[whichContig].name + ", got ", line.Substring(1));
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
                    Console.WriteLine("Truncated allcount file " + experiment.TumorRNAAnalysis.allcountFileName);
                    fileCorrupt = true;
                }

                region.closeRegion();
                outputFile.Close();
                if (fileCorrupt)
                {
                    File.Delete(outputFilename);
                }
                else
                {
                    allcountTimer.Stop();
                    Console.WriteLine("Processed " + experiment.participant.participantId + " in " + (allcountTimer.ElapsedMilliseconds + 500) / 1000 + "s");
                }
            }
        }

        static void Main(string[] args)
        {
            if (args.Count() < 3) {
                Console.WriteLine("usage: RegionalExpression expression_disease_file_directory regionSize <one or more participantIDs>");
                return;
            }

            string expressionDirectory = args[0];
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

            //
            // Add a trailing \ if it doesn't have one
            //
            if (expressionDirectory[expressionDirectory.Count() - 1] != '\\')
            {
                expressionDirectory = expressionDirectory + @"\";
            }

            tcgaRecords = ExpressionTools.LoadTCGARecords(null, null, @"\\gcr\scratch\b99\bolosky\tcga-all.xml");
            ExpressionTools.LoadTCGARecordsForLocalRealigns(tcgaRecords, null, @"\\gcr\scratch\b99\bolosky\realigns.txt");
            ExpressionTools.LoadTCGAAdditionalMetadata(tcgaRecords, @"\\bolosky\f$\sequence\reads\tcga\tcgaAdditionalMetadata.txt");
          
            participants = ExpressionTools.BuildParticipantData(tcgaRecords, out allSamples, @"\\gcr\scratch\b99\bolosky\clinical");

            sampleToParticipantIDMap = ExpressionTools.CreateSampleToParticipantIDMap(tcgaRecords);
            // do we really need these?  This is very slow and memory consumptive.  ExpressionTools.AddAllMAFFilesToParticipants(participants, sampleToParticipantIDMap, @"\\gcr\scratch\b99\bolosky\mafs\");
            experiments = ExpressionTools.LoadExperimentsFromFile(@"\\gcr\scratch\b99\bolosky\experiments.txt", participants, tcgaRecords);

            participantToExperimentMapping = ExpressionTools.BuildParticipantToExperimentMapping(experiments);

 
            string expressionDisease = "";

            var peopleByDisease = new Dictionary<string, List<ExpressionTools.Experiment>>();

            for (int i = 2; i < args.Count(); i++)  // for each person we're processing
            {
                ExpressionTools.Participant participant;
                if (!participants.ContainsKey(args[i]))
                {
                    Console.WriteLine("Couldn't find participant " + args[i] + ", are you sure it's a correct participant ID?");
                    return;
                }

                participant = participants[args[i]];
                var experiment = participantToExperimentMapping[participant.participantId];
                if (experiment.TumorRNAAnalysis == null || experiment.TumorRNAAnalysis.allcountFileName == null) {
                    Console.WriteLine("Participant " + participant.participantId + " doesn't appear to have an allcount file");
                    continue;
                }

                if (!peopleByDisease.ContainsKey(experiment.disease_abbr))
                {
                    peopleByDisease.Add(experiment.disease_abbr, new List<ExpressionTools.Experiment>());
                }

                peopleByDisease[experiment.disease_abbr].Add(experiment);
            }

            //
            // Now run a parallel phase for each disease.
            //
            foreach (var entry in peopleByDisease) {
                //
                // First, load the expression for that disease.
                //

                var experiment = entry.Value[0];
                Console.Write("Loading expression file for " + experiment.disease_abbr + "...");
                var timer = new Stopwatch();
                timer.Start();

                expression = null;  // Let the garbage collector get rid of the previous one while we're loading the next

                expression = ExpressionTools.LoadExpressionFile(expressionDirectory + "expression_" + experiment.disease_abbr);   // This can take forever!
                expressionDisease = experiment.disease_abbr;

                timer.Stop();
                Console.WriteLine((timer.ElapsedMilliseconds + 500) / 1000 + "s");

                int totalNumberOfExperiments = experiments.Count();
                timer.Reset();
                timer.Start();

                var threads = new List<Thread>();
                for (int i = 0; i < Environment.ProcessorCount; i++) {
                    threads.Add(new Thread(() => ProcessPeople(entry.Value)));
                }

                threads.ForEach(t => t.Start());
                threads.ForEach(t => t.Join());

                timer.Stop();
                Console.WriteLine("Processed " + totalNumberOfExperiments + " experiments in " + (timer.ElapsedMilliseconds + 500) / 1000 + "seconds");

            }
        }
    }
}
