using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Threading;

namespace FilterVCF
{
    class Program
    {
        public class WorkLeft
        {
            public int nQueued;
            public int nFinished;
        }

        class WorkItem
        {
            public WorkItem(StreamWriter outputFile_, StreamWriter scriptFile_, string disease_abbr_, string participantID_, string reference_assembly_, string normalDNAAnalysisID_,
            string normalDNA_, string normalRNAAnalysisID_, string normalRNA_, string tumorRNAAnalysisID_, string tumorRNA_, string filename_, WorkLeft workLeft_)
            {
                outputFile = outputFile_;
                scriptFile = scriptFile_;
                disease_abbr = disease_abbr_;
                participantID = participantID_;
                reference_assembly = reference_assembly_;
                normalDNAAnalysisID = normalDNAAnalysisID_;
                normalDNA = normalDNA_;
                normalRNAAnalysisID = normalRNAAnalysisID_;
                normalRNA = normalRNA_;
                tumorRNAAnalysisID = tumorRNAAnalysisID_;
                tumorRNA = tumorRNA_;
                filename = filename_;
                workLeft = workLeft_;
            }

            public void Execute()
            {
                ProcessOneFile(outputFile, scriptFile, disease_abbr, participantID, reference_assembly, normalDNAAnalysisID,
                    normalDNA, normalRNAAnalysisID, normalRNA, tumorRNAAnalysisID, tumorRNA, filename, workLeft);
            }

            StreamWriter outputFile;
            StreamWriter scriptFile;
            string disease_abbr;
            string participantID;
            string reference_assembly;
            string normalDNAAnalysisID; 
            string normalDNA;
            string normalRNAAnalysisID;
            string normalRNA;
            string tumorRNAAnalysisID;
            string tumorRNA;
            string filename;
            WorkLeft workLeft;
        }

        static void WorkerThread(object workQueue_)
        {
            List<WorkItem> workQueue = (List<WorkItem>)workQueue_;
            while (true)
            {
                WorkItem workItem;
                lock (workQueue)
                {
                    if (workQueue.Count() == 0)
                    {
                        return;
                    }

                    workItem = workQueue[0];
                    workQueue.RemoveAt(0);
                }

                workItem.Execute();
            }
        }

        static string DefaultExperimentsFilename = @"f:\temp\expression\experiments.txt";

        static void PrintUsage()
        {
            Console.WriteLine("usage: FilterVCF outputFilename outputScriptFileName {experimentsFile}");
            Console.WriteLine("ExperimentsFile is " + DefaultExperimentsFilename + " by default");
        }

        static void ProcessOneFile(StreamWriter outputFile, StreamWriter scriptFile, string disease_abbr, string participantID, string reference_assembly, string normalDNAAnalysisID, 
            string normalDNA, string normalRNAAnalysisID, string normalRNA, string tumorRNAAnalysisID, string tumorRNA, string filename, WorkLeft workLeft)
        {
            var reader = new StreamReader(filename);

            string line;
            while (null != (line = reader.ReadLine()))
            {
                if (line.Count() == 0 || line[0] == '#')
                {
                    //
                    // Blank or header line.
                    //
                    continue;
                }

                string[] fields = line.Split('\t');

                if (fields.Count() < 8) continue;

                //
                // Ignore variants with low quality.
                //
                double quality = Convert.ToDouble(fields[5]);
                if (quality < 30)
                {
                    continue;
                }

                //
                // Ignore variants that aren't simple SNPs.
                //
                if (fields[3].Count() != 1 || fields[4].Count() != 1 || fields[3] == "." || fields[4] == ".")
                {
                    continue;
                }

                //
                // Ignore homozygous variants (here defined as allele frequency between 1/3 and 2/3).
                //
                string[] info = fields[7].Split(';');
                var afs = info.Where(s => s.Count() > 3 & s[0] == 'A' && s[1] == 'F' && s[2] == '=');
                if (afs.Count() != 1)
                {
                    Console.WriteLine("Couldn't find AF field in " + line);
                    continue;
                }

                double af = 0;
                foreach (var afstring in afs)
                {
                    //
                    // There is exactly one of these.
                    //
                    af = Convert.ToDouble(afstring.Substring(3));
                }

                if (af < .333 || af > .667)
                {
                    continue;
                }

                // Find if it's in the range we care about.

                if (fields[0] != "9" && fields[0].ToLower() != "chr9")
                {
                    continue;
                }

                int pos = Convert.ToInt32(fields[1]);
                if (!((pos >= 21967750 && pos <= 21971207) || (pos >= 21968573 && pos <= 21968770) || (pos >= 21970900 && pos <= 21971207) || (pos >= 21974402 && pos <= 21974826)))
                {
                    continue;
                }

                string chromPrefix;
                if (reference_assembly == "NCBI36_BCCAGSC_variant".ToLower() || reference_assembly == "grch37-lite" || reference_assembly == "grch37" || reference_assembly == "hg19_broad_variant" 
                    || reference_assembly == "grch37-lite_wugsc_variant_1" || reference_assembly == "grch37-lite_wugsc_variant_2"
                    || reference_assembly == "grch37-lite-+-hpv_redux-build" || reference_assembly == "HS37D5".ToLower() || reference_assembly == "NCBI36_WUGSC_variant".ToLower())
                {
                    chromPrefix = " ";
                }
                else
                {
                    chromPrefix = " chr";
                }

                string extractedReadsFilenameBase = disease_abbr + @"\" + participantID + "-CDKN2A-chr9-" + (pos - 100) + "-" + (pos + 100) + "-RNA-";
                string extractedReadsTumorFilename = extractedReadsFilenameBase + "tumor";
                string extractedReadsNormalFilename;
                lock (scriptFile)
                {
                    if ("" == normalRNA)
                    {
                        extractedReadsNormalFilename = "";
                    }
                    else
                    {
                        extractedReadsNormalFilename = extractedReadsFilenameBase + "normal";
                        scriptFile.WriteLine(@"samtools view " + normalRNA + " " + chromPrefix + @"9:" + (pos - 100) + "-" + (pos + 100) + " > " + extractedReadsNormalFilename);
                    }
                    scriptFile.WriteLine(@"samtools view " + tumorRNA + " " + chromPrefix + @"9:" + (pos - 100) + "-" + (pos + 100) + " > " + extractedReadsTumorFilename);

                    outputFile.WriteLine(participantID + "\t" + disease_abbr + "\t" + reference_assembly + "\t" + normalDNAAnalysisID + "\t" + normalDNA + "\t" + normalRNAAnalysisID + "\t" + normalRNA + "\t" + tumorRNAAnalysisID + "\t" + tumorRNA + "\t" +
                        extractedReadsNormalFilename + "\t" + extractedReadsTumorFilename + "\t" + line);
                }
            }

            lock (workLeft)
            {
                workLeft.nFinished++;
                Console.WriteLine("" + workLeft.nFinished + ":" + filename);
                if (workLeft.nFinished == workLeft.nQueued)
                {
                    Monitor.PulseAll(workLeft);
                }
            }
        }

        static void Main(string[] args)
        {

            string experimentsFilename = DefaultExperimentsFilename;
            if (args.Count() != 2 && args.Count() != 3)
            {
                PrintUsage();
                return;
            }

            if (args.Count() == 3)
            {
                experimentsFilename = args[2];
            }

            StreamWriter outputFile = new StreamWriter(args[0]);
            StreamWriter scriptFile = new StreamWriter(args[1]);


            outputFile.WriteLine("ParticipantID\tDisease Abbr\tReference Assembly\tNormal DNA Analysis ID\tNormal DNA Filename\tNormal RNA Analysis ID\tNormal RNA Filename\tTumor RNA AnanlysisID\tTumor RNA Filename\tExtracted Reads Normal Filename\tExtracted Reads Tumor Filename\tChrom\tPos\tID\tRef\tAlt\tQual\tFilter\tInfo\tFormat");
            string[] experiments = File.ReadAllLines(experimentsFilename);

            WorkLeft workLeft = new WorkLeft();

            workLeft.nQueued = experiments.Count();
            workLeft.nFinished = 0;

            List<WorkItem> workQueue = new List<WorkItem>();

            foreach (var experiment in experiments)
            {
                string[] fields = experiment.Split('\t');
                if (fields.Count() < 12 || fields[0] == "disease_abbr" || fields[11] == "")
                {
                    lock (workLeft)
                    {
                        workLeft.nQueued--;
                    }
                    continue;
                }

                workQueue.Add(new WorkItem(outputFile, scriptFile, fields[0], fields[2], fields[1], fields[5], fields[9], fields[6], fields[10], fields[3], fields[7], fields[11], workLeft));
            }

            List<Thread> threads = new List<Thread>();

            for (int i = 0; i < 15; i++)
            {
                var start = new ParameterizedThreadStart(WorkerThread);
                Thread thread = new Thread(start);
                threads.Add(thread);
                thread.Start(workQueue);
            }

            lock (workLeft)
            {
                while (workLeft.nFinished < workLeft.nQueued)
                {
                    Monitor.Wait(workLeft);
                }
            }

            outputFile.Close();
            scriptFile.Close();

        }
    }
}
