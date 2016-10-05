using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ExpressionLib;
using System.IO;
using System.Threading;
using System.Diagnostics;


namespace SelectGermlineVariants
{
    class Program
    {

        const long granularity = 10000; // How often to select a variant if one is available.
        const long isolationDistance = 150; // How many bases around a selected variant must match the germline exclusively.

        class CandidateVariant
        {
            public CandidateVariant(string line_, long locus_, double odds_)
            {
                line = line_;
                locus = locus_;
                odds = odds_;
            }

            public string line;
            public long locus;
            public double odds;
        }

        static void EmitBestCandidate(StreamWriter outputFile, string chromosome, List<CandidateVariant> liveCandidates)
        {
            if (liveCandidates.Count() == 0)
            {
                return;
            }

            CandidateVariant bestCandiate = liveCandidates[0];

            for (int i = 0; i < liveCandidates.Count(); i++) // Yes, I know there's some clever c# way to do this with linq, but it would take me longer to look up than to write this loop
            {
                if (liveCandidates[i].odds > bestCandiate.odds)
                {
                    bestCandiate = liveCandidates[i];
                }
            }

            outputFile.WriteLine(chromosome + "\t" + bestCandiate.locus + "\t" + bestCandiate.line);
        }

        static void EliminateTooNearCandidates(List<CandidateVariant> liveCandidates, long locus)
        {
            var candidatesToEliminate = new List<CandidateVariant>();
            foreach (var candidate in liveCandidates)
            {
                if (candidate.locus + isolationDistance >= locus)
                {
                    candidatesToEliminate.Add(candidate);
                }
            }

            foreach (var candidate in candidatesToEliminate)
            {
                liveCandidates.Remove(candidate);
            }
        }

        static void ProcessRuns(List<string> workItems)
        {
            while (true)
            {
                string vcfPathname;
                lock (workItems)
                {
                    if (workItems.Count() == 0)
                    {
                        return;
                    }

                    vcfPathname = workItems[0];
                    workItems.RemoveAt(0);

                    Console.WriteLine("" + workItems.Count() + " vcf" + (workItems.Count() == 1 ? " remains" : "s remain") + " queued.");
                }

                StreamReader vcfFile = null;
                try
                {
                    vcfFile = new StreamReader(vcfPathname);
                }
                catch (FileNotFoundException)
                {
                    Console.WriteLine("File not found on vcf " + vcfPathname + ".  Skipping.");
                    continue;
                }

                string line;
                while (null != (line = vcfFile.ReadLine()) && line.Count() != 0 && line[0] == '#') {
                    // Skip the header lines
                }

                if (null == line || line.Count() == 0) {
                    Console.WriteLine("Corrupt vcf: missing body: " + vcfPathname);
                    continue;
                }

                var outputFilename = vcfPathname.Substring(0, vcfPathname.LastIndexOf('.')) + ".selectedVariants";
                var outputFile = new StreamWriter(outputFilename);
                outputFile.WriteLine("SelectGermlineVariants v1.0 for input file " + vcfPathname);

                string currentChromosome = "";

                bool badFile = false;

                var liveCandidates = new List<CandidateVariant>();
                var previousGrainsCandidates = new List<CandidateVariant>();
                long lastLocus = -isolationDistance - 1;

                while (null != (line = vcfFile.ReadLine()))
                {
                    var fields = line.Split('\t');

                    if (fields.Count() != 10)
                    {
                        Console.WriteLine("Wrong number of fields (" + fields.Count() + " != 10) in vcf line: '" + line + "' in file " + vcfPathname + ".  Ignoring file.");

                        badFile = true;
                        break;
                    }

                    var infoFields = fields[7].Split(';');

                    var info = new Dictionary<string, string>();
                    foreach (var infoField in infoFields)
                    {
                        var keyValue = infoField.Split('=');
                        if (keyValue.Count() != 2)
                        {
                            Console.WriteLine("Unable to parse info field '" + infoField + " in file " + vcfPathname);
                            badFile = true;
                            break;
                        }

                        info.Add(keyValue[0], keyValue[1]);
                    }

                    if (!info.ContainsKey("AN") || !info.ContainsKey("AC") || !info.ContainsKey("CIGAR") || !info.ContainsKey("DP") || !info.ContainsKey("AF") || !info.ContainsKey("AB") || !info.ContainsKey("ODDS"))
                    {
                        Console.WriteLine("vcf line '" + line + " doesn't contain one or more required info fields.  Skipping file " + vcfPathname);
                        badFile = true;
                    }

                    if (badFile)
                    {
                        break;
                    }

                    bool goodCandidate = true;
                    double alleleFrequency = 0;
                    double alleleBalance = 0;
                    int alleleCount = 0;
                    int alleleNumber = 0;
                    double odds = 0;
                    string cigar = info["CIGAR"];
                    long locus = 0;

                    try
                    {
                        if (info["AF"].Contains(',')) // This happens with multiple alleles (typically, heterozygous with both different from the reference).  Don't try to parse them, just don't select it
                        {
                            goodCandidate = false;
                        }
                        else
                        {
                            alleleFrequency = Convert.ToDouble(info["AF"]);
                            alleleBalance = Convert.ToDouble(info["AB"]);
                            alleleCount = Convert.ToInt32(info["AC"]);
                            alleleNumber = Convert.ToInt32(info["AN"]);
                            odds = Convert.ToDouble(info["ODDS"]);
                        }
                        locus = Convert.ToInt64(fields[1]);

                    }
                    catch (FormatException)
                    {
                        Console.WriteLine("Error parsing info fields in line " + line + " of file " + vcfPathname + ".  Skipping file.");
                        badFile = true;
                        break;
                    }

                    if (fields[0] == currentChromosome && locus < lastLocus)
                    {
                        Console.WriteLine("out-of-order variant " + line + " in file " + vcfPathname + ". Skipping file.");
                        badFile = true;
                        break;
                    }

                    //
                    // Figure out if we've moved into anothe grain, in which case we emit the best candidate we've got for the previous grain.
                    //
                    if (currentChromosome != fields[0] || locus / granularity != lastLocus / granularity)
                    {
                        EmitBestCandidate(outputFile, currentChromosome, previousGrainsCandidates);
                        previousGrainsCandidates = liveCandidates;
                        liveCandidates = new List<CandidateVariant>();

                        if (currentChromosome != fields[0])
                        {
                            //
                            // Starting a new chromosome, so we don't need to hang on to an old
                            // chromosome's candidates to make sure that we don't have any variants too close to
                            // the end of the old grain.
                            //
                            EmitBestCandidate(outputFile, currentChromosome, previousGrainsCandidates);
                            previousGrainsCandidates = new List<CandidateVariant>();
                            currentChromosome = fields[0];
                            lastLocus = -isolationDistance - 1;
                        }
                    }

                    EliminateTooNearCandidates(liveCandidates, locus);
                    EliminateTooNearCandidates(previousGrainsCandidates, locus);

                    goodCandidate = goodCandidate && alleleFrequency == 0.5 && alleleBalance > 0.4 && alleleBalance < 0.6 && alleleCount == 1 && alleleNumber == 2 && cigar == "1X" && odds > 20 && lastLocus + isolationDistance < locus;

                    if (goodCandidate)
                    {
                        liveCandidates.Add(new CandidateVariant(line, locus, odds));
                    }

                    lastLocus = locus;
                } // While we have a VCF line

                if (badFile) {
                    outputFile.Close();
                    File.Delete(outputFilename);
                    continue;
                }
                else
                {
                    EmitBestCandidate(outputFile, currentChromosome, previousGrainsCandidates);
                    EmitBestCandidate(outputFile, currentChromosome, liveCandidates);
                    outputFile.WriteLine("**done**");
                    outputFile.Close();
                }

                vcfFile.Close();
            }
        }

        static void Main(string[] args)
        {
            if (args.Count() == 0)
            {
                Console.WriteLine("usage: SelectGermlineVariants <participant ID>");
                return;
            }
                
            var vcfPathnameByParticipantId = new Dictionary<string, string>();
            List<string> workItems = new List<string>();

            var experimentsFile = new StreamReader(@"\\gcr\scratch\b99\bolosky\experiments.txt");
            string line;

            while (null != (line = experimentsFile.ReadLine()))
            {
                var fields = line.Split('\t');

                vcfPathnameByParticipantId.Add(fields[2], fields[11]);
            }
            
            experimentsFile.Close();
            experimentsFile = null;

            int nValidParticipants = 0;

            foreach (var arg in args)
            {
                if (!vcfPathnameByParticipantId.ContainsKey(arg))
                {
                    Console.WriteLine(arg + " does not appear to be a participant ID; ignoring.");
                    continue;
                }

                if (vcfPathnameByParticipantId[arg] == "")
                {
                    Console.WriteLine(arg + " doesn't appear to have a vcf yet.  Ignoring.");
                    continue;
                }

                workItems.Add(vcfPathnameByParticipantId[arg]);
                nValidParticipants++;
            }


            var timer = new Stopwatch();
            timer.Start();

            var threads = new List<Thread>();
            for (int i = 0; i < Environment.ProcessorCount; i++)
            {
                threads.Add(new Thread(() => ProcessRuns(workItems)));
            }

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            timer.Stop();
            Console.WriteLine("Processed " + nValidParticipants + " participants in " + (timer.ElapsedMilliseconds + 500) / 1000 + "s");
        }
    }
}
