using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.Diagnostics;
using System.Threading;

namespace ExtractASVsForHighASERegions
{
    class Program
    {

        static List<ASETools.ASEMapLine> normalMapLinesToProcess, tumorMapLinesToProcess;

        class AnnotatedVariantAndCase
        {
            public readonly ASETools.AnnotatedVariant annotatedVariant;
            public readonly string caseId;
            public readonly bool normal;
            public readonly bool tumor; // Recall, this is an annotated variant, which may be differently valid for normal or tumor.  Hence, two bools.

            public AnnotatedVariantAndCase(ASETools.AnnotatedVariant annotatedVariant_, string caseId_, bool normal_, bool tumor_)
            {
                annotatedVariant = annotatedVariant_;
                caseId = caseId_;
                normal = normal_;
                tumor = tumor_;
            }

            static public int Compare(AnnotatedVariantAndCase a, AnnotatedVariantAndCase b)
            {
                return ASETools.AnnotatedVariant.CompareByLocus(a.annotatedVariant, b.annotatedVariant);
            }
        } // AnnotatedVariantAndCase

        static int nCasesProcessed = 0;

        static ASETools.Configuration configuration;
        static ASETools.GeneMap geneMap;
        static Dictionary<string, ASETools.ASEMapPerGeneLine> perGeneASEMap;

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            configuration = ASETools.Configuration.loadFromFile(args);

            if (null == configuration)
            {
                Console.WriteLine("Giving up because we were unable to load configuration.");
                return;
            }

            if (configuration.commandLineArgs.Count() > 1)
            {
                Console.WriteLine("usage: ExtractASVsForHighASERegions {ASE cutoff}");
                return;
            }

            double cutoff = 0.5;
            if (configuration.commandLineArgs.Count() == 1)
            {
                cutoff = Convert.ToDouble(configuration.commandLineArgs[0]);    // Yes, I should catch the exception and print the usage message.
            }

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

            if (null == cases)
            {
                Console.WriteLine("Unable to load cases.");
            }

            perGeneASEMap = ASETools.ASEMapPerGeneLine.ReadFromFileToDictionary(configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);

            if (null == perGeneASEMap)
            {
                Console.WriteLine("You must first create the per-gene ASE map in " + configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);
                return;
            }

            var geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(configuration.geneLocationInformationFilename));
            geneMap = new ASETools.GeneMap(geneLocationInformation.genesByName);

            var allMapLines = ASETools.ASEMapLine.ReadFile(configuration.finalResultsDirectory + ASETools.ASEMapFilename);
            if (null == allMapLines)
            {
                return; // It already printed an error message
            }

            normalMapLinesToProcess = allMapLines.Where(x => !x.tumor && x.mean >= cutoff).ToList();
            tumorMapLinesToProcess = allMapLines.Where(x => x.tumor && x.mean >= cutoff).ToList();

            var casesToProcess = cases.Select(x => x.Value).Where(x => x.annotated_selected_variants_filename != "").ToList();
            int nCasesToProcess = casesToProcess.Count();

            var results = new List<AnnotatedVariantAndCase>();

            Console.Write("Processing " + nCasesToProcess + " cases, 1 dot/100 cases: ");

            var threads = new List<Thread>();
            for (int i = 0; i < 1 /*Environment.ProcessorCount*/; i++)
            {
                threads.Add(new Thread(() => WorkerThread(casesToProcess, results)));
            }

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            results.Sort(AnnotatedVariantAndCase.Compare);  // So they come out in the output file in order

            var outputFile = ASETools.CreateStreamWriterWithRetry(@"f:\temp\annotated_variants_in_high_ASE_regions.txt");   // Really, don't use hardcoded filenames.
            outputFile.WriteLine("Case id\ttumor\tnormal\t" + ASETools.AnnotatedVariant.headerString);

            foreach (var variant in results)
            {
                outputFile.WriteLine(variant.caseId + "\t" + variant.tumor + "\t" + variant.normal + "\t" + variant.annotatedVariant.toString());
            }

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine();
            Console.WriteLine("Processed " + nCasesToProcess + " cases in " + ASETools.ElapsedTimeInSeconds(timer));
        } // Main

        static void WorkerThread(List<ASETools.Case> casesToProcess, List<AnnotatedVariantAndCase> results)
        {
            var localResults = new List<AnnotatedVariantAndCase>();
            var localNormalMapLinesToProcess = normalMapLinesToProcess.Where(x => !x.tumor).ToList();    // Cheap deep copy
            var localTumorMapLinesToProcess = normalMapLinesToProcess.Where(x => x.tumor).ToList();    // Cheap deep copy

            while (true)
            {
                ASETools.Case case_;
                lock (casesToProcess)
                {
                    if (casesToProcess.Count() == 0)
                    {
                        results.AddRange(localResults);
                        return;
                    }

                    case_ = casesToProcess[0];
                    casesToProcess.RemoveAt(0);

                    nCasesProcessed++;
                    if (nCasesProcessed % 100 == 0)
                    {
                        Console.Write(".");
                    }
                } // lock

                var variants = ASETools.AnnotatedVariant.readFile(case_.annotated_selected_variants_filename);

                var copyNumber = ASETools.CopyNumberVariation.ReadBothFiles(case_);
 
                bool tumor = false;
                bool normal = false;

                foreach (var variant in variants)
                {
                    normal = variant.IsASECandidate(false, copyNumber, configuration, perGeneASEMap, geneMap) && localNormalMapLinesToProcess.Where(x => x.chromosome == variant.contig && x.locus <= variant.locus && x.locus + 1000000 > variant.locus).Count() > 0;
                    tumor = variant.IsASECandidate(true, copyNumber, configuration, perGeneASEMap, geneMap) && localTumorMapLinesToProcess.Where(x => x.chromosome == variant.contig && x.locus <= variant.locus && x.locus + 1000000 > variant.locus).Count() > 0;

                    if (!normal && !tumor)
                    {
                        continue;
                    }

                    localResults.Add(new AnnotatedVariantAndCase(variant, case_.case_id, normal, tumor));
                }
            } // while true
        }
    }
}
