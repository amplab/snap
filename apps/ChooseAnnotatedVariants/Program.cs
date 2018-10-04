using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;
using System.Diagnostics;

namespace ChooseAnnotatedVariants
{
    class Program
    {
        static ASETools.Configuration configuration;
        static ASETools.ASERepetitiveRegionMap repetitiveRegionMap;
        static PerThreadState globalState = new PerThreadState();
        static ASETools.GeneLocationsByNameAndChromosome geneLocationInformation;
        static ASETools.GeneMap geneMap;
        static Dictionary<string, ASETools.ASEMapPerGeneLine> perGeneASEMap;

        static void Main(string[] args)
        {
            Stopwatch timer = new Stopwatch();
            timer.Start();

            configuration = ASETools.Configuration.loadFromFile(args);
            if (null == configuration)
            {
                Console.WriteLine("Unable to load configuration.");
                return;
            }

            geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(configuration.geneLocationInformationFilename));
            geneMap = new ASETools.GeneMap(geneLocationInformation.genesByName);

            perGeneASEMap = ASETools.ASEMapPerGeneLine.ReadFromFileToDictionary(configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);

            if (null == perGeneASEMap)
            {
                Console.WriteLine("You must first create the per-gene ASE map in " + configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);
                return;
            }

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);
            if (null == cases)
            {
                Console.WriteLine("Unable to load cases");
                return;
            }

            repetitiveRegionMap = ASETools.ASERepetitiveRegionMap.loadFromFile(configuration.redundantChromosomeRegionFilename);
            if (null == repetitiveRegionMap)
            {
                Console.WriteLine("Unable to load repetitive region map.");
                return;
            }

            var listOfCases = cases.Select(_ => _.Value).Where(_ => _.tentative_annotated_selected_variants_filename != "").ToList();
            if (listOfCases.Count() != cases.Count())
            {
                Console.WriteLine("At least one case is missing tentative annotated selected variants.");
                //off for debugging. return;
            }

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Loading tentative variants for", "cases", listOfCases.Count(), out nPerDot);
            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, PerThreadState>(listOfCases, HandleOneCase, FinishUp, null, nPerDot);
            threading.run();
            Console.WriteLine();

            var summaryFilename = configuration.finalResultsDirectory + ASETools.variantSelectionSummaryFilename;
            var summaryFile = ASETools.CreateStreamWriterWithRetry(summaryFilename);
            if (null == summaryFile)
            {
                Console.WriteLine("Unable to open " + summaryFilename + ".  Skipping.");
            } else
            {
                var depthHistogram = new ASETools.PreBucketedHistogram(0, 1000, 1);
                var nearbyHistogram = new ASETools.PreBucketedHistogram(0, 1000, 1);

                foreach (var chr in globalState.variantMap.Select(_ => _.Key).ToList())
                {
                    foreach (var roundedLocus in globalState.variantMap[chr].Select(_ => _.Key).ToList())
                    {
                        foreach (var locus in globalState.variantMap[chr][roundedLocus].Select(_ => _.Key))
                        {
                            depthHistogram.addValue(globalState.variantMap[chr][roundedLocus][locus].Count());
                            int nNearby = 0;
                            for (int candidateRoundedLocus = roundedLocus - treeGranularity; candidateRoundedLocus <= roundedLocus + treeGranularity; candidateRoundedLocus += treeGranularity)
                            {
                                if (!globalState.variantMap[chr].ContainsKey(candidateRoundedLocus))
                                {
                                    continue;
                                }

                                nNearby += globalState.variantMap[chr][candidateRoundedLocus].Where(_ => Math.Abs(_.Key - locus) <= proximity).Select(_ => _.Value.Count()).Sum();
                            } // candidate rounded locus
                            nearbyHistogram.addValue(nNearby);
                        } // locus
                    } // rounded locus
                } // chr

                summaryFile.WriteLine("Distribution of tentative annotated variants by count at the same locus.");
                summaryFile.WriteLine(ASETools.HistogramResultLine.Header());
                depthHistogram.ComputeHistogram().ToList().ForEach(_ => summaryFile.WriteLine(_.ToString()));

                summaryFile.WriteLine("Distribution of tentative annotated variants by count within " + proximity + " of variant.");
                summaryFile.WriteLine(ASETools.HistogramResultLine.Header());
                nearbyHistogram.ComputeHistogram().ToList().ForEach(_ => summaryFile.WriteLine(_.ToString()));

                summaryFile.WriteLine("**done**");
            }

            ASETools.PrintMessageAndNumberBar("Selecting and writing variants for", "cases", listOfCases.Count(), out nPerDot);
            var threading2 = new ASETools.WorkerThreadHelper<ASETools.Case, int>(listOfCases, SelectVariantsForOneCase, null, null, nPerDot);
            threading2.run();

            Console.WriteLine();
            Console.WriteLine("Processed " + listOfCases.Count() + " cases in " + ASETools.ElapsedTimeInSeconds(timer));
        }

        const int treeGranularity = 100;
        const int proximity = 50;

        class PerThreadState
        {
            public Dictionary<string, Dictionary<int, Dictionary<int, List<double>>>> variantMap = new Dictionary<string, Dictionary<int, Dictionary<int, List<double>>>>(); // chr->locus rounded down to tree granularity->locus->VAFs

            public PerThreadState()
            {
                for (int chr = 1; chr <= ASETools.nHumanAutosomes; chr++)
                {
                    variantMap.Add("chr" + chr, new Dictionary<int, Dictionary<int, List<double>>>());
                }

                variantMap.Add("chrX", new Dictionary<int, Dictionary<int, List<double>>>());
                variantMap.Add("chrY", new Dictionary<int, Dictionary<int, List<double>>>());
            }
        }

        static void FinishUp(PerThreadState perThreadState)
        {
            lock (globalState)
            {
                foreach (var chr in perThreadState.variantMap.Select(_ => _.Key).ToList())
                {
                    foreach (var roundedLocus in perThreadState.variantMap[chr].Select(_ => _.Key).ToList())
                    {
                        if (!globalState.variantMap[chr].ContainsKey(roundedLocus))
                        {
                            globalState.variantMap[chr].Add(roundedLocus, new Dictionary<int, List<double>>());
                        }

                        foreach (var locus in perThreadState.variantMap[chr][roundedLocus].Select(_ => _.Key).ToList())
                        {
                            if (!globalState.variantMap[chr][roundedLocus].ContainsKey(locus))
                            {
                                globalState.variantMap[chr][roundedLocus].Add(locus, new List<double>());
                            }

                            globalState.variantMap[chr][roundedLocus][locus].AddRange(perThreadState.variantMap[chr][roundedLocus][locus]);
                        } // locus
                    } // rounded locus
                } // chr
            } // lock

            perThreadState.variantMap = null;
        } // FinishUp

        static void HandleOneCase(ASETools.Case case_, PerThreadState perThreadState)
        {
            var annotatedVariants = ASETools.AnnotatedVariant.readFile(case_.tentative_annotated_selected_variants_filename);

            if (annotatedVariants == null)
            {
                Console.WriteLine("Unable to load tentative annotated variants from " + case_.tentative_annotated_selected_variants_filename);
                return;
            }

            var copyNumber = ASETools.CopyNumberVariation.ReadBothFiles(case_);
            if (null == copyNumber)
            {
                Console.WriteLine("Unable to load copy number for case " + case_.case_id);
                return;
            }

            foreach (var annotatedVariant in annotatedVariants)
            {
                string whyNot;
                if (!annotatedVariant.IsASECandidate(out whyNot, true, copyNumber, configuration, perGeneASEMap, geneMap))
                {
                    continue;
                }

                int roundedLocus = (annotatedVariant.locus / treeGranularity) * treeGranularity;

                if (!perThreadState.variantMap[annotatedVariant.contig].ContainsKey(roundedLocus))
                {
                    perThreadState.variantMap[annotatedVariant.contig].Add(roundedLocus, new Dictionary<int, List<double>>());
                }

                if (!perThreadState.variantMap[annotatedVariant.contig][roundedLocus].ContainsKey(annotatedVariant.locus))
                {
                    perThreadState.variantMap[annotatedVariant.contig][roundedLocus].Add(annotatedVariant.locus, new List<double>());
                }
                perThreadState.variantMap[annotatedVariant.contig][roundedLocus][annotatedVariant.locus].Add(annotatedVariant.GetTumorAltAlleleFraction());
            } // variant
        } // HandleOneCase

        static void SelectVariantsForOneCase(ASETools.Case case_, int state)
        {

        }
    }
}
