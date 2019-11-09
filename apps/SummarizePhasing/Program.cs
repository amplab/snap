using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;

namespace SummarizePhasing
{
    class Program
    {
        static ASETools.CommonData commonData;

        const int maxVariantPairsToTrack = 100;

        class WorkerThreadState
        {
            public Dictionary<int, ASETools.PreBucketedHistogram> distributionOfPhasedPairsBySemanticMutationCount = new Dictionary<int, ASETools.PreBucketedHistogram>();
            public ASETools.PreBucketedHistogram distributionOfPhasedPairs = new ASETools.PreBucketedHistogram(0, maxVariantPairsToTrack, 1);
            public int totalPhasedVariantPairs = 0;
            public Dictionary<string, int> countOfPhasedPairsByChromosome = new Dictionary<string, int>();

            public Dictionary<string, int> countOfPhasedPairsByHugoSymbol = new Dictionary<string, int>();
            public Dictionary<string, int> countOfInconsistentPhasedPairsByHugoSymbol = new Dictionary<string, int>();
            public Dictionary<string, int> countOfPhasedPairsWhereDNAandRNAdisagreeByHugoSymbol = new Dictionary<string, int>();

            public List<string> genesAndCasesWithDNAandRNAdisagreeing = new List<string>();

            public void merge(WorkerThreadState peer)
            {
                for (int somaticMutationCount = 0; somaticMutationCount < 3; somaticMutationCount++)
                {
                    distributionOfPhasedPairsBySemanticMutationCount[somaticMutationCount].merge(peer.distributionOfPhasedPairsBySemanticMutationCount[somaticMutationCount]);
                }

                distributionOfPhasedPairs.merge(peer.distributionOfPhasedPairs);
                totalPhasedVariantPairs += peer.totalPhasedVariantPairs;

                ASETools.MergeDictionaries(countOfPhasedPairsByChromosome, peer.countOfPhasedPairsByChromosome);
                ASETools.MergeDictionaries(countOfPhasedPairsByHugoSymbol, peer.countOfPhasedPairsByHugoSymbol);
                ASETools.MergeDictionaries(countOfInconsistentPhasedPairsByHugoSymbol, peer.countOfInconsistentPhasedPairsByHugoSymbol);
                ASETools.MergeDictionaries(countOfPhasedPairsWhereDNAandRNAdisagreeByHugoSymbol, peer.countOfPhasedPairsWhereDNAandRNAdisagreeByHugoSymbol);
                genesAndCasesWithDNAandRNAdisagreeing.AddRange(peer.genesAndCasesWithDNAandRNAdisagreeing);
            }

            public WorkerThreadState()
            {
                for (int somaticMutationCount = 0; somaticMutationCount < 3; somaticMutationCount++)
                {
                    distributionOfPhasedPairsBySemanticMutationCount.Add(somaticMutationCount, new ASETools.PreBucketedHistogram(0, maxVariantPairsToTrack, 1));
                }

                ASETools.chromosomes.ToList().ForEach(chr => countOfPhasedPairsByChromosome.Add(chr, 0));
            } // ctor
        } // WorkerThreadState

        static WorkerThreadState globalState = new WorkerThreadState();

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);

            if (null == commonData)
            {
                return;
            }


            if (commonData.listOfCases.Any(_ => _.variant_phasing_filename == ""))
            {
                Console.WriteLine("Missing some input files, so the summary will only be partial.");
            }

            var casesToProcess = commonData.listOfCases.Where(_ => _.variant_phasing_filename != "").ToList();
            int nPerDot;

            ASETools.PrintMessageAndNumberBar("Processing", "cases", casesToProcess.Count(), out nPerDot);
            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, WorkerThreadState>(casesToProcess, HandleOneCase, FinishUp, null, nPerDot);
            threading.run();

            var outputFilename = commonData.configuration.finalResultsDirectory + ASETools.PhasingForNearbyVariantsFilename;
            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
            if (null == outputFile)
            {
                Console.WriteLine("Unable to open outputfile " + outputFilename);
                return;
            }

            outputFile.WriteLine("Total of " + globalState.totalPhasedVariantPairs + " phased variant pairs across " + casesToProcess.Count() + " tumors.");
            outputFile.WriteLine();

            outputFile.WriteLine("Distribution of phased pairs by tumor.");
            globalState.distributionOfPhasedPairs.WriteHistogram(outputFile);

            outputFile.WriteLine();
            outputFile.WriteLine("CDF of count of phased pairs per tumor by somatic mutation count");
            ASETools.PreBucketedHistogram.WriteBatchOfHistogramCDFs(outputFile, globalState.distributionOfPhasedPairsBySemanticMutationCount.Select(_ => new KeyValuePair<string, ASETools.PreBucketedHistogram>(_.Key.ToString(), _.Value)).ToList());

            outputFile.WriteLine();
            outputFile.WriteLine("Count of phased pairs by chromosome");
            outputFile.WriteLine("Chromosome\tCount");
            ASETools.chromosomes.ToList().ForEach(chr => outputFile.WriteLine(chr + "\t" + globalState.countOfPhasedPairsByChromosome[chr]));

            outputFile.WriteLine();
            outputFile.WriteLine("Count of phased pairs by gene");
            outputFile.WriteLine("Hugo Symbol\tcount\tCount where DNA and RNA are not internally consistent\tCount where DNA and RNA are internally consistent, but disagree with each other");
            var hugoSymbolsAndCounts = globalState.countOfPhasedPairsByHugoSymbol.ToList();
            hugoSymbolsAndCounts.Sort((a, b) => a.Key.CompareTo(b.Key));
            hugoSymbolsAndCounts.ForEach(_ => outputFile.WriteLine(_.Key + "\t" + _.Value + "\t" + globalState.countOfInconsistentPhasedPairsByHugoSymbol[_.Key] + "\t" + globalState.countOfPhasedPairsWhereDNAandRNAdisagreeByHugoSymbol[_.Key]));
            
            globalState.genesAndCasesWithDNAandRNAdisagreeing.Sort();
            outputFile.WriteLine();
            outputFile.WriteLine("Genes and cases where DNA and RNA are internally consistent but disagree with one another");
            outputFile.WriteLine("Hugo Symbol\tN somatic variants\tchrom\tlocus0\tlocus1\tDNA ref-ref/ref-alt\tDNA alt-ref/alt-alt\tRNA ref-ref/ref-alt\tRNA alt-ref/alt-alt\t" + 
                                 "Case ID\tNormal DNA filename\tTumor DNA filename\tTumor RNA filename");
            globalState.genesAndCasesWithDNAandRNAdisagreeing.ForEach(_ => outputFile.WriteLine(_));


            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine("Processed " + casesToProcess.Count() + " cases in " + ASETools.ElapsedTimeInSeconds(commonData.timer));

        } // Main

        static void HandleOneCase(ASETools.Case case_, WorkerThreadState state)
        {
            var variantPhasing = ASETools.VariantPairPhasing.readFromFile(case_.variant_phasing_filename);

            Enumerable.Range(0, 3).ToList().ForEach(nSomaticMutations => state.distributionOfPhasedPairsBySemanticMutationCount[nSomaticMutations].addValue(variantPhasing.Where(_ => _.nSomaticVariants() == nSomaticMutations).Count()));
            state.distributionOfPhasedPairs.addValue(variantPhasing.Count());

            state.totalPhasedVariantPairs += variantPhasing.Count();

            ASETools.chromosomes.ToList().ForEach(chr => state.countOfPhasedPairsByChromosome[chr] += variantPhasing.Where(_ => _.contig == chr).Count());

            foreach (var variant in variantPhasing)
            {
                foreach (var gene in commonData.geneMap.getGenesMappedTo(variant.contig, variant.locus0))
                {
                    if (!state.countOfPhasedPairsByHugoSymbol.ContainsKey(gene.hugoSymbol))
                    {
                        state.countOfPhasedPairsByHugoSymbol.Add(gene.hugoSymbol, 0);
                        state.countOfInconsistentPhasedPairsByHugoSymbol.Add(gene.hugoSymbol, 0);
                        state.countOfPhasedPairsWhereDNAandRNAdisagreeByHugoSymbol.Add(gene.hugoSymbol, 0);
                    }

                    state.countOfPhasedPairsByHugoSymbol[gene.hugoSymbol]++;
                    if (!(variant.dnaConsistent && variant.rnaConsistent))
                    {
                        state.countOfInconsistentPhasedPairsByHugoSymbol[gene.hugoSymbol]++;
                    }

                    if (variant.dnaConsistent && variant.rnaConsistent && !variant.mutuallyConsistent || gene.hugoSymbol == "COL6A3" && variant.nSomaticVariants() == 1)
                    {
                        state.countOfPhasedPairsWhereDNAandRNAdisagreeByHugoSymbol[gene.hugoSymbol]++;
                        state.genesAndCasesWithDNAandRNAdisagreeing.Add(gene.hugoSymbol + "\t" + variant.nSomaticVariants() + "\t" + variant.contig + "\t" + variant.locus0 + "\t" + variant.locus1 + "\t" +
                            ASETools.ConvertToExcelString(variant.counts[true][true][true] + @"/" + variant.counts[true][true][false]) + "\t" +
                            ASETools.ConvertToExcelString(variant.counts[true][false][true] + @"/" + variant.counts[true][false][false]) + "\t" +
                            ASETools.ConvertToExcelString(variant.counts[false][true][true] + @"/" + variant.counts[false][true][false]) + "\t" +
                            ASETools.ConvertToExcelString(variant.counts[false][false][true] + @"/" + variant.counts[false][false][false]) + "\t" +
                            case_.case_id + "\t" + case_.normal_dna_filename + "\t" + case_.tumor_dna_filename + "\t" + case_.tumor_rna_filename);
                    }
                } // gene covered by this variant
            } // each variant

        } // HandleOneCase

        static void FinishUp(WorkerThreadState state)
        {
            lock (globalState)
            {
                globalState.merge(state);
            } // lock
        } // FinishUp
    }
}
