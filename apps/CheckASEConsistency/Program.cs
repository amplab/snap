using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;
using System.IO;
using System.Diagnostics;
using ASELib;

namespace CheckASEConsistency // Looks for single <case,gene> pairs that have ASE difference of more than 0.5, which might indicate something funny going on.
{
    class Program
    {
        static ASETools.Configuration configuration;
        static List<ASETools.Case> casesToProcess = new List<ASETools.Case>();
        static ASETools.GeneMap geneMap;
        static StreamWriter outputFile;
        static int nCasesProcessed = 0;
        static int nFunnyInstances = 0;
        static int nFunnyWithASomaticMutation = 0;
        static int nInstancesChecked = 0;
        static int nAnnotatedSelectedVariants = 0;

        static Dictionary<string, int> geneToFunnyCount = new Dictionary<string, int>();

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            configuration = ASETools.Configuration.loadFromFile(args);
            if (null == configuration)
            {
                Console.WriteLine("Couldn't load configuration.");
                return;
            }

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

            casesToProcess = cases.Select(x => x.Value).Where(x => x.annotated_selected_variants_filename != null && x.annotated_selected_variants_filename != "").ToList();

            var geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(configuration.geneLocationInformationFilename));
            geneMap = new ASETools.GeneMap(geneLocationInformation.genesByName);

            var outputFilename = configuration.finalResultsDirectory + ASETools.ASEConsistencyFilename;
            outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
            if (null == outputFile)
            {
                Console.WriteLine("Couldn't open output file " + outputFilename);
                return;
            }

            outputFile.WriteLine("Case ID\tcontig\tlocus\tGenes\tn Variants\tminASE\tmaxASE\teach ASE\teach variant is somatic\teach tumor RNA ref count\teach tumor RNA alt count");

            Console.Write("Processing " + casesToProcess.Count() + " of " + cases.Count() + " cases; 1 dot/100 cases:");

            var threads = new List<Thread>();
            for (int i = 0; i < Environment.ProcessorCount; i++)
            {
                threads.Add(new Thread(() => WorkerThread()));
            }

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            Console.WriteLine();    // Newline after progress dots

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine("Gene\t# Instances with funny ASE");
            foreach (var geneEntry in geneToFunnyCount)
            {
                Console.WriteLine(geneEntry.Key + "\t" + geneEntry.Value);
            }

            Console.WriteLine();
            Console.WriteLine("Processed " + nInstancesChecked + " instances from among " + nAnnotatedSelectedVariants + 
                " annotated selected variants, of which " + nFunnyInstances + " were questionable (" + nFunnyWithASomaticMutation + 
                " of those had at least one somatic variant) in " + ASETools.ElapsedTimeInSeconds(timer));
        } // Main

        static void WorkerThread()
        {
            int nInstancesCheckedThisThread = 0;
            int nAnnotatedSelectedVariantsThisThread = 0;

            while (true)
            {
                ASETools.Case case_;
                lock (casesToProcess)
                {
                    if (casesToProcess.Count() == 0)
                    {
                        nInstancesChecked += nInstancesCheckedThisThread;
                        nAnnotatedSelectedVariants += nAnnotatedSelectedVariantsThisThread;
                        return;
                    }

                    case_ = casesToProcess[0];
                    casesToProcess.RemoveAt(0);
                } // lock (casesToProcess)

                var annotatedSelectedVariants = ASETools.AnnotatedVariant.readFile(case_.annotated_selected_variants_filename);
                nAnnotatedSelectedVariantsThisThread += annotatedSelectedVariants.Count();

                annotatedSelectedVariants.Sort(ASETools.AnnotatedVariant.CompareByLocus);   // Puts them in genome order

                if (annotatedSelectedVariants == null)
                {
                    Console.WriteLine("Unable to load annotated selected variants for case " + case_.case_id + " from file " + case_.annotated_selected_variants_filename);
                    continue;
                }

                List<ASETools.GeneLocationInfo> currentGenesMappedTo = new List<ASETools.GeneLocationInfo>();
                List<ASETools.AnnotatedVariant> variantsInThisGroup = new List<ASETools.AnnotatedVariant>();
                for (int i = 0; i < annotatedSelectedVariants.Count(); i++)
                {
                    var variant = annotatedSelectedVariants[i];
                    var genesForThisVariant = geneMap.getGenesMappedTo(variant.contig, variant.locus);
                    if (genesForThisVariant != currentGenesMappedTo)
                    {
                        ProcessSet(case_, currentGenesMappedTo, variantsInThisGroup, ref nInstancesCheckedThisThread);
                        currentGenesMappedTo = genesForThisVariant;
                        variantsInThisGroup = new List<ASETools.AnnotatedVariant>();
                        variantsInThisGroup.Add(variant);
                    }
                }

                ProcessSet(case_, currentGenesMappedTo, variantsInThisGroup, ref nInstancesCheckedThisThread);

                lock (casesToProcess)
                {
                    nCasesProcessed++;
                    if (nCasesProcessed % 100 == 0)
                    {
                        Console.Write(".");
                    }
                } // lock (casesToProcess)
            } // while (true)
        } // WorkerThread

        static void ProcessSet(ASETools.Case case_, List<ASETools.GeneLocationInfo> currentGenesMappedTo, 
            List<ASETools.AnnotatedVariant> variantsInThisGroup, ref int nInstancesCheckedThisThread)
        {
            var variantsToConsider = variantsInThisGroup.Where(x => x.IsASECandidate()).ToList();

            if (variantsToConsider.Count() > 1)
            {
                nInstancesCheckedThisThread++;

                var minASE = variantsToConsider.Select(x => x.tumorRNAReadCounts.AlleleSpecificValue()).Min();
                var maxASE = variantsToConsider.Select(x => x.tumorRNAReadCounts.AlleleSpecificValue()).Max();

                if (maxASE - minASE > 0.5)
                {
                    lock (outputFile)
                    {
                        nFunnyInstances++;
                        if (variantsToConsider.Select(x => x.somaticMutation).Count() > 0)
                        {
                            nFunnyWithASomaticMutation++;
                        }

                        outputFile.Write(case_.case_id + "\t" + variantsToConsider[0].contig + "\t" + variantsToConsider[0].locus + "\t");
                        currentGenesMappedTo.ForEach(x => outputFile.Write(x.hugoSymbol + ","));
                        outputFile.Write("\t" + variantsToConsider.Count() + "\t" + minASE + "\t" + maxASE + "\t");
                        variantsToConsider.ForEach(x => outputFile.Write(x.tumorRNAReadCounts.AlleleSpecificValue() + ","));
                        outputFile.Write("\t");
                        variantsToConsider.ForEach(x => outputFile.Write(x.somaticMutation + ","));
                        outputFile.Write("\t");
                        variantsToConsider.ForEach(x => outputFile.Write(x.tumorRNAReadCounts.nMatchingReference + ","));
                        outputFile.Write("\t");
                        variantsToConsider.ForEach(x => outputFile.Write(x.tumorRNAReadCounts.nMatchingAlt + ","));
                        outputFile.WriteLine();

                        foreach (var gene in currentGenesMappedTo)
                        {
                            if (!geneToFunnyCount.ContainsKey(gene.hugoSymbol))
                            {
                                geneToFunnyCount.Add(gene.hugoSymbol, 0);
                            }
                            geneToFunnyCount[gene.hugoSymbol]++;
                        }
                    } // lock (outputFile)
                } // if ASE differs by more than 0.5
            } // if we have at least two variants in this gene set
        } // ProcessSet


    } // Program
} // CheckASEConsistency
