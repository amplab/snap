using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using ASELib;
using System.Threading;

namespace ASEDistributionBySite
{
    class Program
    {
        static ASETools.GeneMap geneMap;
        static ASETools.Configuration configuration;
        static Dictionary<string, ASETools.ASEMapPerGeneLine> perGeneASEMap;
        static ASETools.Histogram tumorResultsHistogram = new ASETools.Histogram();
        static ASETools.Histogram normalResultsHistogram = new ASETools.Histogram();
        static int nProcessed = 0;


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

            if (configuration.commandLineArgs.Count() != 1)
            {
                Console.WriteLine("usage: ASEDistributionBySite outputFilename");
                return;
            }

            var geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(configuration.geneLocationInformationFilename));
            geneMap = new ASETools.GeneMap(geneLocationInformation.genesByName);
            perGeneASEMap = ASETools.ASEMapPerGeneLine.ReadFromFileToDictionary(configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

            if (null == cases)
            {
                Console.WriteLine("Unable to load cases.");
            }

            var casesToProcess = cases.Where(x => x.Value.annotated_selected_variants_filename != "" && x.Value.tumor_copy_number_filename != "").Select(x => x.Value).ToList();
            int nCases = casesToProcess.Count();
            Console.Write("Processing " + nCases + " cases, 1 dot/100 cases: ");

            var threads = new List<Thread>();
            for (int i = 0; i < Environment.ProcessorCount; i++)
            {
                threads.Add(new Thread(() => WorkerThread(casesToProcess)));
            }

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            var tumorResults = tumorResultsHistogram.ComputeHistogram(0, 1, 0.01);
            var normalResults = normalResultsHistogram.ComputeHistogram(0, 1, 0.01);

            var outputFile = ASETools.CreateStreamWriterWithRetry(args[0]);


            outputFile.WriteLine("Min value\ttumor count\ttumor pdf\ttumor cdf\tnormal count\tnormal pdf\tnormal cdf");

            for (int i = 0; i < tumorResults.Count(); i++)
            {
                outputFile.WriteLine(tumorResults[i].minValue + "\t" + tumorResults[i].count + "\t" + tumorResults[i].pdfValue + "\t" + tumorResults[i].cdfValue + "\t" + normalResults[i].count + "\t" + normalResults[i].pdfValue + "\t" + normalResults[i].cdfValue);
            }

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine();
            Console.WriteLine("Processed " + nCases + " in " + ASETools.ElapsedTimeInSeconds(timer));
        } // Main

        static void WorkerThread(List<ASETools.Case> cases)
        {
            var localTumorHistogram = new ASETools.Histogram();
            var localNormalHistogram = new ASETools.Histogram();

            while (true)
            {
                ASETools.Case case_ = null;

                lock (cases)
                {
                    if (cases.Count() == 0)
                    {
                        tumorResultsHistogram.merge(localTumorHistogram);
                        normalResultsHistogram.merge(localNormalHistogram);
                        return;
                    }

                    nProcessed++;
                    if (nProcessed % 100 == 0)
                    {
                        Console.Write(".");
                    }

                    case_ = cases[0];
                    cases.RemoveAt(0);
                }

                var copyNumber = ASETools.CopyNumberVariation.ReadBothFiles(case_);
                var annotatedSelectedVariants = ASETools.AnnotatedVariant.readFile(case_.annotated_selected_variants_filename).Where(x => !x.somaticMutation).ToList();

                if (annotatedSelectedVariants == null)
                {
                    Console.WriteLine("Unable to read annotated selected variants from " + case_.annotated_selected_variants_filename);
                    continue;
                }

                foreach (var variant in annotatedSelectedVariants)
                {
                    if (variant.IsASECandidate(true, copyNumber, configuration, perGeneASEMap, geneMap))
                    {
                        localTumorHistogram.addValue(variant.GetTumorAlleleSpecificExpression());
                    }

                    if (variant.IsASECandidate(false, copyNumber, configuration, perGeneASEMap, geneMap))
                    {
                        localNormalHistogram.addValue(variant.GetNormalAlleleSpecificExpression());
                    }
                }

            } // while true
        }
    }
}
