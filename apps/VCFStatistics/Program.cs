using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.Diagnostics;

namespace VCFStatistics
{
    class Program
    {
        static Dictionary<bool, Dictionary<string, ASETools.PreBucketedHistogram>> resultsByChromosome = new Dictionary<bool, Dictionary<string, ASETools.PreBucketedHistogram>>(); // selected -> chromosome (or "other") -> histogram
        static Dictionary<bool, ASETools.PreBucketedHistogram> resultsByCase = new Dictionary<bool, ASETools.PreBucketedHistogram>();   // selected -> histogram

        static Dictionary<bool, Dictionary<string, ASETools.RunningMeanAndStdDev>> statsByChromosome = new Dictionary<bool, Dictionary<string, ASETools.RunningMeanAndStdDev>>(); // selected -> chromosome (or "other") -> stats
        static Dictionary<bool, ASETools.RunningMeanAndStdDev> statsByCase = new Dictionary<bool, ASETools.RunningMeanAndStdDev>(); // selected-> stats


        static List<string> chromosomesAndOther = new List<string>();

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            var configuration = ASETools.Configuration.loadFromFile(args);
            if (null == configuration)
            {
                Console.WriteLine("Unable to load configuration.");
                return;
            }

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname).Select(_ => _.Value).ToList();
            if (null == cases)
            {
                Console.WriteLine("Unable to load cases.");
                return;
            }

            if (cases.Any(_ => _.vcf_statistics_filename == "" || _.annotated_selected_variants_filename == ""))
            {
                Console.WriteLine("Some cases are missing VCF statistics and/or annotated selected variants.");
                return;
            }

            const int maxBucketForChromosome = 500000;
            const int incrementForChromosome = 2500;

            const int maxBucketForCase = 7000000;
            const int incrementForCase = 25000;

            ASETools.chromosomes.ToList().ForEach(chr => chromosomesAndOther.Add(chr));
            chromosomesAndOther.Add("other");

            ASETools.BothBools.ToList().ForEach(selected => resultsByChromosome.Add(selected, new Dictionary<string, ASETools.PreBucketedHistogram>()));
            ASETools.BothBools.ToList().ForEach(selected => chromosomesAndOther.ForEach(chr => resultsByChromosome[selected].Add(chr, new ASETools.PreBucketedHistogram(0, maxBucketForChromosome, incrementForChromosome))));
            ASETools.BothBools.ToList().ForEach(selected => resultsByCase.Add(selected, new ASETools.PreBucketedHistogram(0, maxBucketForCase, incrementForCase)));

            ASETools.BothBools.ToList().ForEach(selected => statsByChromosome.Add(selected, new Dictionary<string, ASETools.RunningMeanAndStdDev>()));
            ASETools.BothBools.ToList().ForEach(selected => chromosomesAndOther.ForEach(chr => statsByChromosome[selected].Add(chr, new ASETools.RunningMeanAndStdDev())));
            ASETools.BothBools.ToList().ForEach(selected => statsByCase.Add(selected, new ASETools.RunningMeanAndStdDev()));
 
            int nPerDot;
            var casesToRun = cases.Where(_ => _.vcf_filename != "").ToList();

            ASETools.PrintMessageAndNumberBar("Processing", "cases", casesToRun.Count(), out nPerDot);
            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(casesToRun, HandleOneCase, null, null, nPerDot);
            threading.run();

            string outputFilename = configuration.finalResultsDirectory + ASETools.VCFStatisticsFilename; 
            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
            if (null == outputFile)
            {
                Console.WriteLine("Unable to open output file " + outputFilename);
                return;
            }

            foreach (var selected in ASETools.BothBools) {

                outputFile.WriteLine("Distribution of germline variants (of any quality) by case" + (selected ? " for annotated selected variants" : ""));
                resultsByCase[selected].WriteHistogram(outputFile);

                outputFile.WriteLine();
                outputFile.WriteLine("CDFs of count of germline variants (of any quality) by case by chromosome" + (selected ? " for annotated selected variants" : ""));
                ASETools.PreBucketedHistogram.WriteBatchOfHistogramCDFs(outputFile, resultsByChromosome[selected].ToList());

                outputFile.WriteLine();
                outputFile.WriteLine("PDFs of count of germline variants (of any quality) by case by chromosome" + (selected ? " for annotated selected variants" : ""));
                ASETools.PreBucketedHistogram.WriteBatchOfHistogramPDFs(outputFile, resultsByChromosome[selected].ToList());
            } // selected

            outputFile.WriteLine("**done**");
            outputFile.Close();

            //
            // Now do the same thing with selected variants.
            //
            casesToRun = cases.Where(_ => _.tentative_selected_variants_filename != "").ToList();


            Console.WriteLine(ASETools.ElapsedTimeInSeconds(timer));
            
        } // Main

        static string chrToChrOrOther(string chr)
        {
            if (ASETools.chromosomes.Contains(chr))
            {
                return chr;
            }

            return "other";
        }

        static void HandleOneCase(ASETools.Case case_, int state)
        {

            var vcfStatistics = ASETools.VCFStatistics.LoadFromFile(case_.vcf_statistics_filename);
            if (null == vcfStatistics)
            {
                throw new Exception("Unable to load statistics for case " + case_.case_id + " from " + case_.vcf_statistics_filename);
            }

            if (vcfStatistics.vcfCounts["chr1"] < 5000)
            {
                Console.WriteLine("Warning: case " + case_.case_id + " has a VCF file that has fewer than 5000 variants in chr1 in VCF file " + case_.vcf_filename);
            }

            lock (resultsByChromosome[false])
            {
                foreach (var chr in vcfStatistics.vcfCounts.Select(_ => _.Key))
                {
                    resultsByChromosome[false][chr].addValue(vcfStatistics.vcfCounts[chr]);
                    statsByChromosome[false][chr].addValue(vcfStatistics.vcfCounts[chr]);
                }

                resultsByCase[false].addValue(vcfStatistics.totalVariants());
                statsByCase[false].addValue(vcfStatistics.totalVariants());
            } // lock

            if (case_.annotated_selected_variants_filename != "")
            {
                var countsByChromosome = ASETools.AnnotatedVariant.readFile(case_.annotated_selected_variants_filename).GroupByToDict(_ => chrToChrOrOther(_.contig));

                lock (resultsByChromosome[true])
                {
                    int total = 0;
                    foreach (var chr in chromosomesAndOther)
                    {
                        if (!countsByChromosome.ContainsKey(chr))
                        {
                            resultsByChromosome[true][chr].addValue(0);
                            statsByChromosome[true][chr].addValue(0);
                        } else
                        {
                            int n = countsByChromosome[chr].Count();
                            resultsByChromosome[true][chr].addValue(n);
                            statsByChromosome[true][chr].addValue(n);

                            total += n;
                        }
                    } // chr


                    resultsByCase[true].addValue(total);
                    statsByCase[true].addValue(total);
                }// lock
            }

        } // HandleOneCase
    }
}
