using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.Diagnostics;
using System.IO;

namespace UniparentalDisomy
{
    class Program
    {
        static ASETools.CommonData commonData;

        static StreamWriter outputFile;

        static ASETools.PreBucketedHistogram globalTumorHistogram = new ASETools.PreBucketedHistogram(0, 1, 0.01);
        static Dictionary<string, ASETools.PreBucketedHistogram> globalTumorPerChromosomeHistograms = new Dictionary<string, ASETools.PreBucketedHistogram>();
        static Dictionary<string, ASETools.PreBucketedHistogram> globalPerDiseaseHistograms = new Dictionary<string, ASETools.PreBucketedHistogram>();
        static Dictionary<string, ASETools.PreBucketedHistogram> globalPerGeneHistograms = new Dictionary<string, ASETools.PreBucketedHistogram>();

        const int minReads = 10;

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);

            if (null == commonData) return;

            if (commonData.listOfCases.Any(_ => _.tentative_annotated_selected_variants_filename == "" || _.tumor_copy_number_file_id != "" && _.tumor_copy_number_filename == ""))     // We use tentative, because the final selection process tosses cases where the tumor DNA is out of balance, which is the signal here
            {
                Console.WriteLine("Not all cases have input files.  Generate/download them and rerun this.");
                return;
            }

            foreach (var chromosome in ASETools.chromosomes)
            {
                globalTumorPerChromosomeHistograms.Add(chromosome, new ASETools.PreBucketedHistogram(0, 1, 0.01));
            }

            foreach (var disease in commonData.diseases)
            {
                globalPerDiseaseHistograms.Add(disease, new ASETools.PreBucketedHistogram(0, 1, 0.01));
            }

            outputFile = ASETools.CreateStreamWriterWithRetry(commonData.configuration.finalResultsDirectory + ASETools.UniparentalDisomyFilename);
            if (null == outputFile) return;

            var histogramsFile = ASETools.CreateStreamWriterWithRetry(commonData.configuration.finalResultsDirectory + ASETools.UniparentalDisomyHistogramsFilename);
            if (null == histogramsFile)
            {
                return;
            }

                outputFile.WriteLine("Case ID\tdisease\tchromosome\tn Germline Variant Sites\tn Tested Normal\tn Tested Tumor\tfrac 0.9 one allele normal\tfrac 0.9 one allele tumor");

            var casesToUse = commonData.listOfCases.Where(_ => _.tumor_copy_number_filename != "").ToList();

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "cases", casesToUse.Count(), out nPerDot);
            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(casesToUse, HandleOneCase, null, null, nPerDot);
            threading.run();

            outputFile.WriteLine("**done**");
            outputFile.Close();

            histogramsFile.WriteLine("Overall distribution");
            histogramsFile.WriteLine(ASETools.HistogramResultLine.Header());
            globalTumorHistogram.ComputeHistogram().ToList().ForEach(_ => histogramsFile.WriteLine(_));

            foreach (var chromosome in ASETools.chromosomes)
            {
                histogramsFile.WriteLine();
                histogramsFile.WriteLine(chromosome + " (n = " + globalTumorPerChromosomeHistograms[chromosome].count() + ")");
                histogramsFile.WriteLine(ASETools.HistogramResultLine.Header());
                globalTumorPerChromosomeHistograms[chromosome].ComputeHistogram().ToList().ForEach(_ => histogramsFile.WriteLine(_));
            }

            foreach (var disease in commonData.diseases)
            {
                histogramsFile.WriteLine();
                histogramsFile.WriteLine(disease + " (n = " + globalPerDiseaseHistograms[disease].count() + ")");
                histogramsFile.WriteLine(ASETools.HistogramResultLine.Header());
                globalPerDiseaseHistograms[disease].ComputeHistogram().ToList().ForEach(_ => histogramsFile.WriteLine(_));
            }

            histogramsFile.WriteLine("**done**");
            histogramsFile.Close();

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(commonData.timer));
        } // Main

        static void HandleOneCase(ASETools.Case case_, int state)
        {
            var annotatedSelectedVariants = ASETools.AnnotatedVariant.readFile(case_.tentative_annotated_selected_variants_filename);
            if (annotatedSelectedVariants == null)
            {
                throw new Exception("Unable to read tentative annotated selected variants for case " + case_.case_id + " from file " + case_.annotated_selected_variants_filename);
            }

            var tumorCopyNumber = ASETools.CopyNumberVariation.ReadFile(case_.tumor_copy_number_filename);
            if (null == tumorCopyNumber)
            {
                throw new Exception("Unable to read tumor copy number file for case " + case_.case_id + " from file " + case_.tumor_copy_number_filename);
            }

            var allChromosomesHistogram = new ASETools.PreBucketedHistogram(0, 1, 0.01);

            var perGeneHistograms = new Dictionary<string, ASETools.PreBucketedHistogram>();
            var outputLines = new List<string>();

            foreach (var chromosome in ASETools.chromosomes)
            {
                var perChromosomeHistogram = new ASETools.PreBucketedHistogram(0, 1, 0.01);

                var variantsForThisChromosome = annotatedSelectedVariants.Where(_ => _.contig == chromosome && !_.somaticMutation).ToList();

                var nVariants = variantsForThisChromosome.Count();
                var normalUseful = variantsForThisChromosome.Where(_ => _.normalDNAReadCounts.usefulReads() >= minReads).ToList();
                var tumorUseful = variantsForThisChromosome.Where(_ => _.tumorDNAReadCounts.usefulReads() >= minReads && 
                                        tumorCopyNumber.Where(cnv => cnv.OverlapsLocus(ASETools.chromosomeNameToNonChrForm(chromosome), _.locus, _.locus + 1)).Where(r => Math.Abs(r.Segment_Mean) > 1.0).ToList().Count() > 0).ToList(); // Exclude CNVs

                if (tumorUseful.Count() < 10)
                {
                    continue;
                }

                var tumor90PercentFraction = (double)tumorUseful.Where(_ => _.tumorDNAReadCounts.AlleleSpecificValue() >= .8).Count() / tumorUseful.Count();

                perChromosomeHistogram.addValue(tumor90PercentFraction);
                allChromosomesHistogram.addValue(tumor90PercentFraction);

                outputLines.Add(case_.case_id + "\t" + case_.disease() + "\t" + chromosome + "\t" + nVariants + "\t" + normalUseful.Count() + "\t" + tumorUseful.Count() + "\t" +
                    ((double)normalUseful.Where(_ => _.normalDNAReadCounts.AlleleSpecificValue() >= .8).Count() / normalUseful.Count()) + "\t" + tumor90PercentFraction);

                lock (globalTumorPerChromosomeHistograms)
                {
                    globalTumorPerChromosomeHistograms[chromosome].merge(perChromosomeHistogram);                        
                }
            } // foreach chromosome

            lock (globalTumorHistogram)
            {
                globalTumorHistogram.merge(allChromosomesHistogram);
                globalPerDiseaseHistograms[case_.disease()].merge(allChromosomesHistogram);

                foreach (var per_gene_entry in perGeneHistograms)
                {
                    var hugo_symbol = per_gene_entry.Key;
                    var histogram = per_gene_entry.Value;

                    if (!globalPerDiseaseHistograms.ContainsKey(hugo_symbol))
                    {
                        globalPerGeneHistograms.Add(hugo_symbol, histogram);
                    } else
                    {
                        globalPerGeneHistograms[hugo_symbol].merge(histogram);
                    } 
                } // genes

                foreach (var outputLine in outputLines)
                {
                    outputFile.WriteLine(outputLine);
                }
            } // lock
        } // HandleOneCase
    }
}
