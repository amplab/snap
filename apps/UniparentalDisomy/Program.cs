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
        static Dictionary<string, ASETools.PreBucketedHistogram> globalTumorHistogramByGender = new Dictionary<string, ASETools.PreBucketedHistogram>();
        static Dictionary<string, ASETools.PreBucketedHistogram> globalTumorPerChromosomeHistograms = new Dictionary<string, ASETools.PreBucketedHistogram>();
        static Dictionary<string, ASETools.PreBucketedHistogram> globalPerDiseaseHistograms = new Dictionary<string, ASETools.PreBucketedHistogram>();
        static Dictionary<string, ASETools.PreBucketedHistogram> globalPerGeneHistograms = new Dictionary<string, ASETools.PreBucketedHistogram>();
        //
        // The fraction histograms don't use absolute bases, instead they use fraction of the chromosome.  So, they run from 0 to 1.
        //
        static Dictionary<string, ASETools.PreBucketedHistogram> globalBeginningOfChromosomeFractionByChromosome = new Dictionary<string, ASETools.PreBucketedHistogram>();
        static Dictionary<string, ASETools.PreBucketedHistogram> globalEndOfChromosomeFractionByChromosome = new Dictionary<string, ASETools.PreBucketedHistogram>();
        static Dictionary<string, ASETools.PreBucketedHistogram> globalBeginningOfChromosomeFractionByDisease = new Dictionary<string, ASETools.PreBucketedHistogram>();
        static Dictionary<string, ASETools.PreBucketedHistogram> globalEndOfChromosomeFractionByDisease = new Dictionary<string, ASETools.PreBucketedHistogram>();

        //
        // The distance histograms use absolute distance, and they're all scaled to the largest chromosome.
        //
        static Dictionary<string, ASETools.PreBucketedHistogram> globalBeginningOfChromosomeDistanceByChromosome = new Dictionary<string, ASETools.PreBucketedHistogram>();
        static Dictionary<string, ASETools.PreBucketedHistogram> globalEndOfChromosomeDistanceByChromosome = new Dictionary<string, ASETools.PreBucketedHistogram>();
        static Dictionary<string, ASETools.PreBucketedHistogram> globalBeginningOfChromosomeDistanceByDisease = new Dictionary<string, ASETools.PreBucketedHistogram>();
        static Dictionary<string, ASETools.PreBucketedHistogram> globalEndOfChromosomeDistanceByDisease = new Dictionary<string, ASETools.PreBucketedHistogram>();

        //
        // Soft distance is farthest point that's >= 0.9 and has a mean over all points from the beginning of the chromosome to there of >= 0.9
        //
        static Dictionary<string, ASETools.PreBucketedHistogram> globalSoftDistanceFromBeginningOfChromosome = new Dictionary<string, ASETools.PreBucketedHistogram>();
        static Dictionary<string, ASETools.PreBucketedHistogram> globalSoftDistanceFromEndOfChromosome = new Dictionary<string, ASETools.PreBucketedHistogram>();

        //
        // Distance to first non-homozygous measures the first one that's not, so it overshoots in the same sense that the ordinary distance undershoots.
        //
        static Dictionary<string, ASETools.PreBucketedHistogram> globalBeginningOfChromosomeDistanceToFirstNonHomozygousByChromosome = new Dictionary<string, ASETools.PreBucketedHistogram>();
        static Dictionary<string, ASETools.PreBucketedHistogram> globalEndOfChromosomeDistanceToFirstNonHomozygousByChromosome = new Dictionary<string, ASETools.PreBucketedHistogram>();

        //
        // Rather than thresholding at 90% one allele (0.8 ASV), we look at mean ASV for each chromosome.
        //
        static ASETools.PreBucketedHistogram globalASVHistogram = new ASETools.PreBucketedHistogram(0, 1, 0.01);
        static Dictionary<string, ASETools.PreBucketedHistogram> globalMeanASVByChromosomeHistogram = new Dictionary<string, ASETools.PreBucketedHistogram>();
        static Dictionary<string, ASETools.PreBucketedHistogram> globalBeginningOfChromosomeDistanceToFirstNonHomozygousByChromosomeSoft = new Dictionary<string, ASETools.PreBucketedHistogram>();
        static Dictionary<string, ASETools.PreBucketedHistogram> globalEndOfChromosomeDistanceToFirstNonHomozygousByChromosomeSoft = new Dictionary<string, ASETools.PreBucketedHistogram>();

        //
        // Splitting into with and without a particular gene.
        //
        static Dictionary<string, Dictionary<bool, ASETools.PreBucketedHistogram>> globalWithAndWithoutGeneDistanceSoftBeginningToLastHomozygous = new Dictionary<string, Dictionary<bool, ASETools.PreBucketedHistogram>>(); // Gene -> mutated -> Histogram  (only the chromosome containing the gene)
        static Dictionary<string, Dictionary<bool, ASETools.PreBucketedHistogram>> globalWithAndWithoutGeneDistanceSoftEndToLastHomozygous = new Dictionary<string, Dictionary<bool, ASETools.PreBucketedHistogram>>(); // Gene -> mutated -> Histogram  (only the chromosome containing the gene)
        static Dictionary<string, Dictionary<bool, ASETools.PreBucketedHistogram>> globalWithAndWithoutGeneDistanceSoftBeginningToFirstNonhomozygous = new Dictionary<string, Dictionary<bool, ASETools.PreBucketedHistogram>>(); // Gene -> mutated -> Histogram  (only the chromosome containing the gene)
        static Dictionary<string, Dictionary<bool, ASETools.PreBucketedHistogram>> globalWithAndWithoutGeneDistanceSoftEndToFirstNonhomozygous = new Dictionary<string, Dictionary<bool, ASETools.PreBucketedHistogram>>(); // Gene -> mutated -> Histogram  (only the chromosome containing the gene)

        static Dictionary<int, ASETools.RunningMeanAndStdDev> globalExpressionOfMostCommonAlleleByAgeAtDiagnosis = new Dictionary<int, ASETools.RunningMeanAndStdDev>();

        static Dictionary<string, Dictionary<int, Dictionary<int, int>>> globalGeneByFlankingHomozygousSites = new Dictionary<string, Dictionary<int, Dictionary<int, int>>>();     // Maps hugo symbol -> 0,1,2 (mutation count) -> {0,1,2} flanking homozygous sites -> count

        static Dictionary<string, ASETools.PreBucketedHistogram> globalChr6ByDisease = new Dictionary<string, ASETools.PreBucketedHistogram>();


        static string[] genesToConsider = { "TP53", "VHL", "SAMD9", "KEAP1", "KRAS", "CDKN2A", "SAMD9L", "ATM", "SMAD4", "STK11", "PTEN", "CDH1", "PBRM1", "MUC16", "TTN", "ACTB", "HLA-A", "HLA-B", "HLA-C"};


        static List<string> chromosomesToUse = new List<string>();


        static readonly double largestChromosomeInMegabases = (double)ASETools.chromosomeSizes.Select(_ => _.size).Max()/1000000;

        const int minReads = 10;

        delegate ASETools.ReadCounts ReadCountGetter(ASETools.AnnotatedVariant variant);

        static ReadCountGetter GetReadCount = _ => _.tumorDNAReadCounts;

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);

            if (null == commonData) return;

            var minAgeAtDiagnosis = commonData.clinicalSummariesByPatientId.Where(_ => _.Value.age_at_diagnosis != -1).Select(_ => _.Value.age_at_diagnosis).Min();
            var maxAgeAtDiagnosis = commonData.clinicalSummariesByPatientId.Where(_ => _.Value.age_at_diagnosis != -1).Select(_ => _.Value.age_at_diagnosis).Max();

            for (int i = minAgeAtDiagnosis; i <= maxAgeAtDiagnosis; i++)
            {
                globalExpressionOfMostCommonAlleleByAgeAtDiagnosis.Add(i, new ASETools.RunningMeanAndStdDev());
            }

            ASETools.autosomes.ToList().ForEach(_ => chromosomesToUse.Add(_));
            chromosomesToUse.Add("chrx");

            if (commonData.listOfCases.Any(_ => _.tentative_annotated_selected_variants_filename == "" || _.tumor_copy_number_file_id != "" && _.tumor_copy_number_filename == ""))     // We use tentative, because the final selection process tosses cases where the tumor DNA is out of balance, which is the signal here
            {
                Console.WriteLine("Not all cases have input files.  Generate/download them and rerun this.");
                return;
            }

            foreach (var gender in ASETools.BothGenders)
            {
                globalTumorHistogramByGender.Add(gender, new ASETools.PreBucketedHistogram(0, 1, 0.01));
            }            

            foreach (var geneToConsider in genesToConsider)
            {
                globalWithAndWithoutGeneDistanceSoftBeginningToLastHomozygous.Add(geneToConsider, new Dictionary<bool, ASETools.PreBucketedHistogram>());
                globalWithAndWithoutGeneDistanceSoftEndToLastHomozygous.Add(geneToConsider, new Dictionary<bool, ASETools.PreBucketedHistogram>());
                globalWithAndWithoutGeneDistanceSoftBeginningToFirstNonhomozygous.Add(geneToConsider, new Dictionary<bool, ASETools.PreBucketedHistogram>());
                globalWithAndWithoutGeneDistanceSoftEndToFirstNonhomozygous.Add(geneToConsider, new Dictionary<bool, ASETools.PreBucketedHistogram>());

                foreach (var mutated in ASETools.BothBools)
                {
                    globalWithAndWithoutGeneDistanceSoftBeginningToLastHomozygous[geneToConsider].Add(mutated, new ASETools.PreBucketedHistogram(0, ASETools.chromosomeSizesByName[commonData.geneLocationInformation.genesByName[geneToConsider].chromosome].size / ASETools.Million, 1));
                    globalWithAndWithoutGeneDistanceSoftEndToLastHomozygous[geneToConsider].Add(mutated, new ASETools.PreBucketedHistogram(0, ASETools.chromosomeSizesByName[commonData.geneLocationInformation.genesByName[geneToConsider].chromosome].size / ASETools.Million, 1));
                    globalWithAndWithoutGeneDistanceSoftBeginningToFirstNonhomozygous[geneToConsider].Add(mutated, new ASETools.PreBucketedHistogram(0, ASETools.chromosomeSizesByName[commonData.geneLocationInformation.genesByName[geneToConsider].chromosome].size / ASETools.Million, 1));
                    globalWithAndWithoutGeneDistanceSoftEndToFirstNonhomozygous[geneToConsider].Add(mutated, new ASETools.PreBucketedHistogram(0, ASETools.chromosomeSizesByName[commonData.geneLocationInformation.genesByName[geneToConsider].chromosome].size / ASETools.Million, 1));
                }

                globalGeneByFlankingHomozygousSites.Add(geneToConsider, new Dictionary<int, Dictionary<int, int>>());
                foreach (var mutationCount in ASETools.ZeroOneTwo)
                {
                    globalGeneByFlankingHomozygousSites[geneToConsider].Add(mutationCount, new Dictionary<int, int>());
                    for (int i = 0; i <= 2; i++)
                    {
                        globalGeneByFlankingHomozygousSites[geneToConsider][mutationCount].Add(i, 0);
                    }
                }
            }

            foreach (var chromosome in chromosomesToUse)
            {
                globalTumorPerChromosomeHistograms.Add(chromosome, new ASETools.PreBucketedHistogram(0, 1, 0.01));
                globalBeginningOfChromosomeFractionByChromosome.Add(chromosome, new ASETools.PreBucketedHistogram(0, 1, 0.01));
                globalEndOfChromosomeFractionByChromosome.Add(chromosome, new ASETools.PreBucketedHistogram(0, 1, 0.01));
                globalBeginningOfChromosomeDistanceByChromosome.Add(chromosome, new ASETools.PreBucketedHistogram(0, largestChromosomeInMegabases, 1));
                globalEndOfChromosomeDistanceByChromosome.Add(chromosome, new ASETools.PreBucketedHistogram(0, largestChromosomeInMegabases, 1));
                globalBeginningOfChromosomeDistanceToFirstNonHomozygousByChromosome.Add(chromosome, new ASETools.PreBucketedHistogram(0, largestChromosomeInMegabases, 1));
                globalEndOfChromosomeDistanceToFirstNonHomozygousByChromosome.Add(chromosome, new ASETools.PreBucketedHistogram(0, largestChromosomeInMegabases, 1));
                globalMeanASVByChromosomeHistogram.Add(chromosome, new ASETools.PreBucketedHistogram(0, 1, 0.01));
                globalSoftDistanceFromBeginningOfChromosome.Add(chromosome, new ASETools.PreBucketedHistogram(0, largestChromosomeInMegabases, 1));
                globalSoftDistanceFromEndOfChromosome.Add(chromosome, new ASETools.PreBucketedHistogram(0, largestChromosomeInMegabases, 1));
                globalBeginningOfChromosomeDistanceToFirstNonHomozygousByChromosomeSoft.Add(chromosome, new ASETools.PreBucketedHistogram(0, largestChromosomeInMegabases, 1));
                globalEndOfChromosomeDistanceToFirstNonHomozygousByChromosomeSoft.Add(chromosome, new ASETools.PreBucketedHistogram(0, largestChromosomeInMegabases, 1));
            }

            foreach (var disease in commonData.diseases)
            {
                globalPerDiseaseHistograms.Add(disease, new ASETools.PreBucketedHistogram(0, 1, 0.01));
                globalBeginningOfChromosomeFractionByDisease.Add(disease, new ASETools.PreBucketedHistogram(0, 1, 0.01));
                globalEndOfChromosomeFractionByDisease.Add(disease, new ASETools.PreBucketedHistogram(0, 1, 0.01));
                globalBeginningOfChromosomeDistanceByDisease.Add(disease, new ASETools.PreBucketedHistogram(0, largestChromosomeInMegabases, 1));
                globalEndOfChromosomeDistanceByDisease.Add(disease, new ASETools.PreBucketedHistogram(0, largestChromosomeInMegabases, 1));
                globalChr6ByDisease.Add(disease, new ASETools.PreBucketedHistogram(0, 1, 0.01));
            }

            outputFile = ASETools.CreateStreamWriterWithRetry(commonData.configuration.finalResultsDirectory + ASETools.UniparentalDisomyFilename);
            if (null == outputFile) return;

            var histogramsFile = ASETools.CreateStreamWriterWithRetry(commonData.configuration.finalResultsDirectory + ASETools.UniparentalDisomyHistogramsFilename);
            if (null == histogramsFile)
            {
                return;
            }

            outputFile.Write("Case ID\tgender\tvital status\tage at diagnosis\tdays to death\tdisease\tchromosome\tn Germline Variant Sites\tn Tested Normal\tn Tested Tumor\tfrac 0.9 one allele normal\tfrac 0.9 one allele tumor" +
                "\tfraction of bases consecutive from beginning with all at least 0.9 one allele\tfraction of bases consecutive from end with all at least 0.9 one allele\tn Consecutive tumor beginning\tn Consecutive tumor end\tMean ASV tumor" +
                "\tbeginning fraction soft consecutive\tend fraction soft consecutive");
            foreach (var gene in genesToConsider)
            {
                outputFile.Write("\t" + gene + " mutation count");
            }
            outputFile.WriteLine("\tmap\tspatial map");


            var casesToUse = commonData.listOfCases.Where(_ => _.tentative_asv_without_cnvs_filename != "").ToList();

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "cases", casesToUse.Count(), out nPerDot);
            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(casesToUse, HandleOneCase, null, null, nPerDot);
            threading.run();

            outputFile.WriteLine("**done**");
            outputFile.Close();

            histogramsFile.WriteLine("Genes by mutation count and fraction with neither flanking site homozygous");
            histogramsFile.WriteLine("hugo symbol\t0 mutations fraction with no homozygous flanking sites\t1 mutations fraction with no homozygous flanking sites\t>1 mutations fraction with no homozygous flanking sites\tp 0 mutations differs from 1 mutatio by Welch's T Test (two sided), Bonferroni corrected");
            foreach (var geneToConsider in genesToConsider)
            {
                var thisGene = globalGeneByFlankingHomozygousSites[geneToConsider];

                histogramsFile.Write(geneToConsider);
                foreach (var mutationCount in ASETools.ZeroOneTwo)
                {
                    histogramsFile.Write("\t" + (double)thisGene[mutationCount][0] / thisGene[mutationCount].Select(_ => _.Value).Sum());
                }


                List<double>[] distributions = new List<double>[2]; // index is mutation count
                for (int mutationCount = 0; mutationCount < 2; mutationCount++)
                {
                    distributions[mutationCount] = new List<double>();
                    for (int flankingMutationCount = 0; flankingMutationCount <= 2; flankingMutationCount++)
                    {
                        for (int i = 0; i < thisGene[mutationCount][flankingMutationCount]; i++)
                        {
                            if (flankingMutationCount == 0)
                            {
                                distributions[mutationCount].Add(0);    // There MUST be a better way to do this!
                            }
                            else
                            {
                                distributions[mutationCount].Add(1);
                            }
                        }
                    }
                }
                histogramsFile.WriteLine("\t" + ASETools.WelchsTTest.TwoSidedTTest(distributions[0], distributions[1]) * genesToConsider.Count()); // Multiplying by the count of genes is the Bonferroni correction
            }

            histogramsFile.WriteLine();

            histogramsFile.WriteLine("Genes by flanking homozygous count");
            histogramsFile.Write("Hugo Symbol");
            foreach (var mutationCount in ASETools.ZeroOneTwo)
            {
                for (int nFlankingHeterzygous = 0; nFlankingHeterzygous <= 2; nFlankingHeterzygous++)
                {
                    histogramsFile.Write("\t" + ASETools.ZeroOneManyString(mutationCount) + " mutations n with " + nFlankingHeterzygous + " flanking heterozygous");
                }

                for (int nFlankingHeterzygous = 0; nFlankingHeterzygous <= 2; nFlankingHeterzygous++)
                {
                    histogramsFile.Write("\t" + ASETools.ZeroOneManyString(mutationCount) + " mutations % of this mutation count with " + nFlankingHeterzygous + " flanking heterozygous");
                }
            }

            histogramsFile.WriteLine();

            foreach (var geneToConsider in genesToConsider)
            {
                var thisGene = globalGeneByFlankingHomozygousSites[geneToConsider];

                histogramsFile.Write(geneToConsider);

                foreach (var nMutations in ASETools.ZeroOneTwo)
                {
                    foreach (var nFlankingHeterozygous in ASETools.ZeroOneTwo)
                    {
                        histogramsFile.Write("\t" + thisGene[nMutations][nFlankingHeterozygous]);
                    }

                    var nInClass = thisGene[nMutations].Select(_ => _.Value).Sum();

                    foreach (var nFlankingHeterozygous in ASETools.ZeroOneTwo)
                    {
                        histogramsFile.Write("\t" + (double)thisGene[nMutations][nFlankingHeterozygous] / nInClass);
                    }
                }

                histogramsFile.WriteLine(); 
            } // geneToConsider

            histogramsFile.WriteLine();
            histogramsFile.WriteLine("Summary Of CDFs ");
            var headersAndHistograms = new List<KeyValuePair<string, ASETools.PreBucketedHistogram>>();
            headersAndHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("overall", globalTumorHistogram));

            chromosomesToUse.ToList().ForEach(_ => headersAndHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>(_, globalTumorPerChromosomeHistograms[_])));

            ASETools.PreBucketedHistogram.WriteBatchOfHistogramCDFs(histogramsFile, headersAndHistograms);

            writePerDiseaseHistograms("Summary by disease", histogramsFile, globalPerDiseaseHistograms);
            //writePerDiseaseHistograms("Distance at beginning of chromosome by disease", histogramsFile, globalBeginningOfChromosomeDistanceByDisease);
            //writePerDiseaseHistograms("Distance at end of chromosome by disease", histogramsFile, globalEndOfChromosomeDistanceByDisease);

            //writePerChromsomeHistograms("Distance at beginning of chromosome by chromosome", histogramsFile, globalBeginningOfChromosomeDistanceByChromosome);
            //writePerChromsomeHistograms("Distance from end of chromosome by chromosome", histogramsFile, globalEndOfChromosomeDistanceByChromosome);

            //writePerChromsomeHistograms("Distance at beginning of chromosome to first measured non-homozygous locus", histogramsFile, globalBeginningOfChromosomeDistanceToFirstNonHomozygousByChromosome);
            //writePerChromsomeHistograms("Distance from end of chromosome to first measured non-homozygous locus", histogramsFile, globalEndOfChromosomeDistanceToFirstNonHomozygousByChromosome);

            //writePerDiseaseHistograms("Fraction of chromosome at beginning of chromosome by disease", histogramsFile, globalBeginningOfChromosomeFractionByDisease);
            //writePerDiseaseHistograms("Fraction of chromosome at end of chromosome by disease", histogramsFile, globalEndOfChromosomeFractionByDisease);

            //writePerChromsomeHistograms("Fraction of chromosome at beginning of chromosome by chromosome", histogramsFile, globalBeginningOfChromosomeFractionByChromosome);
            //writePerChromsomeHistograms("Fraction of chromosome at end of chromosome by chromosome", histogramsFile, globalEndOfChromosomeFractionByChromosome);

            headersAndHistograms = new List<KeyValuePair<string, ASETools.PreBucketedHistogram>>();
            ASETools.BothGenders.ToList().ForEach(_ => headersAndHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>(_, globalTumorHistogramByGender[_])));
            ASETools.PreBucketedHistogram.WriteBatchOfHistogramCDFs(histogramsFile, headersAndHistograms);

            writePerChromsomeHistograms("Distance at beginning of chromosome (soft) by chromosome", histogramsFile, globalSoftDistanceFromBeginningOfChromosome);
            writePerChromsomeHistograms("Distance at beginning of chromosome to first non-homozygous locus (soft) by chromosome", histogramsFile, globalBeginningOfChromosomeDistanceToFirstNonHomozygousByChromosomeSoft);

            writePerChromsomeHistograms("Distance of chromosome at end of chromosome (soft) by chromosome", histogramsFile, globalSoftDistanceFromEndOfChromosome);
            writePerChromsomeHistograms("Distance at end of chromosome to first non-homozygous locus (soft) by chromosome", histogramsFile, globalEndOfChromosomeDistanceToFirstNonHomozygousByChromosomeSoft);

            writePerChromsomeHistograms("Mean ASV by chromosome", histogramsFile, globalMeanASVByChromosomeHistogram);


            foreach (var geneToConsider in genesToConsider)
            {
                headersAndHistograms = new List<KeyValuePair<string, ASETools.PreBucketedHistogram>>();
                headersAndHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("not mutated to last homozygous", globalWithAndWithoutGeneDistanceSoftBeginningToLastHomozygous[geneToConsider][false]));
                headersAndHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("mutated to last homozygous", globalWithAndWithoutGeneDistanceSoftBeginningToLastHomozygous[geneToConsider][true]));
                headersAndHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("not mutated to first non-homozygous", globalWithAndWithoutGeneDistanceSoftBeginningToFirstNonhomozygous[geneToConsider][false]));
                headersAndHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("mutated to first non-homozygous", globalWithAndWithoutGeneDistanceSoftBeginningToFirstNonhomozygous[geneToConsider][true]));

                var chromosome = commonData.geneLocationInformation.genesByName[geneToConsider].chromosome;
                var chromosomeSize = ASETools.chromosomeSizesByName[chromosome].size;
                var minLocus = commonData.geneLocationInformation.genesByName[geneToConsider].minLocus;
                var maxLocus = commonData.geneLocationInformation.genesByName[geneToConsider].maxLocus;

                histogramsFile.WriteLine(geneToConsider + " at " + chromosome + ":" + ASETools.NumberWithCommas(minLocus) + "-" + ASETools.NumberWithCommas(maxLocus) + " from beginning.");
                ASETools.PreBucketedHistogram.WriteBatchOfHistogramCDFs(histogramsFile, headersAndHistograms);

                headersAndHistograms = new List<KeyValuePair<string, ASETools.PreBucketedHistogram>>();
                headersAndHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("not mutated to last homozygous", globalWithAndWithoutGeneDistanceSoftEndToLastHomozygous[geneToConsider][false]));
                headersAndHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("mutated to last homozygous", globalWithAndWithoutGeneDistanceSoftEndToLastHomozygous[geneToConsider][true]));
                headersAndHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("not mutated to first non-homozygous", globalWithAndWithoutGeneDistanceSoftEndToFirstNonhomozygous[geneToConsider][false]));
                headersAndHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("mutated to first non-homozygous", globalWithAndWithoutGeneDistanceSoftEndToFirstNonhomozygous[geneToConsider][true]));

                histogramsFile.WriteLine(geneToConsider + " at " + chromosome + ":" + ASETools.NumberWithCommas(minLocus) + "-" + ASETools.NumberWithCommas(maxLocus) + 
                    "(" + ASETools.NumberWithCommas(chromosomeSize - maxLocus) + "-" + ASETools.NumberWithCommas(chromosomeSize - minLocus) + " from end).");
                ASETools.PreBucketedHistogram.WriteBatchOfHistogramCDFs(histogramsFile, headersAndHistograms);
            }

            histogramsFile.WriteLine();


            histogramsFile.WriteLine("Overall ASV");
            globalASVHistogram.WriteHistogram(histogramsFile);

            writePerDiseaseHistograms("Chromosome 6 by disease", histogramsFile, globalChr6ByDisease);

            histogramsFile.WriteLine();
            histogramsFile.WriteLine("Mean expression of more common allele at heterozygous germline sites that are still heterozygous in the tumor by age.");
            histogramsFile.WriteLine("age\tn\tmean\tstd deviation");
            for (int i = minAgeAtDiagnosis; i <= maxAgeAtDiagnosis; i++)
            {
                var x = globalExpressionOfMostCommonAlleleByAgeAtDiagnosis[i];
                histogramsFile.WriteLine(i + "\t" + x.getCount() + "\t" + x.getMeanAndStdDev().mean + "\t" + x.getMeanAndStdDev().stddev);
            }


            histogramsFile.WriteLine("**done**");
            histogramsFile.Close();

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(commonData.timer));
        } // Main

        static void writePerDiseaseHistograms(string headerString, StreamWriter histogramsFile, Dictionary<string, ASETools.PreBucketedHistogram> histograms)
        {
            histogramsFile.WriteLine(headerString);
            ASETools.PreBucketedHistogram.WriteBatchOfHistogramCDFs(histogramsFile, commonData.diseases.Select(_ => new KeyValuePair<string, ASETools.PreBucketedHistogram>(_, histograms[_])).ToList());
            histogramsFile.WriteLine();
        }

        static void writePerChromsomeHistograms(string headerString, StreamWriter histogramsFile, Dictionary<string, ASETools.PreBucketedHistogram> histograms)
        {
            histogramsFile.WriteLine(headerString);
            ASETools.PreBucketedHistogram.WriteBatchOfHistogramCDFs(histogramsFile, chromosomesToUse.Select(_ => new KeyValuePair<string, ASETools.PreBucketedHistogram>(_, histograms[_])).ToList());
            histogramsFile.WriteLine();
        }
        static void HandleOneCase(ASETools.Case case_, int state)
        {
            var annotatedSelectedVariants = ASETools.AnnotatedVariant.readFile(case_.tentative_asv_without_cnvs_filename);
            if (annotatedSelectedVariants == null)
            {
                throw new Exception("Unable to read tentative annotated selected variants for case " + case_.case_id + " from file " + case_.annotated_selected_variants_filename);
            }

            var perGeneHistograms = new Dictionary<string, ASETools.PreBucketedHistogram>();
            var outputLines = new List<string>();

            var aseCandidates = annotatedSelectedVariants.Where(_ => !_.somaticMutation && _.IsASECandidate(true, null, commonData)).ToList();   // Don't need copy number, since our input is pre-filtered
            var nASECandidates = aseCandidates.Count();
            var totalASE = aseCandidates.Select(_ => _.tumorRNAReadCounts.AlleleSpecificValue()).Sum();

            foreach (var chromosome in chromosomesToUse)
            {
                if (chromosome == "chrx" && (!commonData.clinicalSummariesByPatientId.ContainsKey(case_.case_id) || commonData.clinicalSummariesByPatientId[case_.case_id].getGender() != "female"))
                {
                    //
                    // Skip X in men, since it's haploid.
                    //
                    continue;
                }



                var variantsForThisChromosome = annotatedSelectedVariants.Where(_ => _.contig.ToLower() == chromosome && !_.somaticMutation).ToList();
                var chromosomeSize = ASETools.chromosomeSizesByName[chromosome].size;



                var nVariants = variantsForThisChromosome.Count();
                var normalUseful = variantsForThisChromosome.Where(_ => _.normalDNAReadCounts.usefulReads() >= minReads).ToList();
                var tumorUseful = variantsForThisChromosome.Where(_ => GetReadCount(_).usefulReads() >= minReads).ToList(); // Exclude CNVs
                var nTumorUseful = tumorUseful.Count();

                if (nTumorUseful < 10)
                {
                    continue;
                }

                tumorUseful.Sort((x, y) => x.locus.CompareTo(y.locus));
       
                var tumor90PercentFraction = (double)tumorUseful.Where(_ => GetReadCount(_).AlleleSpecificValue() >= .8).Count() / nTumorUseful;

                //
                // I suspect there's some more C# way of doing this, but something that's not O(n^2) isn't occurring to me, so a loop it is.
                //
                var nTumor90PercentConsecutiveBeginning = nTumorUseful; // The loop falls through if it's all >=90%, so this is the default.
                double fractionOfBasesTumor90PercentConsecutiveBeginning = 1;
                for (int i = 0; i < nTumorUseful; i++)
                {
                    if (GetReadCount(tumorUseful[i]).AlleleSpecificValue() < 0.8)
                    {
                        nTumor90PercentConsecutiveBeginning = i;

                        if (i == 0)
                        {
                            fractionOfBasesTumor90PercentConsecutiveBeginning = 0;
                        } else
                        {
                            fractionOfBasesTumor90PercentConsecutiveBeginning = (double)tumorUseful[i - 1].locus / chromosomeSize;
                        }
                        break;
                    }
                }

                var nTumor90PercentConsecutiveEnd = nTumorUseful; // The loop falls through if it's all >=90%, so this is the default.
                double fractionOfBasesTumor90PercentConsecutiveEnd = 1;
                for (int i = nTumorUseful - 1; i >= 0; i--)
                {
                    if (GetReadCount(tumorUseful[i]).AlleleSpecificValue() >= 0.8)
                    {
                        continue;
                    }

                    int n = nTumorUseful - i - 1;
                    nTumor90PercentConsecutiveEnd = n;

                    if (i == nTumorUseful - 1)
                    {
                        fractionOfBasesTumor90PercentConsecutiveEnd = 0;
                    } else
                    {
                        fractionOfBasesTumor90PercentConsecutiveEnd = ((double)chromosomeSize - tumorUseful[i + 1].locus) / chromosomeSize;
                    }
                    break;
                }

                //
                // Now compute the soft ends.
                //
                int indexOfSoftBeginningInclusive = -1;
                double totalASVSoFar = 0;
                for (int i = 0; i < nTumorUseful; i++)
                {
                    double asv = GetReadCount(tumorUseful[i]).AlleleSpecificValue();
                    totalASVSoFar += asv;
                    if (asv >= 0.8 && totalASVSoFar / (i + 1) >= 0.8)
                    {
                        indexOfSoftBeginningInclusive = i;
                    }
                }
                if (indexOfSoftBeginningInclusive == 0)
                {
                    indexOfSoftBeginningInclusive = -1; // Skip singletons at the end of the chromosome.
                }

                double distanceSoftBeginningInclusiveInMegabases;
                if (indexOfSoftBeginningInclusive == -1)
                {
                    distanceSoftBeginningInclusiveInMegabases = 0;
                } else
                {
                    distanceSoftBeginningInclusiveInMegabases = (double)tumorUseful[indexOfSoftBeginningInclusive].locus / ASETools.Million;
                }


                double distanceSoftBeginningFirstNonhomozygousMegabases;
                if (indexOfSoftBeginningInclusive >= nTumorUseful - 1)
                {
                    distanceSoftBeginningFirstNonhomozygousMegabases = (double)chromosomeSize / ASETools.Million;
                } else
                {
                    distanceSoftBeginningFirstNonhomozygousMegabases = (double)tumorUseful[indexOfSoftBeginningInclusive + 1].locus / ASETools.Million;
                }

                int indexOfSoftEndInclusive = nTumorUseful;
                totalASVSoFar = -1;
                for (int i = nTumorUseful - 1; i >= 0; i--)
                {
                    var asv = GetReadCount(tumorUseful[i]).AlleleSpecificValue();
                    totalASVSoFar += asv;
                    if (asv >= 0.8 && totalASVSoFar / (nTumorUseful - i) >= 0.8)
                    {
                        indexOfSoftEndInclusive = i;
                    }
                }
                if (indexOfSoftEndInclusive == nTumorUseful - 1)
                {
                    indexOfSoftEndInclusive = nTumorUseful;// Skip singletons at the end of the chromosome.
                }

                double distanceSoftEndInclusiveInMegabases;
                if (indexOfSoftEndInclusive == nTumorUseful)
                {
                    distanceSoftEndInclusiveInMegabases = 0;
                } else 
                {
                    distanceSoftEndInclusiveInMegabases = (double)(chromosomeSize - tumorUseful[indexOfSoftEndInclusive].locus) / ASETools.Million;
                }

                double distanceToSoftEndFirstNonHomozygousInMegabases;
                if (indexOfSoftEndInclusive <= 0)
                {
                    distanceToSoftEndFirstNonHomozygousInMegabases = (double)chromosomeSize / ASETools.Million;
                } else
                {
                    distanceToSoftEndFirstNonHomozygousInMegabases = ((double)chromosomeSize - tumorUseful[indexOfSoftEndInclusive - 1].locus) / ASETools.Million;
                }

                var flankingMutationsByGene = new Dictionary<string, int>();    // Maps gene->flanking mutations.  -1 means "not in this chromosome"

                foreach (var geneToConsider in genesToConsider)
                {
                    if (commonData.geneLocationInformation.genesByName[geneToConsider].chromosome != chromosome)
                    {
                        continue;
                    }

                    var minLocus = commonData.geneLocationInformation.genesByName[geneToConsider].minLocus;
                    bool previousMutationIsHomozygous = true;   // Start true, so that if the first variant site is beyond the start of the gene we treat it as double homozygous.
                    bool foundBothLoci = false;
                    for (int i = 0; i < nTumorUseful; i++)
                    {
                        if (tumorUseful[i].locus > minLocus)
                        {
                            int nHomozygous = GetReadCount(tumorUseful[i]).AlleleSpecificValue() >= 0.8 ? 1 : 0;
                            if (previousMutationIsHomozygous)
                            {
                                nHomozygous++;
                            }
                            flankingMutationsByGene.Add(geneToConsider, nHomozygous);
                            foundBothLoci = true;
                            break;
                        }

                        previousMutationIsHomozygous = GetReadCount(tumorUseful[i]).AlleleSpecificValue() >= 0.8;
                    }

                    if (!foundBothLoci)
                    {
                        //
                        // The gene is past the last germline variant site.  Treat it as double whatever the last one is.
                        //
                        if (previousMutationIsHomozygous)
                        {
                            flankingMutationsByGene.Add(geneToConsider, 2);
                        } else
                        {
                            flankingMutationsByGene.Add(geneToConsider, 0);
                        }
                    }
                } // genesToConsider


                var meanASV = tumorUseful.Select(_ => GetReadCount(_).AlleleSpecificValue()).Average();

                string outputLine = case_.case_id + "\t";
                if (commonData.clinicalSummariesByPatientId.ContainsKey(case_.case_id))
                {
                    var summary = commonData.clinicalSummariesByPatientId[case_.case_id];
                    outputLine += summary.getGender() + "\t" + summary.vitalStatus + "\t" + summary.age_at_diagnosis + "\t" + summary.days_to_death;
                } else
                {
                    outputLine += "unknown\tunknown\tunknown\tunknown";
                }

                outputLine += "\t" +
                    case_.disease() + "\t" + chromosome + "\t" + nVariants + "\t" + normalUseful.Count() + "\t" + nTumorUseful + "\t" +
                    ((double)normalUseful.Where(_ => _.normalDNAReadCounts.AlleleSpecificValue() >= .8).Count() / normalUseful.Count()) + "\t" + tumor90PercentFraction + "\t" +
                    fractionOfBasesTumor90PercentConsecutiveBeginning + "\t" + fractionOfBasesTumor90PercentConsecutiveEnd +
                    "\t" + nTumor90PercentConsecutiveBeginning + "\t" + nTumor90PercentConsecutiveEnd + "\t" + meanASV + "\t" + distanceSoftBeginningInclusiveInMegabases + "\t" + distanceSoftEndInclusiveInMegabases + "\t";

                foreach (var gene in genesToConsider)
                {
                    outputLine += annotatedSelectedVariants.Where(_ => _.somaticMutation && _.Hugo_symbol == gene).Count() + "\t";
                }
                outputLine += "=\"";
                for (int i = 0; i < nTumorUseful; i++)
                {
                    if (GetReadCount(tumorUseful[i]).AlleleSpecificValue() >= 0.8)
                    {
                        outputLine += "+";
                    } else
                    {
                        outputLine += "-";
                    }
                }
                outputLine += "\"\t=\"";

                var chromosomeSizeInMegabases = (ASETools.chromosomeSizesByName[chromosome].size + 999999) / 1000000;   // Rounded up
                var centromereInMegabases = (ASETools.chromosomeSizesByName[chromosome].centromere + 999999) / 1000000;   // Rounded up
                for (int i = 0; i < chromosomeSizeInMegabases; i++)
                {
                    if (i == centromereInMegabases)
                    {
                        outputLine += "C";
                    }

                    var inThisMegabase = tumorUseful.Where(_ => _.locus / 1000000 == i).ToList();

                    int n = inThisMegabase.Count();
                    if (n == 0)
                    {
                        outputLine += " ";
                        continue;
                    }

                    var meanInThisMegabase = inThisMegabase.Select(_ => GetReadCount(_).AlleleSpecificValue()).Average();
                    if (meanInThisMegabase < 0.4)
                    {
                        outputLine += "-";
                    } else if (meanInThisMegabase < 0.7)
                    {
                        outputLine += "~";
                    } else
                    {
                        outputLine += "+";
                    }
                } // for each megabase

                outputLine += "\"";

                outputLines.Add(outputLine);

                lock (globalTumorPerChromosomeHistograms)
                {
                    globalTumorPerChromosomeHistograms[chromosome].addValue(tumor90PercentFraction);
                    
                    if (commonData.clinicalSummariesByPatientId.ContainsKey(case_.case_id))
                    {
                        globalTumorHistogramByGender[commonData.clinicalSummariesByPatientId[case_.case_id].getGender()].addValue(tumor90PercentFraction);
                    }

                    globalBeginningOfChromosomeFractionByChromosome[chromosome].addValue(fractionOfBasesTumor90PercentConsecutiveBeginning);
                    globalEndOfChromosomeFractionByChromosome[chromosome].addValue(fractionOfBasesTumor90PercentConsecutiveEnd);
                    globalBeginningOfChromosomeFractionByDisease[case_.disease()].addValue(fractionOfBasesTumor90PercentConsecutiveBeginning);
                    globalEndOfChromosomeFractionByDisease[case_.disease()].addValue(fractionOfBasesTumor90PercentConsecutiveEnd);
                    globalBeginningOfChromosomeDistanceByChromosome[chromosome].addValue(fractionOfBasesTumor90PercentConsecutiveBeginning * chromosomeSize / ASETools.Million);
                    globalEndOfChromosomeDistanceByChromosome[chromosome].addValue(fractionOfBasesTumor90PercentConsecutiveEnd * chromosomeSize / ASETools.Million);
                    globalBeginningOfChromosomeDistanceByDisease[case_.disease()].addValue(fractionOfBasesTumor90PercentConsecutiveBeginning * chromosomeSize / ASETools.Million);
                    globalEndOfChromosomeDistanceByDisease[case_.disease()].addValue(fractionOfBasesTumor90PercentConsecutiveEnd * chromosomeSize / ASETools.Million);
                    if (nTumor90PercentConsecutiveBeginning == nTumorUseful)
                    {
                        // 
                        // a 100% homozygous chromosome
                        //
                        globalBeginningOfChromosomeDistanceToFirstNonHomozygousByChromosome[chromosome].addValue((double)ASETools.chromosomeSizesByName[chromosome].size / ASETools.Million);
                        globalEndOfChromosomeDistanceToFirstNonHomozygousByChromosome[chromosome].addValue((double)ASETools.chromosomeSizesByName[chromosome].size / ASETools.Million);
                    }
                    else
                    {
                        globalBeginningOfChromosomeDistanceToFirstNonHomozygousByChromosome[chromosome].addValue((double)tumorUseful[nTumor90PercentConsecutiveBeginning].locus / ASETools.Million);
                        globalEndOfChromosomeDistanceToFirstNonHomozygousByChromosome[chromosome].addValue(((double)ASETools.chromosomeSizesByName[chromosome].size - tumorUseful[nTumorUseful - nTumor90PercentConsecutiveBeginning -1].locus) / ASETools.Million);
                    }

                    globalPerDiseaseHistograms[case_.disease()].addValue(tumor90PercentFraction);
                    globalTumorHistogram.addValue(tumor90PercentFraction);

                    globalASVHistogram.addValue(meanASV);
                    globalMeanASVByChromosomeHistogram[chromosome].addValue(meanASV);

                    globalSoftDistanceFromBeginningOfChromosome[chromosome].addValue(distanceSoftBeginningInclusiveInMegabases);
                    globalBeginningOfChromosomeDistanceToFirstNonHomozygousByChromosomeSoft[chromosome].addValue(distanceSoftBeginningFirstNonhomozygousMegabases);
                    globalSoftDistanceFromEndOfChromosome[chromosome].addValue(distanceSoftEndInclusiveInMegabases);
                    globalEndOfChromosomeDistanceToFirstNonHomozygousByChromosomeSoft[chromosome].addValue(distanceToSoftEndFirstNonHomozygousInMegabases);

                    if (chromosome == "chr6")
                    {
                        globalChr6ByDisease[case_.disease()].addValue(tumor90PercentFraction);
                    }

                    foreach (var geneToConsider in genesToConsider)
                    {
                        if (chromosome != commonData.geneLocationInformation.genesByName[geneToConsider].chromosome)
                        {
                            continue;
                        }

                        var mutationCount = annotatedSelectedVariants.Where(_ => _.somaticMutation && _.Hugo_symbol == geneToConsider).Count();

                        globalWithAndWithoutGeneDistanceSoftBeginningToLastHomozygous[geneToConsider][mutationCount > 0].addValue(distanceSoftBeginningInclusiveInMegabases);
                        globalWithAndWithoutGeneDistanceSoftEndToLastHomozygous[geneToConsider][mutationCount > 0].addValue(distanceSoftEndInclusiveInMegabases);
                        globalWithAndWithoutGeneDistanceSoftBeginningToFirstNonhomozygous[geneToConsider][mutationCount > 0].addValue(distanceSoftBeginningFirstNonhomozygousMegabases);
                        globalWithAndWithoutGeneDistanceSoftEndToFirstNonhomozygous[geneToConsider][mutationCount > 0].addValue(distanceToSoftEndFirstNonHomozygousInMegabases);

                        if (flankingMutationsByGene.ContainsKey(geneToConsider))
                        {
                            globalGeneByFlankingHomozygousSites[geneToConsider][ASETools.ZeroOneMany(mutationCount)][flankingMutationsByGene[geneToConsider]]++;
                        }
                    }
                }
            } // foreach chromosome

            lock (globalPerGeneHistograms)
            {
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

                int age_at_diagnosis;
                if (nASECandidates != 0 && commonData.clinicalSummariesByPatientId.ContainsKey(case_.case_id) && (age_at_diagnosis = commonData.clinicalSummariesByPatientId[case_.case_id].age_at_diagnosis) != -1)
                {
                    globalExpressionOfMostCommonAlleleByAgeAtDiagnosis[age_at_diagnosis].addValue(.5 + totalASE / nASECandidates / 2);   // The math converts from ASV to fraction of most common allele.
                }

                foreach (var outputLine in outputLines)
                {
                    outputFile.WriteLine(outputLine);
                }
            } // lock
        } // HandleOneCase
    }
}
