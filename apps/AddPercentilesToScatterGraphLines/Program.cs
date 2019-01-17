using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;
using System.Diagnostics;

namespace AddPercentilesToScatterGraphLines
{
    class Program
    {
        static void Main(string[] args)
        {
            var commonData = ASETools.CommonData.LoadCommonData(args);
            if (null == commonData)
            {
                return;
            }

            if (commonData.configuration.commandLineArgs.Count() != 2)
            {
                Console.WriteLine("usage: AddPercentilesToScatterGraphLines chromosome disease");
                return;
            }

            var disease = commonData.configuration.commandLineArgs[1];
            var chromosome = ASETools.chromosomeNameToNonChrForm(commonData.configuration.commandLineArgs[0]);

            if (!commonData.expressionDistributionByChromosomeMap.map[disease].ContainsKey(chromosome) || commonData.expressionDistributionByChromosomeMap.map[disease][chromosome] == "" || !File.Exists(commonData.expressionDistributionByChromosomeMap.map[disease][chromosome]))
            {
                Console.WriteLine("Expression distribution by chromosome file doesn't exist.");
                return;
            }

            Console.Write("Loading scatter graph lines...");
            var scatterGraphLines = ASETools.GeneScatterGraphLine.LoadAllGeneScatterGraphLines(commonData.configuration.geneScatterGraphsDirectory, false, "*");
            Console.WriteLine(ASETools.ElapsedTimeInSeconds(commonData.timer));

            if (scatterGraphLines == null)
            {
                Console.WriteLine("Unable to load scatter graph lines");
                return;
            }

            var subTimer = new Stopwatch();
            subTimer.Start();
            Console.Write("Loading RNA mapped base counts...");
            var tumorRNAMappedBaseCounts = new Dictionary<string, ASETools.MappedBaseCount>();
            foreach (var case_ in commonData.listOfCases)
            {
                if (case_.tumor_rna_mapped_base_count_filename != "")
                {
                    tumorRNAMappedBaseCounts.Add(case_.case_id, ASETools.MappedBaseCount.readFromFile(case_.tumor_rna_mapped_base_count_filename));
                }
            }
            Console.WriteLine(ASETools.ElapsedTimeInSeconds(subTimer));

            var linesForThisDiseaseAndChromosome = scatterGraphLines.Where(_ => _.disease == disease && ASETools.chromosomeNameToNonChrForm(_.Chromosome) == ASETools.chromosomeNameToNonChrForm(chromosome)).ToList();
            var linesByStartPosition = linesForThisDiseaseAndChromosome.GroupByToDict(_ => _.Start_Position);

            var outputFile = ASETools.CreateStreamWriterWithRetry(commonData.configuration.geneScatterGraphsLinesWithPercentilesDirectory + ASETools.GeneScatterGraphLinesWithPercentilesPrefix + disease + "_" + ASETools.chromosomeNameToNonChrForm(chromosome));
            if (null == outputFile)
            {
                Console.WriteLine("Unable to open output file.");
                return;
            }

            int nBasesPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "bases", ASETools.chromosomeSizesByName[chromosome].size, out nBasesPerDot);
            int lastBaseProcessed = 0;

            outputFile.WriteLine(ASETools.GeneScatterGraphLine.percentileHeaderLine);

            var filename = commonData.expressionDistributionByChromosomeMap.map[disease][chromosome];
            ASETools.ExpressionDistribution.ReadFromFile(filename, _ => CheckPositionForLine(_, linesByStartPosition, outputFile, tumorRNAMappedBaseCounts, ref lastBaseProcessed, nBasesPerDot)); // Doing it this way avoids loading the whole thing into memory at once, which sometimes causes out-of-memory problems.

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine();
            Console.WriteLine("Total run time " + ASETools.ElapsedTimeInSeconds(commonData.timer));
        } // Main

        static void CheckPositionForLine(ASETools.ExpressionDistribution expressionDistribution, Dictionary<int, List<ASETools.GeneScatterGraphLine>> linesByStartPosition, StreamWriter outputFile, Dictionary<string, ASETools.MappedBaseCount> tumorRNAMappedBaseCounts, ref int lastBaseProcessed, int nBasesPerDot)
        {
            if (lastBaseProcessed / nBasesPerDot != expressionDistribution.locus / nBasesPerDot)
            {
                Console.Write(".");
            }

            lastBaseProcessed = expressionDistribution.locus;

            if (!linesByStartPosition.ContainsKey(expressionDistribution.locus))
            {
                return;
            }

            foreach (var scatterGraphLine in linesByStartPosition[expressionDistribution.locus])
            {
                if (scatterGraphLine.tumorRNAReadCounts != null && tumorRNAMappedBaseCounts.ContainsKey(scatterGraphLine.case_id))
                {
                    scatterGraphLine.tumorRNAFracAltPercentile = new double[11];
                    scatterGraphLine.tumorRNAFracRefPercentile = new double[11];
                    scatterGraphLine.tumorRNAFracAllPercentile = new double[11];

                    for (int j = 0; j < 11; j++)
                    {
                        if (expressionDistribution.getPercentile(j * 10) == 0)
                        {
                            scatterGraphLine.tumorRNAFracAltPercentile[j] = double.NegativeInfinity;
                            scatterGraphLine.tumorRNAFracRefPercentile[j] = double.NegativeInfinity;
                            scatterGraphLine.tumorRNAFracAllPercentile[j] = double.NegativeInfinity;
                        }
                        else
                        {
                            scatterGraphLine.tumorRNAFracAltPercentile[j] = ((double)scatterGraphLine.tumorRNAReadCounts.nMatchingAlt / tumorRNAMappedBaseCounts[scatterGraphLine.case_id].mappedBaseCount) / expressionDistribution.getPercentile(j * 10);
                            scatterGraphLine.tumorRNAFracRefPercentile[j] = ((double)scatterGraphLine.tumorRNAReadCounts.nMatchingReference / tumorRNAMappedBaseCounts[scatterGraphLine.case_id].mappedBaseCount) / expressionDistribution.getPercentile(j * 10);
                            scatterGraphLine.tumorRNAFracAllPercentile[j] = (((double)scatterGraphLine.tumorRNAReadCounts.totalReads()) / tumorRNAMappedBaseCounts[scatterGraphLine.case_id].mappedBaseCount) / expressionDistribution.getPercentile(j * 10);
                        }
                    }

                    outputFile.WriteLine(scatterGraphLine.toPercentileLine());
                }
            }// foreach line at this locus
        }// CheckPositionForLine
    } // Program
} // namespace AddPercentilesToScatterGraphs
