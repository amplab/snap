using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.Diagnostics;

namespace ExpressionByGene
{
    class Program
    {
        static ASETools.CommonData commonData;
        static ASETools.ExpressionFile expressionFile;
        static ASETools.BasesInCodingAndKnownExpressionRegions basesInCodingAndKnownExpressionRegions;

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);

            if (commonData == null)
            {
                //
                // It already printed an error message.
                //
                return;
            }

            if (commonData.configuration.commandLineArgs.Count() != 1)
            {
                Console.WriteLine("usage: ExpressionByGene disease_name");
                return;
            }

            var disease = commonData.configuration.commandLineArgs[0];

            basesInCodingAndKnownExpressionRegions = ASETools.BasesInCodingAndKnownExpressionRegions.LoadFromFile(commonData.configuration.basesInKnownCodingRegionsDirectory + ASETools.basesInKnownCodingRegionsPrefix + disease + ".txt");
            if (null == basesInCodingAndKnownExpressionRegions)
            {
                Console.WriteLine("Unable to load bases in known coding regions.");
                return;
            }

            var expressionFilename = commonData.configuration.expressionFilesDirectory + ASETools.Expression_filename_base + disease;
            if (!File.Exists(expressionFilename))
            {
                Console.WriteLine("expression file " + expressionFilename + " doesn't exist.");
                return;
            }

            Console.Write("Loading expression file...");
            try
            {
                expressionFile = new ASETools.ExpressionFile();
                expressionFile.LoadFromFile(expressionFilename);
            } catch
            {
                Console.WriteLine("Error loading expression file " + expressionFilename);
                return;
            }

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(commonData.timer));

            var selectedGenes = ASETools.SelectedGene.LoadFromFile(commonData.configuration.selectedGenesFilename).GroupByToDict(x => x.Hugo_Symbol);
 
            var casesToProcess = commonData.listOfCases.Where(x => x.disease() == disease && x.expression_by_gene_filename == "").ToList();
            int casesPerDot = casesToProcess.Count() > 100 ? 10 : 1;
            Console.WriteLine("Processing " + casesToProcess.Count() + " cases, 1 dot/" + casesPerDot + " cases:");
            ASETools.PrintNumberBar(casesToProcess.Count() / casesPerDot);

            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(casesToProcess, HandleOneCase, null, null, casesPerDot);

            threading.run();

            Console.WriteLine();
            Console.WriteLine(ASETools.ElapsedTimeInSeconds(commonData.timer));
        } // Main

        static void HandleOneCase(ASETools.Case case_, int state)
        {
            if (case_.tumor_rna_allcount_filename == "" || case_.tumor_rna_mapped_base_count_filename == "")
            {
                Console.WriteLine("Case " + case_.case_id + " doesn't have a tumor rna allcount file or tumor rna mapped base count file.  Skipping.");
                return;
            }

            long tumorMappedBaseCount = (ASETools.MappedBaseCount.readFromFile(case_.tumor_rna_mapped_base_count_filename).mappedBaseCount);

            var expressionByGene = new Dictionary<string, double>();
            foreach(var hugo_symbol in commonData.geneLocationInformation.genesByName.Select(x => x.Key).Where(x => basesInCodingAndKnownExpressionRegions.isGeneKnown(x)).ToList())
            {
                expressionByGene.Add(hugo_symbol, 0);
            }

            var allcountReader = new ASETools.AllcountReader(case_.tumor_rna_allcount_filename);

            if (!allcountReader.openFile())
            {
                Console.WriteLine("Unable to open allcount file " + case_.tumor_rna_allcount_filename);
                return;
            }

            allcountReader.ReadAllcountFile((x, y, z) => processBase(x, y, z, expressionByGene, tumorMappedBaseCount));

            var outputFilename = ASETools.GetDirectoryFromPathname(case_.tumor_rna_allcount_filename) + @"\" + case_.case_id + ASETools.expressionByGeneExtension;
            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
            if (null == outputFile)
            {
                Console.WriteLine("Unable to open output file " + outputFilename + ".  SKipping.");
                return;
            }

            outputFile.WriteLine("Hugo Symbol\tFraction of mean expression");

            foreach (var geneEntry in expressionByGene)
            {
                var hugo_symbol = geneEntry.Key;
                var totalExpression = geneEntry.Value;

                var geneLocationInformation = commonData.geneLocationInformation.genesByName[hugo_symbol];

                outputFile.WriteLine(ASETools.ConvertToExcelString(hugo_symbol) + "\t" + (totalExpression / basesInCodingAndKnownExpressionRegions.basesForGene(hugo_symbol)));
            }


            outputFile.WriteLine("**done**");
            outputFile.Close();
        } // HandleOneCase

        static void processBase(string contigName, int location, int currentMappedReadCount, Dictionary<string, double> expressionByGene, long mappedBaseCount)
        {
            ASETools.MeanAndStdDev meanAndStdDev;
            if (expressionFile.getValue(contigName, location, out meanAndStdDev))
            {
                if (meanAndStdDev.mean == 0)
                {
                    return;
                }

                foreach (var gene in commonData.geneMap.getGenesMappedTo(contigName, location))
                {
                    if (expressionByGene.ContainsKey(gene.hugoSymbol) && gene.containsLocusInCodingAndKnownExpressionRegion(contigName, location, expressionFile))
                    {
                        expressionByGene[gene.hugoSymbol] += ((double)currentMappedReadCount / mappedBaseCount) / meanAndStdDev.mean;
                    }
                } // foreach gene that covers this locus (and is selected)
            } // If we have expression values for this locus

        } // processBase
    } // Program
} // ExpressionByGene
