using ASELib;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ComputeGeneExpressionFraction
{
    internal class Program
    {
        static ASETools.CommonData commonData;

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);
            if (commonData == null)
            {
                return;
            }

            if (commonData.configuration.commandLineArgs.Any(_ => !commonData.cases.ContainsKey(_) && !_.ToLower().EndsWith(".gz")))
            {
                Console.WriteLine("At least one of the input args isn't a valid case ID or allcount.gz file");
                return;
            }

            int nItemsPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "Cases", commonData.configuration.commandLineArgs.Count(), out nItemsPerDot);

            var threading = new ASETools.WorkerThreadHelper<string, int>(commonData.configuration.commandLineArgs.ToList(), HandleOneCase, null, null, nItemsPerDot);
            threading.run();

            Console.WriteLine("Processed " + commonData.configuration.commandLineArgs.Count() + " cases in " + ASETools.ElapsedTime(commonData.timer));

        } // Main

        static Dictionary<string, double> fractionOfOverallExomeByGene = new Dictionary<string, double>();

        static void HandleOneCase(string caseIDOrAllcountFilename, int state)
        {
            string tumorAllcountFilename = "";
            string normalAllcountFilename = "";

            bool isAllcountFilename = caseIDOrAllcountFilename.ToLower().EndsWith(".gz");

            if (isAllcountFilename)
            {
                normalAllcountFilename = caseIDOrAllcountFilename;
            }
            else
            {
                var case_ = commonData.cases[caseIDOrAllcountFilename];

                if (case_.tumor_rna_allcount_filename == "")
                {
                    Console.WriteLine("case " + caseIDOrAllcountFilename + " doesn't have a tumor allcount file.  Skipping this case.");
                    return;
                }
                tumorAllcountFilename = case_.tumor_rna_allcount_filename;

                if (case_.normal_rna_file_id != "" && case_.normal_rna_allcount_filename == "")
                {
                    Console.WriteLine("Case " + caseIDOrAllcountFilename + " should have a normal RNA allcount file but doesn't.");
                    return;
                }
                normalAllcountFilename = case_.normal_rna_allcount_filename;
            }

            var result = new Dictionary<string, ASETools.GeneExpressionFraction>();
            commonData.geneLocationInformation.genesByName.Select(_ => _.Key).ToList().ForEach(_ => result.Add(_, new ASETools.GeneExpressionFraction(_, 0, 0, 0, 0)));
            var totalBasesByTumor = new Dictionary<bool, long>();

            foreach (var tumor in ASETools.BothBools) {
                totalBasesByTumor.Add(tumor, 0);

                var allcountFilename = tumor ? tumorAllcountFilename : normalAllcountFilename;
                if (allcountFilename == "") // This could be normal for TCGA or tumor for bulk normal marrow RNA
                {
                    continue;
                }

                var allcountReader = new ASETools.AllcountReader(allcountFilename);
                long mappedHQNuclearReads;
                int num_contigs;

                if (!allcountReader.openFile(out mappedHQNuclearReads, out num_contigs))
                {
                    Console.WriteLine("Failed to open allcount file " + allcountFilename);
                    return;
                }

                long totalBasesMapped = 0;

                allcountReader.ReadAllcountFile((x, y, z) => processBase(result, tumor, ref totalBasesMapped, x, y, z));
                allcountReader.Close();

                totalBasesByTumor[tumor] = totalBasesMapped;
            }

            string outputFilename;
            if (isAllcountFilename)
            {
                outputFilename = ASETools.GetDirectoryFromPathname(caseIDOrAllcountFilename) + @"\" + ASETools.GetFileNameFromPathname(caseIDOrAllcountFilename).Substring(0, caseIDOrAllcountFilename.IndexOf('.')) + ASETools.geneExpressionFractionExtension;
            } else
            {
                var case_ = commonData.cases[caseIDOrAllcountFilename];
                outputFilename = ASETools.GetDirectoryFromPathname(tumorAllcountFilename) + @"\" + caseIDOrAllcountFilename + ASETools.geneExpressionFractionExtension;
            }

            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
            if (outputFile == null)
            {
                Console.WriteLine("Unable to open output file " + outputFilename);
                return;
            }

            outputFile.WriteLine("Hugo symbol\tTumor Expression Fraction\tNormal Expression Fraction\tNormalized Tumor Fraction\tNormalized Normal Fraction"); // normal will be 0 if there is no normal RNA
            foreach (var geneExpressionFraction in result.Select(_ => _.Value))
            {
                outputFile.WriteLine(geneExpressionFraction.Hugo_symbol + "\t" + geneExpressionFraction.expressionByTumor[true] / totalBasesByTumor[true] + "\t" + 
                    geneExpressionFraction.expressionByTumor[false] / totalBasesByTumor[false] + "\t" + 
                    geneExpressionFraction.expressionByTumor[true] / totalBasesByTumor[true] / commonData.geneLocationInformation.genesByName[geneExpressionFraction.Hugo_symbol].codingSizeOfLargestIsoform() + "\t" +
                    geneExpressionFraction.expressionByTumor[false] / totalBasesByTumor[false] / commonData.geneLocationInformation.genesByName[geneExpressionFraction.Hugo_symbol].codingSizeOfLargestIsoform());
            }

            outputFile.WriteLine("**done**");
            outputFile.Close();

        } // HandleOneCase

        static void processBase(Dictionary<string, ASETools.GeneExpressionFraction> result, bool tumor, ref long totalBasesMapped, string contigName, int location, int currentMappedReadCount)
        {
            totalBasesMapped += currentMappedReadCount;

            foreach (var geneInfo in commonData.geneMap.getGenesMappedTo(contigName, location))
            {
                result[geneInfo.hugoSymbol].expressionByTumor[tumor] += currentMappedReadCount;
            } // foreach gene mapped here
        } // processBase


    } // Program
} // namespace
