using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.Diagnostics;
using System.Threading;

namespace ComputeBasesInCodingAndKnownExpressionRegions
{
    class Program
    {
        static ASETools.CommonData commonData;
        static StreamWriter outputFile;
        static ASETools.ExpressionFile expressionFile;

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);

            if (null == commonData)
            {
                return;
            }

            if (commonData.configuration.commandLineArgs.Count() != 1)
            {
                Console.WriteLine("usage: ComputeBasesInCodingAndKnownExpressionRegions diseaseName");
                return;
            }

            var disease = commonData.configuration.commandLineArgs[0];
            if (!ASETools.GetListOfDiseases(commonData.cases).Contains(disease))
            {
                Console.WriteLine(disease + " does not appear to be a valid disease name.");
                return;
            }

            var selectedGenes = ASETools.SelectedGene.LoadFromFile(commonData.configuration.selectedGenesFilename);
            if (null == selectedGenes)
            {
                Console.WriteLine("Unable to load selected genes.");
                return;
            }

            selectedGenes = selectedGenes.Where(x => x.Hugo_Symbol != "Y_RNA").ToList();    // Y_RNA isn't really a gene, but there are many instances of it all over.  We should have excluded it earlier.

            var expressionFilename = commonData.configuration.expressionFilesDirectory + "expression_" + disease;

            Console.Write("Loading expression file...");
            try
            {
                expressionFile = new ASETools.ExpressionFile();
                expressionFile.LoadFromFile(expressionFilename);
            }
            catch
            {
                Console.WriteLine("Error loading expression file " + expressionFilename);
                return;
            }

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(commonData.timer));

            Directory.CreateDirectory(commonData.configuration.basesInKnownCodingRegionsDirectory); // This is a no-op if it already exists

            var outputFilename = commonData.configuration.basesInKnownCodingRegionsDirectory + ASETools.basesInKnownCodingRegionsPrefix + disease + ".txt";
            outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
            if (null == outputFile)
            {
                Console.WriteLine("Unable to open " + outputFilename);
                return;
            }
            outputFile.WriteLine("Hugo Symbol\tBases In Coding And Known Expression Region");

            Console.WriteLine("Processing " + selectedGenes.Count() + " genes, 1 dot/100 genes:");
            ASETools.PrintNumberBar(selectedGenes.Count() / 100);

            var threading = new ASETools.WorkerThreadHelper<ASETools.SelectedGene, int>(selectedGenes, HandleOneGene, null, null, 100);
            threading.run();

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(commonData.timer));
        } // Main

        static void HandleOneGene(ASETools.SelectedGene selectedGene, int state)
        {
            if (!commonData.geneLocationInformation.genesByName.ContainsKey(selectedGene.Hugo_Symbol))
            {
                return;
            }

            if (commonData.geneLocationInformation.genesByName[selectedGene.Hugo_Symbol].codingSizeOfLargestIsoform() > 1000000)
            {
                Console.WriteLine("Absurdly large gene " + selectedGene.Hugo_Symbol + ": " + commonData.geneLocationInformation.genesByName[selectedGene.Hugo_Symbol].codingSizeOfLargestIsoform());
            }

            var nBases = commonData.geneLocationInformation.genesByName[selectedGene.Hugo_Symbol].totalBasesInCodingAndKnownExpressionRegion(expressionFile);

            lock (outputFile)
            {
                outputFile.WriteLine(ASETools.ConvertToExcelString(selectedGene.Hugo_Symbol) + "\t" + nBases);
            }
        } // HandleOneGene
    }
}
