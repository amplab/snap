using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ExpressionLib;
using System.IO;

namespace QuickNDirtyMakeExpressionScript
{
    class Program
    {
        static void Main(string[] args)
        {
            if (args.Count() != 3)
            {
                Console.WriteLine("usage: QuickNDirtyMakeExpressionScript germline_variant_filename mutant_expression_input_file_name expression_generation_script");
                return;
            }

            var inputFile = new StreamReader(args[0]);
            var mutantExpressionInputFile = new StreamWriter(args[1]);
            var samtoolsScriptFile = new StreamWriter(args[2]);

            var tcgaRecords = ExpressionTools.LoadTCGARecords(null, null, @"\\gcr\scratch\b99\bolosky\tcga-all.xml");
            ExpressionTools.LoadTCGARecordsForLocalRealigns(tcgaRecords, null, @"\\gcr\scratch\b99\bolosky\realigns.txt");
            ExpressionTools.LoadTCGAAdditionalMetadata(tcgaRecords, @"\\bolosky\f$\sequence\reads\tcga\tcgaAdditionalMetadata.txt");

            Dictionary<string, ExpressionTools.Sample> allSamples;
            var participants = ExpressionTools.BuildParticipantData(tcgaRecords, out allSamples, @"\\gcr\scratch\b99\bolosky\clinical");

            var experiments = ExpressionTools.LoadExperimentsFromFile(@"\\gcr\scratch\b99\bolosky\experiments.txt", participants, tcgaRecords);

            var experimentsByAnalysisID = new Dictionary<string, ExpressionTools.Experiment>();
            foreach (ExpressionTools.Experiment experiment in experiments)
            {
                experimentsByAnalysisID.Add(experiment.NormalDNAAnalysis.analysis_id, experiment);
                experimentsByAnalysisID.Add(experiment.TumorDNAAnalysis.analysis_id, experiment);
                experimentsByAnalysisID.Add(experiment.TumorRNAAnalysis.analysis_id, experiment);
                if (experiment.NormalRNAAnalysis != null)
                {
                    experimentsByAnalysisID.Add(experiment.NormalRNAAnalysis.analysis_id, experiment);
                }
            }

            //
            // Input file format is gene count followed by count VCF lines with the analysis ID prepended.
            //

            string line;
            while (null != (line = inputFile.ReadLine()))
            {
                var fields = line.Split('\t');
                if (fields.Count() != 2)
                {
                    Console.WriteLine("Bogus gene header " + line);
                    return;
                }

                string gene = fields[0];

                int n = Convert.ToInt32(fields[1]);

                for (int i = 0; i < n; i++ )
                {
                    var vcfPlusAnalysisLine = inputFile.ReadLine();
                    var vcfFields = vcfPlusAnalysisLine.Split('\t');

                    string analysisID = vcfFields[0];
                    if (!experimentsByAnalysisID.ContainsKey(analysisID))
                    {
                        Console.WriteLine("Can't find analysis " + analysisID + " in experiments");
                        continue;
                    }

                    var experiment = experimentsByAnalysisID[analysisID];
                    string participantID = experiment.participant.participantId;
                    string disease = experiment.disease_abbr;
                    string chromosomeWithoutChr; 
                    if (vcfFields[1].Count() < 3 || vcfFields[1].Substring(0,3).ToLower() != "chr")
                    {
                        chromosomeWithoutChr = vcfFields[1];
                    } else {
                        chromosomeWithoutChr = vcfFields[1].Substring(3);
                    }

                    int variantOffset = Convert.ToInt32(vcfFields[2]);
                    int startPos = Math.Max(0, variantOffset - 100);
                    int endPos = variantOffset + 100;

                    string dnaFilename = disease + @"\" + participantID + "-" + gene + "-chr" + chromosomeWithoutChr + "-" + startPos + "-" + endPos;
                    string rnaFilename = dnaFilename + "-RNA";

                    samtoolsScriptFile.WriteLine("samtools view " + experiment.TumorDNAAnalysis.bamFileName + " " + 
                        ExpressionTools.ChromPrefixFromRefassem(experiment.TumorDNAAnalysis.refassemShortName) + chromosomeWithoutChr + ":" + startPos + "-" + endPos +
                        " > " + dnaFilename);
                    samtoolsScriptFile.WriteLine("samtools view " + experiment.TumorRNAAnalysis.bamFileName + " " + 
                        ExpressionTools.ChromPrefixFromRefassem(experiment.TumorRNAAnalysis.refassemShortName) + chromosomeWithoutChr + ":" + startPos + "-" + endPos +
                        " > " + rnaFilename);

                    mutantExpressionInputFile.WriteLine(disease + "\t" + dnaFilename + "\t" + rnaFilename + "\t" + vcfPlusAnalysisLine + "\t" + experiment.TumorRNAAnalysis.refassemShortName);
                }
            }

            mutantExpressionInputFile.Close();
            samtoolsScriptFile.Close();
        }
    }
}
