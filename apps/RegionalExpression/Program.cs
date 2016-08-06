using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ExpressionLib;


namespace RegionalExpression
{
    class Program
    {
        class MutatedGene
        {
            public string chromosome = null;
            public int minOffset = 0x7fffffff;
            public int maxOffset = 0;
            public int mutationCount = 0;
            public string geneName = null;
            public bool broken = false;
        };

        static void Main(string[] args)
        {

            string sourceShare = @"\\msr-genomics-0\d$\";

            if (args.Count() < 1) {
                Console.WriteLine("usage: RegionalExpression <any number of participantIDs");
                return;
            }

            var tcgaRecords = ExpressionTools.LoadTCGARecords(null, null, @"\\gcr\scratch\b99\bolosky\tcga-all.xml");
            ExpressionTools.LoadTCGARecordsForLocalRealigns(tcgaRecords, null, @"\\gcr\scratch\b99\bolosky\realigns.txt");
            ExpressionTools.LoadTCGAAdditionalMetadata(tcgaRecords, @"\\bolosky\f$\sequence\reads\tcga\tcgaAdditionalMetadata.txt");

            Dictionary<string, ExpressionTools.Sample> allSamples;
            var participants = ExpressionTools.BuildParticipantData(tcgaRecords, out allSamples, @"\\gcr\scratch\b99\bolosky\clinical");

            var sampleToParticipantIDMap = ExpressionTools.CreateSampleToParticipantIDMap(tcgaRecords);
            ExpressionTools.AddAllMAFFilesToParticipants(participants, sampleToParticipantIDMap, @"\\gcr\scratch\b99\bolosky\mafs\");
            var experiments = ExpressionTools.LoadExperimentsFromFile(@"\\gcr\scratch\b99\bolosky\experiments.txt", participants, tcgaRecords);

            var participantToExperimentMapping = ExpressionTools.BuildParticipantToExperimentMapping(experiments);

            Dictionary<string, Dictionary<int, ExpressionTools.MeanAndStdDev>> expression = null;
            string expressionDisease = "";

            //
            // Run through the RNA and find mutations.
            //
            var mutationsByGene = new Dictionary<string, MutatedGene>();
            var mutationsByLocation = new Dictionary<string, List<MutatedGene>>();  // List is ordered by minOffset (eventually)

            var allTumorRNA = new StreamReader(sourceShare + @"sequence\reads\tcga\tumor_rna.txt");
            string line;
            while (null != (line = allTumorRNA.ReadLine()))
            {
                var fields = line.Split('\t');
                if (fields[0] == "Cancer_Type")
                {
                    //
                    // Skip the header
                    //
                    continue;
                }

                string cancerType = fields[0].ToLower();

                string variantClassification = fields[10];
                if (variantClassification.ToLower() == "silent")
                {
                    //
                    // Just ignore silent mutations.
                    //
                    continue;
                }

                string geneName = fields[2].ToLower();

                if (geneName == "unknown" || geneName == ".")
                {
                    continue;
                }

                string chrom = fields[6].ToLower();
                if (chrom.Length > 3 && chrom.Substring(0, 3) == "chr")
                {
                    chrom = chrom.Substring(3); // Strip the leading chr if it exists
                }

                if (!mutationsByGene.ContainsKey(geneName))
                {
                    mutationsByGene.Add(geneName, new MutatedGene());
                    mutationsByGene[geneName].geneName = geneName;
                    mutationsByGene[geneName].chromosome = chrom;
                    if (!mutationsByLocation.ContainsKey(chrom))
                    {
                        mutationsByLocation.Add(chrom, new List<MutatedGene>);
                    }
                    mutationsByLocation[chrom].Add(mutationsByGene[geneName]);
                }

                var gene = mutationsByGene[geneName];
                if (!gene.broken && gene.chromosome != chrom) {
                    gene.broken = true;
                    Console.WriteLine("Gene " + geneName + " appears in different chromosomes, " + gene.chromosome + " != " + chrom + ".  Ignoring the gene.");
                    continue;
                }

                int startPosition = Convert.ToInt32(fields[6]);
                int endPosition = Convert.ToInt32(fields[7]);

                if (startPosition > endPosition) {
                    Console.WriteLine("tumor_rna has line with start position " + startPosition + > end position; ignorning")
                }

                gene.maxOffset = Math.Max(gene.maxOffset, fields)

            }



            for (int i = 0; i < args.Count(); i++)  // for each person we're processing
            {
                ExpressionTools.Participant participant;
                if (!participants.ContainsKey(args[i]))
                {
                    Console.WriteLine("Couldn't find participant " + args[i] + ", are you sure it's a correct participant ID?");
                    return;
                }

                participant = participants[args[i]];
                var experiment = participantToExperimentMapping[participant.participantId];
                if (experiment.TumorRNAAnalysis == null || experiment.TumorRNAAnalysis.allcountFileName == null) {
                    Console.WriteLine("Participant " + participant.participantId + " doesn't appear to have an allcount file");
                    continue;
                }

                //
                // If we don't have the right expression data loaded, get it now.
                //
                if (experiment.disease_abbr != expressionDisease)
                {
                    expression = ExpressionTools.LoadExpressionFile(@"\\msr-genomics-0\d$\sequence\reads\tcga\expression\expression_" + experiment.disease_abbr);   // This can take forever!
                    expressionDisease = experiment.disease_abbr;
                }

                //
                // Build up a list of 


            }


        }
    }
}
