using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ExpressionLib;

namespace ProcessVCFs
{
    class Program
    {
        static void Main(string[] args)
        {
            if (args.Count() != 3 && false)
            {
                Console.WriteLine("usage: ProcessVCFs inputFile gencodeFile outputFile");
                Console.WriteLine("inputFile is a list of VCFs to process;they must all be hg18 or hg19 and correspond to the gencode file.");
                return;
            }

            var outputFile = new StreamWriter(args[2]);
            var geneMap = new ExpressionTools.GeneMap(args[1]);


            var tcgaRecords = ExpressionTools.LoadTCGARecords(null, null, @"\\gcr\scratch\b99\bolosky\tcga-all.xml");
            ExpressionTools.LoadTCGARecordsForLocalRealigns(tcgaRecords, null, @"\\gcr\scratch\b99\bolosky\realigns.txt");
            ExpressionTools.LoadTCGAAdditionalMetadata(tcgaRecords);

            Dictionary<string, ExpressionTools.Sample> allSamples;
            var participants = ExpressionTools.BuildParticipantData(tcgaRecords, out allSamples);

            var sampleToParticipantIDMap = ExpressionTools.CreateSampleToParticipantIDMap(tcgaRecords);
            ExpressionTools.AddAllMAFFilesToParticipants(participants, sampleToParticipantIDMap);
            var experiments = ExpressionTools.LoadExperimentsFromFile(@"\\gcr\scratch\b99\bolosky\experiments.txt" ,participants, tcgaRecords);

            //
            // Make a list of genes for which there are any mutations.
            //
            var mutationCountByGene = new Dictionary<string, int>();
            foreach (var experiment in experiments)
            {
                foreach (var maf in experiment.maf) {
                    string hugo_symbol = maf.Hugo_symbol.ToLower();
                    if (!mutationCountByGene.ContainsKey(hugo_symbol))
                    {
                        mutationCountByGene.Add(hugo_symbol, 0);
                    }

                    mutationCountByGene[hugo_symbol]++;
                }
            }

            Console.WriteLine("Have " + mutationCountByGene.Count() + " genes with at least one mutation and " + mutationCountByGene.Where(n => (n.Value >= 10)).Count() + " with at least 10, and " +
                mutationCountByGene.Where(n => (n.Value >= 100)).Count() + " with at least 100.");

            //
            // And one for genes with at least 100 mutations (which we'll use to filter out the noise).
            //
            var genesWithAtLeast100Mutations = mutationCountByGene.Where(n => (n.Value >= 100));

            //
            // Now process the input file
            //
            var inputReader = new StreamReader(args[0]);
            string vcfFilename;
            while (null != (vcfFilename = inputReader.ReadLine()))
            {
                //
                // Each input line is just a VCF filename.
                //
                StreamReader vcfReader = null;
                try
                {
                    vcfReader = new StreamReader(vcfFilename);
                }
                catch (FileNotFoundException)
                {
                    Console.WriteLine("Unable to open VCF file " + vcfFilename + ", skipping");
                    continue;
                }

                string vcfLine;
                bool error = false;
                while (null != (vcfLine = vcfReader.ReadLine()))
                {
                    if (vcfLine.Count() == 0 || vcfLine[0] == '#')
                    {
                        //
                        // Skip this blank or comment line.
                        //
                        continue;
                    }

                    var fields = vcfLine.Split('\t');
                    if (fields.Count() < 6)
                    {
                        Console.WriteLine("Ill-formed line in vcf file " + vcfFilename + ", " + vcfLine);
                        error = true;
                        break;
                    }

                    string chromosome = fields[0];
                    int pos;
                    string reference = fields[3];   // ref is a C# keyword
                    string alt = fields[4];
                    double qual;

                    try
                    {
                        pos = Convert.ToInt32(fields[1]);
                        qual = Convert.ToDouble(fields[5]);
                    } catch (FormatException) {
                        Console.WriteLine("Couldn't parse line from vcf file " + vcfFilename + ", " + vcfLine);
                        error = true;
                        break;
                    }

                    if (qual < 30 || reference.Count() != 1 || alt.Count() != 1)
                    {
                        //
                        // Only care about high quality SNVs
                        continue;
                    }

                }

                if (error)
                {
                    continue;   // Don't write an output for this file.
                }

            } // while we have an input line


            outputFile.Close();
            inputReader.Close();
        } // main
    }
}
