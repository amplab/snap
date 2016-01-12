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

        class MAFAndAssociatedVCFs : IComparable<MAFAndAssociatedVCFs>
        {
            public MAFAndAssociatedVCFs(ExpressionTools.MAFRecord maf_)
            {
                maf = maf_;
            }

            public int CompareTo(MAFAndAssociatedVCFs peer)
            {
                return maf.CompareTo(peer.maf);
            }

            public readonly ExpressionTools.MAFRecord maf;
            public List<string> VCFRecords = new List<string>();
        }

        const int MaxDistanceToRecordVariantsAroundMutations = 200;
        const int minReads = 30;    // There's a similar constant in GenerateScatterGraphs.cs


        class OddsAndOutputLine
        {
            public double odds = 0;
            public string vcfLine = null;
        }

        static void Main(string[] args)
        {
            if (args.Count() != 5 && false)
            {
                Console.WriteLine("usage: ProcessVCFs inputFile gencodeFile interestingGenesInputFile germlineCandidateOutputFile annotatedMAFOutputFile");
                Console.WriteLine("inputFile is a list of VCFs to process; they must all be hg18 or hg19 and correspond to the gencode file.");
                return;
            }

            var geneMap = new ExpressionTools.GeneMap(args[1]);
            var germlineCandidateOutputFile = new StreamWriter(args[3]);
            var annotatedMAFOutputFile = new StreamWriter(args[4]);

            var interestingGenes = new Dictionary<string, string>();
            var interestingGenesReader = new StreamReader(args[2]);
            string line;
            while (null != (line = interestingGenesReader.ReadLine()))
            {
                interestingGenes.Add(line, "");
            }
            interestingGenesReader.Close();
            interestingGenesReader = null;

            var germlineVariantsByGene = new Dictionary<string, Dictionary<string, OddsAndOutputLine>>();    // Maps gene name->(analysisID -> OddsAndOutputLine)


            var tcgaRecords = ExpressionTools.LoadTCGARecords(null, null, @"\\gcr\scratch\b99\bolosky\tcga-all.xml");
            ExpressionTools.LoadTCGARecordsForLocalRealigns(tcgaRecords, null, @"\\gcr\scratch\b99\bolosky\realigns.txt");
            ExpressionTools.LoadTCGAAdditionalMetadata(tcgaRecords, @"\\bolosky\f$\sequence\reads\tcga\tcgaAdditionalMetadata.txt");

            Dictionary<string, ExpressionTools.Sample> allSamples;
            var participants = ExpressionTools.BuildParticipantData(tcgaRecords, out allSamples, @"\\gcr\scratch\b99\bolosky\clinical");

            var sampleToParticipantIDMap = ExpressionTools.CreateSampleToParticipantIDMap(tcgaRecords);
            ExpressionTools.AddAllMAFFilesToParticipants(participants, sampleToParticipantIDMap, @"\\gcr\scratch\b99\bolosky\mafs\");
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
            var genesWithAtLeast500Mutations = new List<string>();
            foreach (var gene in mutationCountByGene.Where(n => (n.Value >= 500))) {
                genesWithAtLeast500Mutations.Add(gene.Key);
            }

            var comparisonMAFRecord = new ExpressionTools.MAFRecord();
            var comparisonAnnotatedMAFRecord = new MAFAndAssociatedVCFs(comparisonMAFRecord);

            int nSavedGermlineVariants = 0;
            int nMultipleGeneHits = 0;
            //
            // Now process the input file
            //
            var inputReader = new StreamReader(args[0]);
            string vcfFilename;
            while (null != (vcfFilename = inputReader.ReadLine()))
            {
                //
                // Each input line is a VCF filename, which is analysisID.vcf.
                //
                var pathComponents = vcfFilename.Split('\\');
                string analysisId = pathComponents[pathComponents.Count() - 1].Split('.')[0].ToLower();
                if (analysisId.Count() != 36)
                {
                    Console.WriteLine("Couldn't parse analysis ID from vcf filename " + vcfFilename + ", got " + analysisId + ", which is the wrong length.");
                    continue;
                }

                if (!tcgaRecords.ContainsKey(analysisId))
                {
                    Console.WriteLine("Can't find TCGA record for vcf " + vcfFilename);
                    continue;
                }

                var sampleID = tcgaRecords[analysisId].aliquot_id;
                if (!sampleToParticipantIDMap.ContainsKey(sampleID))
                {
                    Console.WriteLine("Can't map sample " + sampleID + " to participant.");
                }

                var participant = participants[sampleToParticipantIDMap[sampleID]];
                if (participant.mafs.Count() == 0)
                {
                    Console.WriteLine("Can't find MAFs for participant " + participant.participantId + ", vcf " + vcfFilename);
                    continue;
                }

                //
                // Build a list of MAF records for this participant to which we attach any associated VCFs.
                //
                var annotatedMAFs = new List<MAFAndAssociatedVCFs>();
                foreach (var maf in participant.mafs[0])
                {
                    annotatedMAFs.Add(new MAFAndAssociatedVCFs(maf));
                }
                annotatedMAFs.Sort();

                int nAnnotatedMAFs = annotatedMAFs.Count();

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

                int nVariants = 0;
                int hqVariants = 0;

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

                    nVariants++;

                    var fields = vcfLine.Split('\t');
                    if (fields.Count() < 8)
                    {
                        Console.WriteLine("Ill-formed line in vcf file " + vcfFilename + ", " + vcfLine);
                        error = true;
                        break;
                    }

                    string chromosome = fields[0].ToLower();
                    //
                    // The MAFs never have "chr" prefixes for chromosomes, but some of the reference genomes (and hence VCFs) do.  Strip it here.
                    //
                    if (chromosome.Count() > 3 && chromosome.Substring(0, 3) == "chr")
                    {
                        chromosome = chromosome.Substring(3);
                    }
                    int pos;
                    string reference = fields[3];   // ref is a C# keyword
                    string alt = fields[4];
                    string flagsField = fields[7];
                    string[] flags = flagsField.Split(';');

                    //
                    // The default values for all of the extracted data assure that if we don't see the data, we just reject the variant call.
                    //
                    double qual = 0;
                    int dp = 0;
                    int numalt = 0;
                    int srf = 0;
                    int srr = 0;
                    int saf = 0;
                    int sar = 0;
                    double mqm = 0;
                    double mqmr = 0;
                    double af = 0;
                    double odds = 0;


                    try
                    {
                        pos = Convert.ToInt32(fields[1]);
                        qual = Convert.ToDouble(fields[5]);

                        foreach (var flag in flags)
                        {
                            if (flag.IndexOf(',') != -1)
                            {
                                // It's an A style value, one for each alternate allele, and we don't care about these fields for loci with mutiple non-reference alleles.  So just quit before getting
                                // a format exception.
                                //
                                break;
                            }
                            var nameAndValue = flag.Split('=');
                            if (nameAndValue.Count() != 2)
                            {
                                Console.WriteLine("Unparsable flag " + flag + " in line " + vcfLine + " in file " + vcfFilename);
                                error = true;
                                break;
                            }

                            if (nameAndValue[0] == "DP")
                            {
                                dp = Convert.ToInt32(nameAndValue[1]);
                            }

                            if (nameAndValue[0] == "NUMALT")
                            {
                                numalt = Convert.ToInt32(nameAndValue[1]);
                                if (numalt != 1)
                                {
                                    break;  // We give up here, because some of the subsequent fields won't parse otherwise
                                }
                            }

                            if (nameAndValue[0] == "SRF")
                            {
                                srf = Convert.ToInt32(nameAndValue[1]);
                            }

                            if (nameAndValue[0] == "SRR")
                            {
                                srr = Convert.ToInt32(nameAndValue[1]);
                            }

                            if (nameAndValue[0] == "SAF")
                            {
                                saf = Convert.ToInt32(nameAndValue[1]);
                            }

                            if (nameAndValue[0] == "SAR")
                            {
                                sar = Convert.ToInt32(nameAndValue[1]);
                            }

                            if (nameAndValue[0] == "MQM")
                            {
                                mqm = Convert.ToDouble(nameAndValue[1]);
                            }

                            if (nameAndValue[0] == "MQMR")
                            {
                                mqmr = Convert.ToDouble(nameAndValue[1]);
                            }

                            if (nameAndValue[0] == "AF")
                            {
                                af = Convert.ToDouble(nameAndValue[1]);
                            }

                            if (nameAndValue[0] == "AF")
                            {
                                af = Convert.ToDouble(nameAndValue[1]);
                            }

                            if (nameAndValue[0] == "ODDS")
                            {
                                odds = Convert.ToDouble(nameAndValue[1]);
                            }
                        }

                        if (error)
                        {
                            break;
                        }

                    } catch (FormatException) {
                        Console.WriteLine("Couldn't parse line from vcf file " + vcfFilename + ", " + vcfLine);
                        error = true;
                        break;
                    }


                    //
                    // See which if any MAFS this applies to.
                    //
                    comparisonMAFRecord.Chrom = chromosome;
                    comparisonMAFRecord.Start_position = pos - MaxDistanceToRecordVariantsAroundMutations;  // This could conceivably be negative, which is still OK.

                    int index = annotatedMAFs.BinarySearch(comparisonAnnotatedMAFRecord);
                    if (index < 0)
                    {
                        index = ~index; // This is because BinarySearch returns the bitwise complement of the next larger element if there's not an exact match.
                    }

                    //
                    // If there is more than one exact match, there's no guarantee of which element we'll get, so back up index until it's no longer an exact match.
                    //
                    while (index > 0 && comparisonAnnotatedMAFRecord.CompareTo(annotatedMAFs[index - 1]) == 0)
                    {
                        index--;
                    }

                    //
                    // Index now points at the first MAF not too small to be a candidate for variant.
                    //
                    var genesCoveredByThisVariant = geneMap.GetGenesCovering(chromosome, pos);
                    int initialIndex = index;

                    comparisonAnnotatedMAFRecord.maf.Start_position += 2 * MaxDistanceToRecordVariantsAroundMutations;
                    while (index < nAnnotatedMAFs && annotatedMAFs[index].CompareTo(comparisonAnnotatedMAFRecord) < 1)
                    {
                        annotatedMAFs[index].VCFRecords.Add(vcfLine);
                        index++;
                    }

                    //
                    // The second thing we do with each variant is to see if it's a candidate for the germline variant allele-specific expression part of the study.
                    // For this, we exclude variants that are in the same genes as a somatic mutation, so that we don't wind up measuring allele specific expression
                    // selected for by an oncogenic somatic mutation rather than by the germline variant we're looking at.
                    //
                    double totalReferenceReads = srf + srr;
                    double totalAlternateReads = saf + sar;
                    if (qual < 30 || reference.Count() != 1 || alt.Count() != 1 || dp < minReads || 1 != numalt || srf < 2 || srr < 2 || saf < 2 || sar < 2 || 
                        totalReferenceReads / totalAlternateReads < .5 || totalAlternateReads / totalReferenceReads < .5 || mqm < 30 || mqmr < 30 || af != 0.5 || odds < 30.0)
                    {
                        //
                        // Only care about high quality SNVs.
                        //
                        continue;
                    }

                    index = initialIndex;
                    if (index > 0)
                    {
                        index--;    // Start one further back, because we might be in the middle of the gene.
                    }
                    bool somaticMutationInSameGene = false;

                    while (index < nAnnotatedMAFs && (index < initialIndex || annotatedMAFs[index].CompareTo(comparisonAnnotatedMAFRecord) < 1))
                    {
                        if (genesCoveredByThisVariant.Contains(annotatedMAFs[index].maf.Hugo_symbol))
                        {
                            somaticMutationInSameGene = true;
                            break;
                        }
                        index++;
                    }

                    if (somaticMutationInSameGene)
                    {
                        continue;
                    }

                    hqVariants++;

                    if (genesCoveredByThisVariant.Count() > 1)
                    {
                        nMultipleGeneHits++;
                    }

                    foreach (string gene_ in genesCoveredByThisVariant)
                    {
                        var gene = gene_.ToLower();
                        if (!interestingGenes.ContainsKey(gene))
                        {
                            continue;
                        }

                        if (!germlineVariantsByGene.ContainsKey(gene))
                        {
                            germlineVariantsByGene.Add(gene, new Dictionary<string, OddsAndOutputLine>());
                        }

                        if (!germlineVariantsByGene[gene].ContainsKey(analysisId))
                        {
                            germlineVariantsByGene[gene].Add(analysisId, new OddsAndOutputLine());
                        }

                        if (germlineVariantsByGene[gene][analysisId].odds < odds)
                        {
                            germlineVariantsByGene[gene][analysisId].odds = odds;
                            germlineVariantsByGene[gene][analysisId].vcfLine = vcfLine;
                        }
                        nSavedGermlineVariants++;
                    }
                } // each VCF line

                Console.Write(/*"(" + nVariants + "," + hqVariants + ")"*/ ".");

                vcfReader.Close();

                if (error)
                {
                    continue;   // Don't write an output for this file.
                }

                foreach (var annnotatedMAF in annotatedMAFs)
                {
                    annotatedMAFOutputFile.WriteLine(annnotatedMAF.maf.entire_maf_line);
                    foreach (var annotation in annnotatedMAF.VCFRecords)
                    {
                        annotatedMAFOutputFile.WriteLine("\t" + annotation);
                    }
                }

            } // while we have an input line (i.e., once around this loop for each VCF file)


            int nGenes = 0;
            int nWrittenVariants = 0;
            foreach (var variantByGene in germlineVariantsByGene)
            {
                nGenes++;
                nWrittenVariants += variantByGene.Value.Count();
                germlineCandidateOutputFile.WriteLine(variantByGene.Key + "\t" + variantByGene.Value.Count());
                foreach (var variant in variantByGene.Value)
                {
                    germlineCandidateOutputFile.WriteLine(variant.Key + "\t" + variant.Value.vcfLine);
                }
            }

            germlineCandidateOutputFile.Close();
            annotatedMAFOutputFile.Close();
            inputReader.Close();


        } // main
    }
}
