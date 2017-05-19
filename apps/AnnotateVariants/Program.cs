using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.Diagnostics;
using System.Threading;

namespace AnnotateVariants
{
    class Program
    {
        static ASETools.Genome genome = new ASETools.Genome();

        static void ProcessCases(List<ASETools.Case> casesToProcess, ASETools.ASEConfirguation configuration)
        {
            while (true)
            {
                ASETools.Case case_ = null;

                lock (casesToProcess)
                {
                    if (casesToProcess.Count() == 0)
                    {
                        return;
                    }

                    case_ = casesToProcess[0];
                    casesToProcess.RemoveAt(0);
                }

                if (case_.extracted_maf_lines_filename == "" || case_.selected_variants_filename == "" || case_.tumor_dna_reads_at_selected_variants_filename == "" || case_.tumor_dna_reads_at_selected_variants_index_filename == "" ||
                    case_.tumor_rna_reads_at_selected_variants_filename == "" || case_.tumor_rna_reads_at_selected_variants_index_filename == "" || case_.normal_dna_reads_at_selected_variants_filename == "" || case_.normal_dna_reads_at_selected_variants_index_filename == "")
                {
                    Console.WriteLine("Case " + case_.case_id + " is missing some required data.  Ignoring.");
                    continue;
                }

                var annotatedVariants = new List<ASETools.AnnotatedVariant>();

                var mafLines = ASETools.MAFLine.ReadFile(case_.extracted_maf_lines_filename, case_.maf_file_id, false);

                if (null == mafLines)
                {
                    Console.WriteLine("Case " + case_.case_id + " failed to load extracted MAF lines.  Ignoring.");
                    continue;
                }

                foreach (var mafLine in mafLines)
                {
                    annotatedVariants.Add(AnnotateVariant(true, case_, mafLine.Chromosome, mafLine.Start_Position, mafLine.Match_Norm_Seq_Allele1, mafLine.Tumor_Seq_Allele2, mafLine.Variant_Type, mafLine.getExtractedReadsExtension()));
                }

                var selectedVariants = ASETools.SelectedVariant.LoadFromFile(case_.selected_variants_filename);
                if (null == selectedVariants)
                {
                    Console.WriteLine("Case " + case_.case_id + " failed to load selected variants.  Ignoring.");
                    continue;
                }

                foreach (var selectedVariant in selectedVariants)
                {
                    annotatedVariants.Add(AnnotateVariant(false, case_, selectedVariant.contig, selectedVariant.locus, Convert.ToString(selectedVariant.referenceBase), Convert.ToString(selectedVariant.altBase), "SNP",
                        selectedVariant.getExtractedReadsExtension()));   // All of the variants we selected are SNPs, so it's just a constant
                }

                string outputFilename = ASETools.GetDirectoryFromPathname(case_.selected_variants_filename) + @"\" + case_.case_id + ASETools.annotatedSelectedVariantsExtension;

            }

        }

        static ASETools.AnnotatedVariant AnnotateVariant(bool somatic, ASETools.Case case_, string contig, int start_position, string reference_allele, string alt_allele, string variantType, string subfileExtension)
        {
            ASETools.ReadCounts tumorDNAReadCounts = ComputeReadCounts(case_.tumor_dna_reads_at_selected_variants_filename, contig, start_position, reference_allele, alt_allele, variantType, case_.tumor_dna_file_id + subfileExtension);
            ASETools.ReadCounts tumorRNAReadCounts = ComputeReadCounts(case_.tumor_rna_reads_at_selected_variants_filename, contig, start_position, reference_allele, alt_allele, variantType, case_.tumor_rna_file_id + subfileExtension);
            ASETools.ReadCounts normalDNAReadCounts = ComputeReadCounts(case_.normal_dna_reads_at_selected_variants_filename, contig, start_position, reference_allele, alt_allele, variantType, case_.normal_dna_file_id + subfileExtension);

            if (tumorDNAReadCounts == null || tumorRNAReadCounts == null || normalDNAReadCounts == null)
            {
                Console.WriteLine("Failed to correctly compute annotated variant for case " + case_.case_id);
                throw new FileNotFoundException();
            }

            ASETools.ReadCounts normalRNAReadCounts;
            if (case_.normal_rna_reads_at_selected_variants_filename != "")
            {
                normalRNAReadCounts = ComputeReadCounts(case_.normal_rna_reads_at_selected_variants_filename, contig, start_position, reference_allele, alt_allele, variantType, case_.normal_rna_file_id + subfileExtension);
                if (null == normalDNAReadCounts)
                {
                    Console.WriteLine("Failed to correctly compute annotated variant for normal RNA for case " + case_.case_id);
                }
            } else
            {
                normalRNAReadCounts = null;
            }

            return new ASETools.AnnotatedVariant(somatic, contig, start_position, reference_allele, alt_allele, variantType, tumorDNAReadCounts, tumorRNAReadCounts, normalDNAReadCounts, normalRNAReadCounts);
        }

        static ASETools.ReadCounts ComputeReadCounts(string selectedReadsFilename, string contig, int start_position, string reference_allele, string alt_allele, string variantType, string subfileName)
        {
            int nMatchingRef = 0;
            int nMatchingAlt = 0;
            int nMatchingNeither = 0;
            int nMatchingBoth = 0;

            var consolodatedFile = new ASETools.ConsolodatedFileReader();

            if (!consolodatedFile.open(selectedReadsFilename))
            {
                Console.WriteLine("Unable to open reads at selected variants file " + selectedReadsFilename);
                return null;
            }

            var subfileReader = consolodatedFile.getSubfile(subfileName);
            if (null == subfileReader)
            {
                Console.WriteLine("ComputeReadCounts: no subfile named " + subfileName);
                return null;
            }

            string line;
            while (null != (line = subfileReader.ReadLine()))
            {
                ASETools.SAMLine samLine;

                try
                {
                    samLine = new ASETools.SAMLine(line);
                } catch (FormatException)
                {
                    Console.WriteLine("Unable to parse sam line in extracted reads subfile " + subfileName + ": " + line);
                    subfileReader.Close();
                    return null;
                }

                if (samLine.isUnmapped())
                {
                    //
                    // Probably half of a paired-end read with the other end mapped.  Ignore it.
                    //
                    continue;
                }

                switch (variantType)
                {
                    case "SNP":


                        break;


                    default:
                        Console.WriteLine("Unknown variant type: " + variantType);
                        subfileReader.Close();
                        return null;
                }
            }

            subfileReader.Close();
            return new ASETools.ReadCounts(nMatchingRef, nMatchingAlt, nMatchingNeither, nMatchingBoth);
        }

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            var configuration = ASETools.ASEConfirguation.loadFromFile(args);

            if (null == configuration)
            {
                Console.WriteLine("Giving up because we were unable to load configuration.");
                return;
            }

            if (configuration.commandLineArgs.Count() <2 || configuration.commandLineArgs[0] != "-s" && configuration.commandLineArgs[0] != "-g")
            {
                Console.WriteLine("usage: AnnotateVariants {-s|-g} <case_ids>");
                Console.WriteLine("-s means semantic and -g means germline");
                return;
            }

            bool sematicVariants = configuration.commandLineArgs[0] == "-s";

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

            if (null == cases)
            {
                Console.WriteLine("Unable to load cases file " + configuration.casesFilePathname + ".  You must generate cases before annotating variants.");
            }

            var casesToProcess = new List<ASETools.Case>();
            int nCasesToProcess = casesToProcess.Count();

            for (int i = 1; i < configuration.commandLineArgs.Count(); i++)
            {
                if (!cases.ContainsKey(configuration.commandLineArgs[i]))
                {
                    Console.WriteLine(configuration.commandLineArgs[i] + " does not appear to be a case ID.  Ignoring.");
                } else
                {
                    casesToProcess.Add(cases[configuration.commandLineArgs[i]]);
                }
            }

            genome.load(configuration.indexDirectory);

            var threads = new List<Thread>();

            for (int i = 0; i < Environment.ProcessorCount; i++)
            {
                threads.Add(new Thread(() => ProcessCases(casesToProcess, configuration)));
            }

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            Console.WriteLine("Processed " + nCasesToProcess + " in " + ASETools.ElapsedTimeInSeconds(timer));
        } // Main
    }
}
