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

	enum VariantType {
		somatic,
		germline,
		both
	}

    class Program
    {
        static ASETools.Genome genome = new ASETools.Genome();

        static int failedCases = 0;
        static int nCasesProcessed = 0;
        static int nCasesToProcess;

        static void ProcessCases(List<ASETools.Case> casesToProcess, ASETools.Configuration configuration, VariantType variantType)
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

                var timer = new Stopwatch();
                timer.Start();

                var annotatedVariants = new List<ASETools.AnnotatedVariant>();

                try
                {
                    if (!variantType.Equals(VariantType.germline))
                    {
                        var mafLines = ASETools.MAFLine.ReadFile(case_.extracted_maf_lines_filename, case_.maf_file_id, false);

                        if (null == mafLines)
                        {
                            Console.WriteLine("Case " + case_.case_id + " failed to load extracted MAF lines.  Ignoring.");
                            continue;
                        }

                        foreach (var mafLine in mafLines)
                        {
                            try
                            {
                                annotatedVariants.Add(AnnotateVariant(mafLine.Hugo_Symbol, true, case_, mafLine.Chromosome, 
                                                                        mafLine.Start_Position, mafLine.Match_Norm_Seq_Allele1, mafLine.Tumor_Seq_Allele2, 
                                                                        mafLine.Variant_Type, mafLine.Variant_Classification, mafLine.getExtractedReadsExtension()));
                            } catch (FileNotFoundException f)
                            {
                                //
                                // This can happen when there are so many reads that it's > 2GB.  We just skip annotating those variants.  Unfortunately, it can also happen
                                // when the file actually doesn't exist.  So, check to see if it's the first one.
                                //
                                if (annotatedVariants.Count() == 0)
                                {
                                    throw f;
                                } else
                                {
                                    Console.WriteLine("Skipping annotating somatic variant that had too many reads.  Case " + case_.case_id + ", hugo symbol " + mafLine.Hugo_Symbol + ", chrom " + mafLine.Chromosome + " start " + mafLine.Start_Position);
                                    continue;
                                }
                            }
                        }
                    } // somatic variants

                    int nSomaticVariants = annotatedVariants.Count();

                    if (!variantType.Equals(VariantType.somatic))
                    {
                        var selectedVariants = ASETools.SelectedVariant.LoadFromFile(case_.selected_variants_filename);
                        if (null == selectedVariants)
                        {
                            Console.WriteLine("Case " + case_.case_id + " failed to load selected variants.  Ignoring.");
                            continue;
                        }

                        foreach (var selectedVariant in selectedVariants)
                        {
                            try
                            {
                                annotatedVariants.Add(AnnotateVariant("", false, case_, selectedVariant.contig, selectedVariant.locus, Convert.ToString(selectedVariant.referenceBase), Convert.ToString(selectedVariant.altBase), "SNP", "",
                                    selectedVariant.getExtractedReadsExtension()));   // All of the variants we selected are SNPs, so it's just a constant
                            } catch (FileNotFoundException f)
                            {
                                //
                                // This can happen when there are so many reads that it's > 2GB.  We just skip annotating those variants.  Unfortunately, it can also happen
                                // when the file actually doesn't exist.  So, check to see if it's the first one.
                                //
                                if (annotatedVariants.Count() == nSomaticVariants)
                                {
                                    throw f;
                                }
                                else
                                {
                                    Console.WriteLine("Skipping annotating germline variant that had too many reads.  Case " + case_.case_id + ", chrom " + selectedVariant.contig + " start " + selectedVariant.locus);
                                    continue;
                                }
                            }
                        }
                    } // germline variants
                } catch {
                    Console.WriteLine("Error handling case " + case_.case_id + ".  Skipping it.");
                    Interlocked.Increment(ref failedCases);
                    continue;
                }

				// write out annotated selected variants
                string outputFilename = ASETools.GetDirectoryFromPathname(case_.selected_variants_filename) + @"\" + case_.case_id + ASETools.annotatedSelectedVariantsExtension;

                lock (casesToProcess)
                {
                    nCasesProcessed++;
                    Console.WriteLine("Writing " + annotatedVariants.Count() + " annotated variants to file " + outputFilename + ", it took " + ASETools.ElapsedTimeInSeconds(timer) + " to compute. " + nCasesProcessed + "/" + nCasesToProcess);
                }
				ASETools.AnnotatedVariant.writeFile(outputFilename, annotatedVariants);

			} // while true (grabbing cases from the work queue)
        } // ProcessCases

        static ASETools.AnnotatedVariant AnnotateVariant(string Hugo_symbol, bool somatic, ASETools.Case case_, string contig, int start_position, string reference_allele, string alt_allele, string variantType, string variantClassification, string subfileExtension)
        {
            ASETools.ReadCounts tumorDNAReadCounts = ASETools.ReadCounts.ComputeReadCounts(case_.tumor_dna_reads_at_selected_variants_filename, contig, start_position, reference_allele, alt_allele, variantType, case_.tumor_dna_file_id + subfileExtension, genome);
            ASETools.ReadCounts tumorRNAReadCounts = ASETools.ReadCounts.ComputeReadCounts(case_.tumor_rna_reads_at_selected_variants_filename, contig, start_position, reference_allele, alt_allele, variantType, case_.tumor_rna_file_id + subfileExtension, genome);
            ASETools.ReadCounts normalDNAReadCounts = ASETools.ReadCounts.ComputeReadCounts(case_.normal_dna_reads_at_selected_variants_filename, contig, start_position, reference_allele, alt_allele, variantType, case_.normal_dna_file_id + subfileExtension, genome);

            if (tumorDNAReadCounts == null || tumorRNAReadCounts == null || normalDNAReadCounts == null)
            {
                Console.WriteLine("Failed to correctly compute annotated variant for case " + case_.case_id);
                throw new FileNotFoundException();
            }

            ASETools.ReadCounts normalRNAReadCounts;
            if (case_.normal_rna_reads_at_selected_variants_filename != "")
            {
                normalRNAReadCounts = ASETools.ReadCounts.ComputeReadCounts(case_.normal_rna_reads_at_selected_variants_filename, contig, start_position, reference_allele, alt_allele, variantType, case_.normal_rna_file_id + subfileExtension, genome);
                if (null == normalDNAReadCounts)
                {
                    Console.WriteLine("Failed to correctly compute annotated variant for normal RNA for case " + case_.case_id);
                }
            } else
            {
                normalRNAReadCounts = null;
            }

            return new ASETools.AnnotatedVariant(Hugo_symbol, somatic, contig, start_position, reference_allele, alt_allele, variantType, variantClassification, tumorDNAReadCounts, tumorRNAReadCounts, normalDNAReadCounts, normalRNAReadCounts);
        }

        //
        // The backward comparison is so the biggest ones go first.
        //
        static int CompareCasesBySelectedReadsFileSizesBackwards(ASETools.Case one, ASETools.Case two)
        {
            return (two.normal_dna_reads_at_selected_variants_size + two.tumor_dna_reads_at_selected_variants_size + two.normal_rna_reads_at_selected_variants_size + two.tumor_rna_reads_at_selected_variants_size).CompareTo(
                one.normal_dna_reads_at_selected_variants_size + one.tumor_dna_reads_at_selected_variants_size + one.normal_rna_reads_at_selected_variants_size + one.tumor_rna_reads_at_selected_variants_size);
        }

        static int Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

			var configuration = ASETools.Configuration.loadFromFile(args);

            if (null == configuration)
            {
                Console.WriteLine("Giving up because we were unable to load configuration.");
                return -1;
            }

			if (configuration.commandLineArgs.Count() < 1)
            {
                Console.WriteLine("usage: AnnotateVariants <-s|-g> <case_ids>");
                Console.WriteLine("-s means somatic, -g means germline. Absense of flag means both will be included.");
                return -1;
            }

            bool somaticVariants = configuration.commandLineArgs[0] == "-s";
			bool germlineVariants = configuration.commandLineArgs[0] == "-g";

			var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

            if (null == cases)
            {
                Console.WriteLine("Unable to load cases file " + configuration.casesFilePathname + ".  You must generate cases before annotating variants.");
            }

            var casesToProcess = new List<ASETools.Case>();

			var startArg = somaticVariants || germlineVariants ? 1 : 0;

            for (int i = startArg; i < configuration.commandLineArgs.Count(); i++)
            {
                if (!cases.ContainsKey(configuration.commandLineArgs[i]))
                {
                    Console.WriteLine(configuration.commandLineArgs[i] + " does not appear to be a case ID.  Ignoring.");
                } else
                {
                    casesToProcess.Add(cases[configuration.commandLineArgs[i]]);
                }
            }

			nCasesToProcess = casesToProcess.Count();

			// get reference 
			genome.load(configuration.indexDirectory);

			var variantType = VariantType.both;

			if (somaticVariants)
				variantType = VariantType.somatic;
			else if (germlineVariants)
				variantType = VariantType.germline;

			var threads = new List<Thread>();

            // Sort the cases to process in reverse order of total selected reads file size, so the slower ones go first in order to avoid 
            // sitting around waiting for one long run at the end.

            casesToProcess.Sort(CompareCasesBySelectedReadsFileSizesBackwards);


            for (int i = 0; i < Environment.ProcessorCount; i++)
			{
				threads.Add(new Thread(() => ProcessCases(casesToProcess, configuration, variantType)));
			}

			threads.ForEach(t => t.Start());
			threads.ForEach(t => t.Join());

			Console.WriteLine("Processed " + nCasesToProcess + " in " + ASETools.ElapsedTimeInSeconds(timer));
            if (failedCases > 0)
            {
                Console.WriteLine("" + failedCases + " of them failed.");
            }

            return failedCases;
        } // Main
    }
}
