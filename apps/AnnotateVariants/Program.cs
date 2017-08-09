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
                            annotatedVariants.Add(AnnotateVariant(mafLine.Hugo_Symbol, true, case_, mafLine.Chromosome, mafLine.Start_Position, mafLine.Match_Norm_Seq_Allele1, mafLine.Tumor_Seq_Allele2, mafLine.Variant_Type, mafLine.Variant_Classification, mafLine.getExtractedReadsExtension()));
                        }
                    } // somatic variants

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
                            annotatedVariants.Add(AnnotateVariant("", false, case_, selectedVariant.contig, selectedVariant.locus, Convert.ToString(selectedVariant.referenceBase), Convert.ToString(selectedVariant.altBase), "SNP", "",
                                selectedVariant.getExtractedReadsExtension()));   // All of the variants we selected are SNPs, so it's just a constant
                        }
                    } // germline variants
                } catch (FileNotFoundException) {
                    Console.WriteLine("Error handling case " + case_.case_id + ".  Skipping it.");
                    Interlocked.Increment(ref failedCases);
                    continue;
                }

				// write out annotated selected variants
                string outputFilename = ASETools.GetDirectoryFromPathname(case_.selected_variants_filename) + @"\" + case_.case_id + ASETools.annotatedSelectedVariantsExtension;

				Console.WriteLine("Writing " + annotatedVariants.Count() + " annotated variants to file " + outputFilename + ", it took " + ASETools.ElapsedTimeInSeconds(timer) + " to compute.");
				ASETools.AnnotatedVariant.writeFile(outputFilename, annotatedVariants);

			} // while true (grabbing cases from the work queue)
        } // ProcessCases

        static ASETools.AnnotatedVariant AnnotateVariant(string Hugo_symbol, bool somatic, ASETools.Case case_, string contig, int start_position, string reference_allele, string alt_allele, string variantType, string variantClassification, string subfileExtension)
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

            return new ASETools.AnnotatedVariant(Hugo_symbol, somatic, contig, start_position, reference_allele, alt_allele, variantType, variantClassification, tumorDNAReadCounts, tumorRNAReadCounts, normalDNAReadCounts, normalRNAReadCounts);
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
                Console.WriteLine("ComputeReadCounts: no subfile named " + subfileName + " in " + selectedReadsFilename);
                return null;
            }

			var padding = 20;
			int[] posArray = new int[padding];
			for (int i = 0; i < padding; i++)
			{
				posArray[i] = start_position - padding/2 + i;
			}

			// get reference and alt sequences
			var reference_bases = posArray
				.Select(r => new Tuple<int, char>(r, genome.getBase(contig, r))).ToDictionary(x => x.Item1, x => x.Item2);
			var alt_bases = posArray
				.Select(r => new Tuple<int, char>(r, genome.getBase(contig, r))).ToDictionary(x => x.Item1, x => x.Item2);

			// If the variant we are looking at is an INS, we have to shift the variant position 
			// by one to match the SamLine insertion location
			int start = variantType == "INS" ? start_position + 1 : start_position;

			// modify alt_bases based on VariantType
			if (variantType == "DEL")
			{
				// replace missing bases with "N" so we can do a direct comparison
				for (var i = 0; i < reference_allele.Length; i++)
				{
					alt_bases[start + i] = 'N';
				}
			}
			else if (variantType != "INS")
			{
				// for NP: replace bases for alt
				for (var i = 0; i < alt_allele.Length; i++)
				{
					alt_bases[start + i] = alt_allele[i];
				}
			}
			else
			{
				// insert INS into dictionary
				KeyValuePair<int, char>[] ins = new KeyValuePair<int, char>[alt_allele.Length];

				var chArr = alt_allele.ToCharArray();
				for (int i = 0; i < chArr.Length; i++)
				{
					ins[i] = new KeyValuePair<int, char>(i + start, chArr[i]);
				}

				alt_bases =
					alt_bases.Where(r => r.Key < start).ToArray().Concat(ins.ToArray())
					.Concat(alt_bases.Where(r => r.Key >= start).Select(r => new KeyValuePair<int, char>(r.Key + alt_allele.Length, r.Value)).ToArray())
					.ToDictionary(x => x.Key, x => x.Value);
			}

			string line;
			while (null != (line = subfileReader.ReadLine()))
			{
				ASETools.SAMLine samLine;

				try
				{
					samLine = new ASETools.SAMLine(line);
				}
				catch (FormatException)
				{
					Console.WriteLine("Unable to parse sam line in extracted reads subfile " + subfileName + ": " + line);
					subfileReader.Close();
					return null;
				}

				if (samLine.isUnmapped() || samLine.mapq < 10)
				{
					//
					// Probably half of a paired-end read with the other end mapped.  Ignore it.
					//
					continue;
				}

				if (contig != samLine.rname)
				{
					// Wrong contig. Ignore it.
					continue;
				}

				// make sure variant can be found on this SamLine
				if (!samLine.mappedBases.ContainsKey(start))
				{
					//
					// This read does not overlap the variant.  Ignore it.
					// The 'before' case differs between variant types, and must be dealt with on a casewise basis
					// below in the switch statement.
					//
					continue;
				}
				if (!(start > 1 + samLine.mappedBases.Keys.Min())  || !(start < samLine.mappedBases.Keys.Max() - 1))
				{
					//
					// There is no padding on the variant. Without it, its hard to make a match call. Ignore.
					//
					continue;
				}

				bool matchesRef = false;
				bool matchesAlt = false;

				switch (variantType)
				{
					case "INS":
					case "DEL":
					case "SNP":
					case "DNP":
					case "TNP":

						// match to ref on both sides of variant
						var refMatchesLeft = matchSequence(reference_bases, samLine.mappedBases, start, false);
						var refMatchesRight = matchSequence(reference_bases, samLine.mappedBases, start, true);

						// match to alt on both sides of variant
						var altMatchesLeft = matchSequence(alt_bases, samLine.mappedBases, start, false);
						var altMatchesRight = matchSequence(alt_bases, samLine.mappedBases, start, true);

						// match criteria: there are some matching bases on both sides, with a total >= 10
						matchesRef = refMatchesLeft > 0 && refMatchesRight > 0 && refMatchesLeft + refMatchesRight >= 10;
						matchesAlt = altMatchesLeft > 0 && altMatchesRight > 0 && altMatchesLeft + altMatchesRight >= 10;

						break;
					default:
						Console.WriteLine("Unknown variant type: " + variantType);
						subfileReader.Close();
						return null;
				}

				// increment matches
				if (matchesRef)
				{
					if (matchesAlt)
					{
						nMatchingBoth++;
					}
					else
					{
						nMatchingRef++;
					}
				}
				else if (matchesAlt)
				{
					nMatchingAlt++;
				}
				else
				{
					nMatchingNeither++;
				}
			} // foreach SAMLine

            subfileReader.Close();
            return new ASETools.ReadCounts(nMatchingRef, nMatchingAlt, nMatchingNeither, nMatchingBoth);
        }
		
		// recursively matches a sequence and counts number of matches until the sequence ends or a mismatch is found
		static int matchSequence(Dictionary<int, char> reference, Dictionary<int, char> sequence, int position, bool goRight, int matchCount = 0, int threshold = 10)
		{
			if (matchCount >= threshold) {
				// We have already seen enough bases that match. Return.
				return matchCount;
			}
			// base case 1: strings have been consumed
			if (!reference.ContainsKey(position) || !sequence.ContainsKey(position))
			{
				// Reached end of sequence. Return. 
				return matchCount;
			}

			bool match = reference[position] == sequence[position];
			if (match)
			{
				if (goRight)
					position += 1;
				else
					position -= 1;

				matchCount += 1;
				return matchSequence(reference, sequence, position, goRight, matchCount);
			}
			else
			{
				// base case 2: hit a mismatch
				return matchCount;
			}


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
                    Console.WriteLine(configuration.commandLineArgs[i] + " does not appear to be a caseID.  Ignoring.");
                } else
                {
                    casesToProcess.Add(cases[configuration.commandLineArgs[i]]);
                }
            }

			int nCasesToProcess = casesToProcess.Count();

			// get reference 
			genome.load(configuration.indexDirectory);

			var variantType = VariantType.both;

			if (somaticVariants)
				variantType = VariantType.somatic;
			else if (germlineVariants)
				variantType = VariantType.germline;

			var threads = new List<Thread>();

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
