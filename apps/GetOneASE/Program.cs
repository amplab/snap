using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;

namespace GetOneASE
{
    class Program
    {
        static void Main(string[] args)
        {
            var configuration = ASETools.Configuration.loadFromFile(args);

            if (null == configuration)
            {
                Console.WriteLine("Giving up because we were unable to load configuration.");
                return;
            }

            if (configuration.commandLineArgs.Count() < 2 || configuration.commandLineArgs.Count() % 2 != 0)
            {
                Console.WriteLine("usage: GetOneASE {caseId chromosome}");
                return;
            }

            var geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(configuration.geneLocationInformationFilename));
            var geneMap = new ASETools.GeneMap(geneLocationInformation.genesByName);
            var perGeneASEMap = ASETools.ASEMapPerGeneLine.ReadFromFileToDictionary(configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

            if (null == cases)
            {
                Console.WriteLine("Unable to load cases.");
            }

            Console.WriteLine("CaseId\tContig\ttumor read count\tuseful tumor read count\tlocus\tASE");

            for (int i = 0; i < configuration.commandLineArgs.Count(); i += 2)
            {
                if (!cases.ContainsKey(args[i]))
                {
                    Console.WriteLine("Couldn't find case " + args[i]);
                    return;
                }

                var case_ = cases[args[i]];

                if (case_.annotated_selected_variants_filename == "" || case_.tumor_copy_number_filename == "")
                {
                    Console.WriteLine("Case is missing annotated selected variants and/or tumor copy number");
                }

                var annotatedSelectedVariants = ASETools.AnnotatedVariant.readFile(case_.annotated_selected_variants_filename);
                var tumorCopyNumberVariation = ASETools.CopyNumberVariation.ReadFile(case_.tumor_copy_number_filename, case_.tumor_copy_number_file_id).Where(r => Math.Abs(r.Segment_Mean) > 1.0).ToList();

                annotatedSelectedVariants = annotatedSelectedVariants.Where(x => x.contig == args[i + 1] && x.IsASECandidate(true, tumorCopyNumberVariation, configuration, perGeneASEMap, geneMap) && !x.somaticMutation).ToList();
                annotatedSelectedVariants.Sort();

                foreach (var variant in annotatedSelectedVariants)
                {
                    Console.WriteLine(case_.case_id + "\t" + variant.contig + "\t" + variant.tumorDNAReadCounts.totalReads() + "\t" + variant.tumorDNAReadCounts.usefulReads() + "\t" + variant.locus + "\t" + variant.GetAlleleSpecificExpression(true));
                }
                Console.WriteLine();
            } // for each input set
        }
    }
}
