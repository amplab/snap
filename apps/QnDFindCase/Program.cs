using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.Diagnostics;

namespace QnDFindCase
{
    class Program
    {

        static int tp53Min;
        static int tp53Max;
        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            var configuration = ASETools.Configuration.loadFromFile(args);
            if (null == configuration)
            {
                Console.WriteLine("Unable to load configuration.");
                return;
            }

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname).Select(x => x.Value).Where(x => x.annotated_selected_variants_filename != "").ToList();

            if (null == cases)
            {
                Console.WriteLine("Unable to load cases.");
                return;
            }

            var aseCorrection = ASETools.ASECorrection.LoadFromFile(configuration.finalResultsDirectory + ASETools.ASECorrectionFilename);
            if (null == aseCorrection)
            {
                Console.WriteLine("Unable to load ASE correction");
                return;
            }

            var geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(configuration.geneLocationInformationFilename));
            var geneMap = new ASETools.GeneMap(geneLocationInformation.genesByName);

            tp53Min = geneLocationInformation.genesByName["TP53"].minLocus;
            tp53Max = geneLocationInformation.genesByName["TP53"].maxLocus;

            var threading = new ASETools.ASVThreadingHelper<int>(cases, aseCorrection, (x, y) => true, handleOneCase, null, null, 100);

            Console.Write("Processing " + cases.Count() + " cases, 1 dot/100: ");

            threading.run();

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(timer));
        }  // Main
        static void handleOneCase(ASETools.Case case_, int state, List<ASETools.AnnotatedVariant> variantsForThisCase, ASETools.ASECorrection aseCorrection, Dictionary<bool, List<ASETools.CopyNumberVariation>> copyNumber)
        {
            if (variantsForThisCase.Where(x => x.somaticMutation && x.Hugo_symbol == "TP53" && x.CausesNonsenseMediatedDecay()).Count() == 0)
            {
                return;
            }

            if (variantsForThisCase.Where(x => !x.somaticMutation && x.contig == "chr17" && x.locus > tp53Min && x.locus < tp53Max).Count() == 0)
            {
                return;
            }

            Console.WriteLine("Case " + case_.case_id + " has both a nonsense mediated decay causing somatic mutation and a germline site in TP53.");
        } // handleOneCase

    } // Program


}
