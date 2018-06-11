using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.Diagnostics;
using System.IO;

namespace QnDFindCase
{
    class Program
    {

        static int tp53Min;
        static int tp53Max;
        static ASETools.Configuration configuration;
        static ASETools.GeneMap geneMap;
        static Dictionary<string, ASETools.ASEMapPerGeneLine> perGeneASEMap;
        static StreamWriter outputFile;

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            configuration = ASETools.Configuration.loadFromFile(args);
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
            geneMap = new ASETools.GeneMap(geneLocationInformation.genesByName);

            perGeneASEMap = ASETools.ASEMapPerGeneLine.ReadFromFileToDictionary(configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);

            if (null == perGeneASEMap)
            {
                Console.WriteLine("You must first create the per-gene ASE map in " + configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);
                return;
            }

            tp53Min = geneLocationInformation.genesByName["TP53"].minLocus;
            tp53Max = geneLocationInformation.genesByName["TP53"].maxLocus;

            outputFile = ASETools.CreateStreamWriterWithRetry(@"\temp\tp53_vaf.txt");
            outputFile.WriteLine("case ID\tgermline VAF\tsomatic mutation count");

            var threading = new ASETools.ASVThreadingHelper<int>(cases, aseCorrection, (x, y) => true, handleOneCase, null, null, 100);

            Console.WriteLine("Processing " + cases.Count() + " cases, 1 dot/100: ");
            ASETools.PrintNumberBar(cases.Count() / 100);

            threading.run();

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(timer));

            outputFile.WriteLine("**done**");
            outputFile.Close();
        }  // Main

        static void handleOneCase(ASETools.Case case_, int state, List<ASETools.AnnotatedVariant> variantsForThisCase, ASETools.ASECorrection aseCorrection, Dictionary<bool, List<ASETools.CopyNumberVariation>> copyNumber)
        {
            var tp53mutations = variantsForThisCase.Where(x => x.contig == "chr17" && x.locus >= tp53Min && x.locus <= tp53Max).ToList();

            var somatic = tp53mutations.Where(x => x.somaticMutation && !x.isSilent()).ToList();

            var germline = tp53mutations.Where(x => !x.somaticMutation && x.IsASECandidate(true, copyNumber, configuration, perGeneASEMap, geneMap)).ToList();

            if (germline.Count() != 0)
            {
                lock(outputFile)
                {
                    outputFile.WriteLine(case_.case_id + "\t" + germline[0].GetAltAlleleFraction(true) + "\t" + somatic.Count());
                }
            }
        } // handleOneCase

    } // Program


}
