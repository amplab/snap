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


            var threading = new ASETools.ASVThreadingHelper<int>(cases, aseCorrection, (x, y) => true, handleOneCase, null, null, 100);

            Console.WriteLine("Processing " + cases.Count() + " cases, 1 dot/100: ");
            ASETools.PrintNumberBar(cases.Count() / 100);

            threading.run();

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(timer));
        }  // Main

        static void handleOneCase(ASETools.Case case_, int state, List<ASETools.AnnotatedVariant> variantsForThisCase, ASETools.ASECorrection aseCorrection, Dictionary<bool, List<ASETools.CopyNumberVariation>> copyNumber)
        {
            var p53Mutations = variantsForThisCase.Where(x => x.Hugo_symbol == "TP53" && !x.isSilent()).ToList();
            if (p53Mutations.Count() != 1)
            {
                return;
            }

            if (!p53Mutations[0].IsASECandidate(true, copyNumber, configuration, perGeneASEMap, geneMap))
            {
                return;
            }

            if (p53Mutations[0].isSilent())
            {
                return;
            }

            if (p53Mutations[0].GetAltAlleleFraction(true) < 0.75)
            {
                return;
            }

            if (p53Mutations[0].getReadCount(true, true).usefulReads() < 30 || 
                p53Mutations[0].getReadCount(true, false).usefulReads() < 30 || 
                p53Mutations[0].getReadCount(false, true).usefulReads() < 30) {
                return;
            }

            lock (configuration)
            {
                Console.WriteLine(case_.case_id + " VAF: " + p53Mutations[0].GetAltAlleleFraction(true) + " mutation type: " + p53Mutations[0].variantClassification + ", locus: " + p53Mutations[0].locus);
            }

        } // handleOneCase

    } // Program


}
