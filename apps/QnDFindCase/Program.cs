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

        static int cdkn2aMin;
        static int cdkn2aMax;
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

            cdkn2aMin = geneLocationInformation.genesByName["CDKN2A"].minLocus;
            cdkn2aMax = geneLocationInformation.genesByName["CDKN2A"].maxLocus;

            outputFile = ASETools.CreateStreamWriterWithRetry(@"\temp\cdkn2a_vaf.txt");
            outputFile.WriteLine("case ID\tgermline VAF (only ASE candidates)\tsomatic mutation count\tsomatic VAF (mean of all ASE candidate mutations)");

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
            var cdkn2amutations = variantsForThisCase.Where(x => x.contig == "chr9" && x.locus >= cdkn2aMin && x.locus <= cdkn2aMax).ToList();

            var somatic = cdkn2amutations.Where(x => x.somaticMutation && !x.isSilent()).ToList();

            var germline = cdkn2amutations.Where(x => !x.somaticMutation && x.IsASECandidate(true, copyNumber, configuration, perGeneASEMap, geneMap)).ToList();

            if (somatic.Count() == 0 && germline.Count() == 0)
            {
                return;
            }

            lock (outputFile)
            {
                outputFile.Write(case_.case_id + "\t");
                if (germline.Count() > 0)
                {
                    outputFile.Write(germline[0].GetAltAlleleFraction(true));
                } else
                {
                    outputFile.Write("*");
                }

                outputFile.Write("\t" + somatic.Count() + "\t");

                if (somatic.Where(x => x.IsASECandidate(true, copyNumber, configuration, perGeneASEMap, geneMap)).Count() > 0)
                {
                    outputFile.WriteLine(somatic.Where(x => x.IsASECandidate(true, copyNumber, configuration, perGeneASEMap, geneMap)).Select(x => x.GetAltAlleleFraction(true)).Average());
                } else
                {
                    outputFile.WriteLine("*");
                }
            }
        } // handleOneCase

    } // Program


}
