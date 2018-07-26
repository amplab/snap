using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.Diagnostics;

namespace ASEScatter
{
    class Program
    {
        static ASETools.Configuration configuration;
        static ASETools.GeneMap geneMap;
        static Dictionary<string, ASETools.ASEMapPerGeneLine> perGeneASEMap;
        static StreamWriter outputFile;
        static List<ASETools.GeneScatterGraphLine> p53ScatterGraphLines;


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

            p53ScatterGraphLines = ASETools.GeneScatterGraphLine.LoadAllGeneScatterGraphLines(configuration.geneScatterGraphsDirectory, true, "TP53").Where(x => x.Variant_Classification != "Silent").ToList();
            if (null == p53ScatterGraphLines || p53ScatterGraphLines.Count() == 0)
            {
                Console.WriteLine("Unable to load p53 scatter graph lines.");
            }

            outputFile = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory + ASETools.Configuration.GlobalScatterGraphFilename);
            if (null == outputFile)
            {
                Console.WriteLine("Unable to open output file");
                return;
            }

            outputFile.WriteLine("case id\tcancer type\tn Germline Variant ASE Candidate Sites\tmutation count\tmean ASE\tTP53 classification");


            var threading = new ASETools.ASVThreadingHelper<int>(cases, aseCorrection, (x, y) => true, handleOneCase, null, null, 100);

            Console.WriteLine("Processing " + cases.Count() + " cases, 1 dot/100: ");
            ASETools.PrintNumberBar(cases.Count() / 100);

            threading.run();

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(timer));
        } // Main

        static void handleOneCase(ASETools.Case case_, int state, List<ASETools.AnnotatedVariant> variantsForThisCase, ASETools.ASECorrection aseCorrection, Dictionary<bool, List<ASETools.CopyNumberVariation>> copyNumber)
        {
            int nMutations = variantsForThisCase.Where(x => x.somaticMutation).Count();
            var aseCandidates = variantsForThisCase.Where(x => !x.somaticMutation && x.IsASECandidate(true, copyNumber, configuration, perGeneASEMap, geneMap)).ToList();
            int nASECandidates = aseCandidates.Count();

            if (nASECandidates == 0)
            {
                return;
            }

            double meanASE = aseCandidates.Select(x => x.GetAlleleSpecificExpression(true, aseCorrection)).Average();
            string p53Classification = "";

            var p53Mutations = p53ScatterGraphLines.Where(x => x.case_id == case_.case_id).ToList();
            if (p53Mutations.Count() == 0)
            {
                p53Classification = "0";
            } else if (p53Mutations.Count() > 1)
            {
                p53Classification = "2";
            } else
            {
                var line = p53Mutations[0];
                if (line.tumorDNAReadCounts.usefulReads() < configuration.minDNAReadCoverage)
                {
                    p53Classification = "Insufficient Data";
                } else if (line.tumorDNAReadCounts.AltFraction() < 0.4)
                {
                    p53Classification = "Subclone";
                } else if (line.tumorDNAReadCounts.AltFraction() > 0.6)
                {
                    p53Classification = "Loss of Het";
                } else if (line.tumorRNAReadCounts.usefulReads() < configuration.minRNAReadCoverage)
                {
                    p53Classification = "Insufficient Data";
                } else if (ASETools.NonsenseMediatedDecayCausingVariantClassifications.Contains(line.Variant_Classification))
                {
                    p53Classification = "Nonsense Mediated Decay";
                } else if (line.tumorRNAReadCounts.AltFraction() > 0.6)
                {
                    p53Classification = ">60% Mutant";
                } else if (line.tumorRNAReadCounts.AltFraction() < 0.4)
                {
                    p53Classification = "<40% Mutant";
                } else
                {
                    p53Classification = "Even Expression";
                }
            }

            lock (outputFile)
            {
                outputFile.WriteLine(case_.case_id + "\t" + case_.disease() + "\t" + nASECandidates + "\t" + nMutations + "\t"  +  meanASE + "\t" + p53Classification);
            }
        } // handleOneCase
    }
}
