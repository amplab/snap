using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using ASELib;
using System.Threading;

namespace ComputePerCaseASE
{
    class Program
    {
        static ASETools.GeneMap geneMap;
        static ASETools.Configuration configuration;
        static Dictionary<bool, Dictionary<string, ASETools.ASEMapPerGeneLine>> perGeneASEMap;

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            configuration = ASETools.Configuration.loadFromFile(args);

            if (null == configuration)
            {
                Console.WriteLine("Giving up because we were unable to load configuration.");
                return;
            }

            var geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(configuration.geneLocationInformationFilename));
            geneMap = new ASETools.GeneMap(geneLocationInformation.genesByName);
            perGeneASEMap = ASETools.ASEMapPerGeneLine.ReadFromFileToDictionary(configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

            if (null == cases)
            {
                Console.WriteLine("Unable to load cases.");
            }

            var casesToProcess = cases.Where(x => x.Value.per_case_ase_filename == "" && x.Value.annotated_selected_variants_filename != "" && x.Value.tumor_copy_number_filename != "").Select(x => x.Value).ToList();
            int nCases = casesToProcess.Count();
            Console.Write("Processing " + nCases + " cases, 1 dot/100 cases: ");

            var threads = new List<Thread>();
            for (int i = 0; i < Environment.ProcessorCount; i++)
            {
                threads.Add(new Thread(() => WorkerThread(casesToProcess)));
            }

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            Console.WriteLine();
            Console.WriteLine("Processed " + nCases + " in " + ASETools.ElapsedTimeInSeconds(timer));
        }

        static int nCasesProcessed = 0;

        class Measurements
        {
            public int nNormal = 0, nTumor = 0;
            public double normal = 0, tumor = 0;

            public string getASE(bool forTumor)
            {
                if (forTumor)
                {
                    if (nTumor == 0)
                    {
                        return "*";
                    }
                    return Convert.ToString(tumor / nTumor);
                } else
                {
                    if (nNormal == 0)
                    {
                        return "*";
                    }
                    return Convert.ToString(normal / nNormal);
                }
            }
        }

        static void WorkerThread(List<ASETools.Case> casesToProcess)
        {
            while (true)
            {
                ASETools.Case case_;

                lock (casesToProcess)
                {
                    if (casesToProcess.Count() == 0)
                    {
                        return;
                    }

                    case_ = casesToProcess[0];
                    casesToProcess.RemoveAt(0);

                    nCasesProcessed++;
                    if (nCasesProcessed % 100 == 0)
                    {
                        Console.Write('.');
                    }
                } // lock

                var annotatedSelectedVariants = ASETools.AnnotatedVariant.readFile(case_.annotated_selected_variants_filename);
                var tumorCopyNumberVariation = ASETools.CopyNumberVariation.ReadFile(case_.tumor_copy_number_filename, case_.tumor_copy_number_file_id).Where(r => Math.Abs(r.Segment_Mean) > 1.0).ToList();
                List<ASETools.CopyNumberVariation> normalCopyNumberVariation = null;
                if (case_.normal_copy_number_filename != "")
                {
                    normalCopyNumberVariation = ASETools.CopyNumberVariation.ReadFile(case_.normal_copy_number_filename, case_.normal_copy_number_file_id).Where(r => Math.Abs(r.Segment_Mean) > 1.0).ToList();
                }

                if (annotatedSelectedVariants == null)
                {
                    Console.WriteLine("Unable to read annotated selected variants from " + case_.annotated_selected_variants_filename);
                    continue;
                }


                var overall = new Measurements();
                var perChromosome = new Measurements[ASETools.nHumanNuclearChromosomes]; 
                for (int i = 0; i < ASETools.nHumanNuclearChromosomes; i++)
                {
                    perChromosome[i] = new Measurements();
                }

                foreach (var variant in annotatedSelectedVariants)
                {
                    int index = ASETools.ChromosomeNameToIndex(variant.contig);

                    if (!variant.somaticMutation && variant.IsASECandidate(false, normalCopyNumberVariation, configuration, perGeneASEMap, geneMap))
                    {
                        overall.nNormal++;
                        overall.normal += variant.GetNormalAlleleSpecificExpression();
                        if (index != -1) {
                            perChromosome[index].nNormal++;
                            perChromosome[index].normal += variant.GetNormalAlleleSpecificExpression();
                        }
                    }

                    if (!variant.somaticMutation && variant.IsASECandidate(true, tumorCopyNumberVariation, configuration, perGeneASEMap, geneMap))
                    {
                        overall.nTumor++;
                        overall.tumor += variant.GetTumorAlleleSpecificExpression();
                        if (index != -1)
                        {
                            perChromosome[index].nTumor++;
                            perChromosome[index].tumor += variant.GetTumorAlleleSpecificExpression();
                        }
                    }
                } // foreach variant

                var outputFilename = ASETools.GetDirectoryFromPathname(case_.annotated_selected_variants_filename) + @"\" + case_.case_id + ASETools.perCaseASEExtension;
                var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
                if (outputFile == null)
                {
                    Console.WriteLine("Unable to open output file " + outputFilename + ".  Skipping.");
                    continue;
                }

                outputFile.Write("Normal ASE\tTumor ASE");
                for (int i = 0; i < ASETools.nHumanNuclearChromosomes; i++)
                {
                    outputFile.Write("\tNormal " + ASETools.ChromosomeIndexToName(i) + " ASE\tTumor " + ASETools.ChromosomeIndexToName(i) + "ASE");
                }
                outputFile.WriteLine();

                outputFile.Write(overall.getASE(false) + "\t" + overall.getASE(true));
                for (int i = 0; i < ASETools.nHumanNuclearChromosomes; i++)
                {
                    outputFile.Write("\t" + perChromosome[i].getASE(false) + "\t" + perChromosome[i].getASE(true));
                }
                outputFile.WriteLine();
                outputFile.WriteLine("**done**");
                outputFile.Close();

            } // while true
        } // WorkerThread
    }
}
