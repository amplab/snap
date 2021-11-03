using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;

namespace MakeCasePairing
{
    class Program
    {
        static ASETools.CommonData commonData;

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);
            if (null == commonData)
            {
                return;
            }

            var caseIdsToRun = commonData.configuration.commandLineArgs.ToList();
            var casesToRun = commonData.listOfCases.Where(_ => caseIdsToRun.Contains(_.case_id)).ToList();

            var casesWithoutEnoughData = casesToRun.Where(_ => _.normal_dna_filename == "" || _.tumor_dna_filename == "" || _.tumor_rna_filename == "" ||
                                                (_.normal_rna_file_id != "" && (_.normal_rna_filename == "" || _.normal_rna_reads_at_tentative_selected_variants_filename == ""))).ToList();

            if (casesWithoutEnoughData.Count() != 0)
            {
                Console.Write("Missing some input data for cases:");
                casesWithoutEnoughData.ForEach(_ => Console.Write(" " + _.case_id));
                Console.WriteLine();
                Console.WriteLine("Skipping them.");

                casesToRun = casesToRun.Where(_ => !casesWithoutEnoughData.Contains(_)).ToList();
            }

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "cases", casesToRun.Count(), out nPerDot);
            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(casesToRun, HandleOneCase, null, null, nPerDot);
            threading.run();


            Console.WriteLine(ASETools.ElapsedTimeInSeconds(commonData.timer));
        } // Main

        static void HandleOneCase(ASETools.Case case_, int state)
        {
            var pairedness = new Dictionary<bool, Dictionary<bool, bool>>();

            foreach (var tumor in ASETools.BothBools)
            {
                pairedness.Add(tumor, new Dictionary<bool, bool>());

                foreach (var dna in ASETools.BothBools)
                {
                    var bamFilename = case_.getBAMFilename(tumor, dna);
                    if (bamFilename != "")
                    {
                        var flagstatOutput = ASETools.RunProcessAndGetOutput(commonData.configuration.binariesDirectory + "samtools.exe", "flagstat " + bamFilename);

                        if (flagstatOutput.Count() < 4 || !flagstatOutput[0].Contains("in total") || !flagstatOutput[3].Contains("paired in sequencing") || !flagstatOutput[0].Contains(" ") || !flagstatOutput[3].Contains(" "))
                        {
                            Console.WriteLine("Too few lines of output or bad output format for flagstat for BAM file " + bamFilename);
                            flagstatOutput.ForEach(_ => Console.WriteLine(_));
                            return;
                        } // flagstat has enough lines

                        var total = Convert.ToInt64(flagstatOutput[0].Substring(0, flagstatOutput[0].IndexOf(' ')));
                        if (total < 10000)
                        {
                            Console.WriteLine(bamFilename + " has too few total reads: " + total);
                            flagstatOutput.ForEach(_ => Console.WriteLine(_));
                            return;
                        }

                        var pairedInSequencing = Convert.ToInt64(flagstatOutput[3].Substring(0, flagstatOutput[3].IndexOf(' ')));
                        if (pairedInSequencing < 0 || pairedInSequencing > total)
                        {
                            Console.WriteLine(bamFilename + " has bogus paired in sequencing " + pairedInSequencing);
                            flagstatOutput.ForEach(_ => Console.WriteLine(_));
                            return;
                        }

                        pairedness[tumor].Add(dna, pairedInSequencing * 100 >= total);  // At least 1% paired means it's paired.
                    } else
                    {
                        // Sometimes there's no normal RNA.  We can't get here in other cases (and if by some reason we do, we'll just crash).
                        pairedness[false].Add(false, false);
                    } // if there's a BAM file
                } // foreach DNA
            } // foreach tumor

            var casePairedness = new ASETools.CasePairedness(case_.case_id, pairedness);

            var outputFilename = ASETools.ShareFromPathname(case_.normal_dna_filename) + @"\gdc\derived_files\" + case_.case_id + @"\" + case_.case_id + ASETools.casePairednessExtension;
            ASETools.CasePairedness.WriteToFile(outputFilename, casePairedness);
        } // HandleOneCase
    } // Program
} // namespace MakeCasePairing
