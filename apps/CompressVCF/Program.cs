using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.Diagnostics;

namespace CompressVCF
{
    class Program
    {
        static ASETools.Configuration configuration;
        static Dictionary<string, ASETools.Case> cases;

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            configuration = ASETools.Configuration.loadFromFile(args);

            if (null == configuration)
            {
                return;
            }

            cases = ASETools.Case.LoadCases(configuration.casesFilePathname);
            if (null == cases)
            {
                return;
            }

            if (configuration.commandLineArgs.Count() < 1 )
            {
                Console.WriteLine("Usage: CompressVCF {-configration configuration} [caseId]");
                return;
            }

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "VCFs", configuration.commandLineArgs.Count(), out nPerDot);
            var threading = new ASETools.WorkerThreadHelper<string, int>(configuration.commandLineArgs.ToList(), HandleOneCase, null, null, nPerDot);
            threading.run();

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(timer));
        } // Main

        static void HandleOneCase(string caseId, int state)
        {
            if (!cases.ContainsKey(caseId))
            {
                Console.WriteLine(caseId + " is not a valid case ID.");
                return;
            }

            var case_ = cases[caseId];
            if (case_.vcf_filename == "")
            {
                Console.WriteLine("Case " + caseId + " does not have a VCF file to compress.");
                return;
            }

            var commandLine = "/c " + configuration.binariesDirectory + "cat " + case_.vcf_filename + " | " + configuration.binariesDirectory + "gzip.exe -9" + " > " + case_.vcf_filename + ".gz";
            ASETools.RunAndWaitForProcess(configuration.binariesDirectory + "cmd.exe", commandLine);

        } // HandleOneCase
    }
}
