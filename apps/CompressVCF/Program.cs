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

        static ASETools.CommonData commonData;

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);

            if (null == commonData)
            {
                return;
            }

            if (commonData.configuration.commandLineArgs.Count() <1 )
            {
                Console.WriteLine("Usage: CompressVCF {-configration configuration} [caseId]");
                return;
            }

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "VCFs", commonData.configuration.commandLineArgs.Count(), out nPerDot);
            var threading = new ASETools.WorkerThreadHelper<string, int>(commonData.configuration.commandLineArgs.ToList(), HandleOneCase, null, null, nPerDot);
            threading.run();

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(commonData.timer));
        } // Main

        static void HandleOneCase(string caseId, int state)
        {
            if (!commonData.cases.ContainsKey(caseId))
            {
                Console.WriteLine(caseId + " is not a valid case ID.");
                return;
            }

            var case_ = commonData.cases[caseId];
            if (case_.vcf_filename == "")
            {
                Console.WriteLine("Case " + caseId + " does not have a VCF file to compress.");
                return;
            }

            var commandLine = "/c " + commonData.configuration.binariesDirectory + "cat " + case_.vcf_filename + " | " + commonData.configuration.binariesDirectory + "gzip.exe -9" + " > " + case_.vcf_filename + ".gz";
            ASETools.RunAndWaitForProcess(commonData.configuration.binariesDirectory + "cmd.exe", commandLine);

        } // HandleOneCase
    }
}
