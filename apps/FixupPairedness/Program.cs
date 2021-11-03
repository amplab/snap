using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;

// The initial cut of the CasePairedness class left out the **done** line.  Fix that by reading them and writing them back out with the fixed-up version of ASELib.dll

namespace FixupPairedness
{
    class Program
    {

        static void Main(string[] args)
        {
            var commonData = ASETools.CommonData.LoadCommonData(args);

            if (commonData == null)
            {
                return;
            }

            var casesToRun = commonData.listOfCases.Where(_ => _.case_pairedness_filename != "").ToList();

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "cases", casesToRun.Count(), out nPerDot);
            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(casesToRun, HandleOneCase, null, null, nPerDot);
            threading.run();

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(commonData.timer));
        } // Main

        static void HandleOneCase(ASETools.Case case_, int state)
        {
            var pairedness = ASETools.CasePairedness.ReadFromFile(case_.case_pairedness_filename);

            ASETools.CasePairedness.WriteToFile(case_.case_pairedness_filename, pairedness);    // Re-write it with the **done** this time.
        }
    } // Program
} // namespace
