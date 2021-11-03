using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;

namespace ConsolodateCasePairedness
{
    class Program
    {
        static ASETools.CommonData commonData;
        static List<ASETools.CasePairedness> allPairedness = new List<ASETools.CasePairedness>();

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);

            if (null == commonData)
            {
                return;
            } // failed to load common data


            if (commonData.listOfCases.Any(_ => _.case_pairedness_filename == ""))
            {
                Console.WriteLine("Not all cases have pairedness computed.");
                return;
            }

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "cases", commonData.listOfCases.Count(), out nPerDot);
            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(commonData.listOfCases, HandleOneCase, null, null, nPerDot);
            threading.run();

            var outputFilename = commonData.configuration.finalResultsDirectory + ASETools.ConsolodatedCasePairednessFilename;
            ASETools.CasePairedness.WriteToFile(outputFilename, allPairedness);

            Console.WriteLine("Consolodated " + commonData.listOfCases.Count() + " cases in " + ASETools.ElapsedTimeInSeconds(commonData.timer));

        } // Main

        static void HandleOneCase(ASETools.Case case_, int state)
        {
            var pairedness = ASETools.CasePairedness.ReadFromFile(case_.case_pairedness_filename);

            lock (allPairedness)
            {
                allPairedness.Add(pairedness);
            } // lock
        } // HandleOneCase
    } // Program
} // namespace
