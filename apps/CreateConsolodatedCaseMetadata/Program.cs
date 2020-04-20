using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.CodeDom;

namespace CreateConsolodatedCaseMetadata
{
    class Program
    {
        static List<ASETools.CaseMetadata> allMetadata = new List<ASETools.CaseMetadata>();
        static ASETools.CommonData commonData = null;
        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);
            if (null == commonData)
            {
                return;
            }

            if (commonData.listOfCases.Any(_ => _.case_metadata_filename == ""))
            {
                Console.WriteLine("Not all cases have metadata computed.");
                return;
            }

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "cases", commonData.listOfCases.Count(), out nPerDot);
            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(commonData.listOfCases, HandleOneCase, null, null, nPerDot);
            threading.run();

            var outputFilename = commonData.configuration.finalResultsDirectory + ASETools.ConsolodatedCaseMetadataFilename;
            ASETools.CaseMetadata.WriteToFile(outputFilename, allMetadata);

            Console.WriteLine("Consolodated " + commonData.listOfCases.Count() + " cases in " + ASETools.ElapsedTimeInSeconds(commonData.timer));

        } // Main

        static void HandleOneCase(ASETools.Case case_, int state)
        {
            var caseMetadata = ASETools.CaseMetadata.ReadFromFile(case_.case_metadata_filename);
            if (null == caseMetadata)
            {
                throw new Exception("Unable to read case metadata for case " + case_.case_id + " from " + case_.case_metadata_filename);
            }

            lock (allMetadata)
            {
                allMetadata.Add(caseMetadata);
            }
        } // HandleOneCase
    }
}
