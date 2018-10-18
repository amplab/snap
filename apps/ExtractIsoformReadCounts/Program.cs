using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;

namespace ExtractIsoformReadCounts
{
    class Program
    {
        static ASETools.CommonData commonData;
        static ASETools.IsoformMap isoformMap;

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);

            if (null == commonData)
            {
                return;
            }

            if (commonData.configuration.commandLineArgs.Count() == 0)
            {
                Console.WriteLine("usage: ExtractIsoformReads <caseId>");
                return;
            }

            isoformMap = new ASETools.IsoformMap(commonData.geneLocationInformation.genesByName);

            var casesToProcess = new List<ASETools.Case>();
            foreach (var caseId in commonData.configuration.commandLineArgs)
            {
                if (!commonData.cases.ContainsKey(caseId))
                {
                    Console.WriteLine(caseId + " does not appear to be a case ID");
                    return;
                }

                casesToProcess.Add(commonData.cases[caseId]);
            }


            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "cases", casesToProcess.Count(), out nPerDot);
            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(casesToProcess, HandleOneCase, null, null, nPerDot);
            threading.run();

            Console.WriteLine("Processed " + casesToProcess.Count() + " cases in " + ASETools.ElapsedTimeInSeconds(commonData.timer));
        } // Main

        static void HandleOneCase(ASETools.Case case_, int state)
        {
            var readCounts = new Dictionary<string, Dictionary<bool, int>>();
            int defaultNormal = case_.normal_rna_allcount_filename == "" ? -1 : 0;  // -1 when we don't have data.

            foreach (var isoform in isoformMap.getAllIsoforms())
            {
                readCounts.Add(isoform.ucscId, new Dictionary<bool, int>());
                readCounts[isoform.ucscId].Add(true, 0);
                readCounts[isoform.ucscId].Add(false, defaultNormal);
            }

            foreach (var tumor in ASETools.BothBools)
            {
                var allcountFilename = tumor ? case_.tumor_rna_allcount_filename : case_.normal_rna_allcount_filename;
                if (allcountFilename == "")
                {
                    continue;   // Most cases don't have normal RNA.  Just skip them.
                }

                var allcountReader = new ASETools.AllcountReader(case_.tumor_rna_allcount_filename);

                if (!allcountReader.openFile())
                {
                    Console.WriteLine("Unable to open allcount file " + case_.tumor_rna_allcount_filename);
                    return;
                }

                allcountReader.ReadAllcountFile((contigName, location, currentMappedReadCount) => processBase(tumor, readCounts, contigName, location, currentMappedReadCount));
            } // tumor

            var allReadCounts = readCounts.Select(_ => new ASETools.IsoformReadCounts(_.Key, _.Value[true], _.Value[false])).ToList();

            ASETools.IsoformReadCounts.WriteToFile(ASETools.GetDirectoryFromPathname(case_.tumor_rna_allcount_filename) + @"\" + case_.case_id + ASETools.isoformReadCountsExtension, allReadCounts);
        } // HandleOneCase

        static void processBase(bool tumor,  Dictionary<string, Dictionary<bool, int>> readCounts, string contigName, int location, int currentMappedReadCount)
        {
            foreach (var isoform in isoformMap.getIsoformsMappedTo(contigName, location))
            {
                readCounts[isoform.ucscId][tumor] += currentMappedReadCount;
            } // foreach isoform
        } // processBase
    }
}
