using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;
using ASELib;

namespace ComputeExpressionDistribution
{
    class Program
    {

        static ASETools.ThrottledParallelQueue<ASETools.Case> [] resultQueues;

        static ASETools.CommonData commonData;

        class ReadDepth
        {
            public readonly string contig;
            public readonly int location;
            public readonly double fractionOfHQMappedBases;

            public ReadDepth(string contig_, int location_, double fractionOfHQMappedBases_)
            {
                contig = contig_;
                location = location_;
                fractionOfHQMappedBases = fractionOfHQMappedBases_;
            }
        }

        class CaseAndResultQueue
        {
            public readonly ASETools.Case case_;
            public readonly ASETools.ThrottledParallelQueue<ReadDepth> queue = new ASETools.ThrottledParallelQueue<ReadDepth>(1024,1);

            public CaseAndResultQueue(ASETools.Case case__)
            {
                case_ = case__
            }
        }

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);
            if (null == commonData)
            {
                return;
            }

            if (commonData.configuration.commandLineArgs.Count() != 1 || !commonData.cases.Any(x => x.Value.disease() == commonData.configuration.commandLineArgs[0]))
            {
                Console.WriteLine("usage: ComputeExpressionDistribution disease");
                return;
            }

            var disease = commonData.configuration.commandLineArgs[0];

            if (commonData.cases.Select(x => x.Value).Any(x => x.disease() == disease && x.tumor_rna_allcount_filename == ""))
            {
                Console.WriteLine("Some cases don't have a tumor RNA allcount file.");
                return;
            }

            var casesToProcess = commonData.cases.Select(x => x.Value).Where(x => x.disease() == disease).Select(x => new CaseAndResultQueue(x)).ToList();

            var resultQueues = new ASETools.ThrottledParallelQueue<ReadDepth>[casesToProcess.Count()];

            var threading = new ASETools.WorkerThreadHelper<CaseAndResultQueue, int>(casesToProcess, HandleOneCase, null, null);

            // One thread per case.  We read all of the values for a given locus of all of the threads at one time, write them to the file and
            // repeat, so we don't have enormous memory requirements.

            threading.start(casesToProcess.Count());

            var outputFile = ASETools.CreateStreamWriterWithRetry


 
        }
    }
}
