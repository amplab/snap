using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;
using ASELib;

namespace AllSitesReadDepthDistribution
{
    class Program
    {

        static ASETools.Configuration configuration;
        static Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>> readDepth = new Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>>();  // tumor, DNA in that order

        const int maxReadDepthInHistogram = 2000;
        const int nCasesPerDot = 100;

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            configuration = ASETools.Configuration.loadFromFile(args);
            if (null == configuration)
            {
                Console.WriteLine("Unable to load configuration");
                return;
            }

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname).Select(x => x.Value).ToList();   
            var nCases = cases.Count();

            foreach (var tumor in ASETools.BothBools)
            {
                readDepth.Add(tumor, new Dictionary<bool, ASETools.PreBucketedHistogram>());
                foreach (var dna in ASETools.BothBools)
                {
                    readDepth[tumor].Add(dna, new ASETools.PreBucketedHistogram(0, maxReadDepthInHistogram, 1, "read depth"));
                } // dna
            } // tumor

            Console.Write("Processing " + nCases + " cases, " + nCasesPerDot + " cases/dot: ");

            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(cases, HandleOneCase, null, null, nCasesPerDot);
            threading.run();

            Console.WriteLine();

            var outputFile = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory + ASETools.AllSitesReadDepthFilename);

            foreach (var tumor in ASETools.BothBools)
            {
                foreach (var dna in ASETools.BothBools)
                {
                    outputFile.WriteLine("Read depth distribution at all sites tumor: " + tumor + " dna: " + dna);
                    outputFile.WriteLine(ASETools.HistogramResultLine.Header());
                    readDepth[tumor][dna].ComputeHistogram().ToList().ForEach(x => outputFile.WriteLine(x.ToString()));
                    outputFile.WriteLine();
                } // dna
            } // tumor

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine("Processed " + nCases + " cases in " + ASETools.ElapsedTimeInSeconds(timer));
        }

        static void HandleOneCase(ASETools.Case case_, int unused)
        {
            foreach(var tumor in ASETools.BothBools)
            {
                foreach (var dna in ASETools.BothBools)
                {
                    var inputFilename = case_.getAllcountFilename(tumor, dna);
                    if ("" == inputFilename) continue;  // Some normal RNA files don't exist.  Skip them.

                    var allcountReader = new ASETools.AllcountReader(inputFilename);
                    if (!allcountReader.openFile())
                    {
                        Console.WriteLine("Unable to open allcount file " + inputFilename);
                        continue;
                    }

                    var histogram = new ASETools.PreBucketedHistogram(0, maxReadDepthInHistogram, 1, "local read depth");

                    allcountReader.ReadAllcountFile((x, y, z) => processBase(histogram, x, y, z));

                    lock (readDepth[tumor][dna])
                    {
                        readDepth[tumor][dna].merge(histogram);
                    }
                } // dna
            } // tumor
        } // handle one case

        static void processBase(ASETools.PreBucketedHistogram histogram, string contigName, int location, int currentMappedReadCount)
        {
            histogram.addValue(currentMappedReadCount);
        }
    }
}
