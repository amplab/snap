using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.Security.AccessControl;

namespace SummarizeVenn
{
    class Program
    {
        class ConcordanceState
        {
            public Dictionary<ASETools.VariantType, Dictionary<ASETools.AlignerSet, ASETools.PreBucketedHistogram>> weighted = new Dictionary<ASETools.VariantType, Dictionary<ASETools.AlignerSet, ASETools.PreBucketedHistogram>>();
            public Dictionary<ASETools.VariantType, Dictionary<ASETools.AlignerSet, ASETools.PreBucketedHistogram>> raw = new Dictionary<ASETools.VariantType, Dictionary<ASETools.AlignerSet, ASETools.PreBucketedHistogram>>();

            public ConcordanceState()
            {
                foreach (var variantType in ASETools.EnumUtil.GetValues<ASETools.VariantType>())
                {
                    weighted.Add(variantType, new Dictionary<ASETools.AlignerSet, ASETools.PreBucketedHistogram>());
                    raw.Add(variantType, new Dictionary<ASETools.AlignerSet, ASETools.PreBucketedHistogram>());

                    foreach (var alignerSet in ASETools.allAlignerSets)
                    {
                        weighted[variantType].Add(alignerSet, new ASETools.PreBucketedHistogram(0, 1, 0.01));
                        raw[variantType].Add(alignerSet, new ASETools.PreBucketedHistogram(0, 1000000, 1000));
                    }
                }
            } // ctor

            public void merge(ConcordanceState peer)
            {
                foreach (var variantType in ASETools.EnumUtil.GetValues<ASETools.VariantType>())
                {
                    foreach (var alignerSet in ASETools.allAlignerSets)
                    {
                        weighted[variantType][alignerSet].merge(peer.weighted[variantType][alignerSet]);
                        raw[variantType][alignerSet].merge(peer.raw[variantType][alignerSet]);
                    }
                }
            } // merge

        } // ConcordanceState

        static ConcordanceState globalState = new ConcordanceState();

        static ASETools.Configuration configuration;
        static ASETools.VariantCaller variantCaller;
        static bool tumor;

        static void Main(string[] args)
        {
            configuration = ASETools.Configuration.loadFromFile(args);
            if (configuration == null)
            {
                return;
            }

            if (configuration.commandLineArgs.Count() != 2)
            {
                Console.WriteLine("Usage: SummarizeVenn VariantCaller TumorOrNormal");
                return;
            }

            if (!ASETools.tumorToString.ContainsValue(configuration.commandLineArgs[1]))
            {
                Console.WriteLine("Second parameter must be either " + ASETools.tumorToString[true] + " or " + ASETools.tumorToString[false]);
                return;
            }

            if (!ASETools.variantCallerName.ContainsValue(configuration.commandLineArgs[0]))
            {
                Console.Write("First parameter must be a variant caller name, one of:");
                ASETools.alignerName.Select(_ => _.Value).ToList().ForEach(_ => Console.Write(" " + _));
                Console.WriteLine();
                return;
            }

            tumor = configuration.commandLineArgs[1] == ASETools.tumorToString[true];
            variantCaller = ASETools.variantCallerName.Where(_ => _.Value == configuration.commandLineArgs[0]).ToList()[0].Key;

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);
            if (null == cases)
            {
                Console.WriteLine("Unable to load cases");
                return;
            }

            var listOfCases = cases.Select(_ => _.Value).ToList();

            if (listOfCases.Any(_ => _.perVariantCaller[variantCaller][tumor].venn_filename == ""))
            {
                Console.WriteLine("Some cases are missing Venn run.");
                //return;
            }

            var casesToRun = listOfCases.Where(_ => _.perVariantCaller[variantCaller][tumor].venn_filename != "").ToList();

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "cases", casesToRun.Count(), out nPerDot);
            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, ConcordanceState>(casesToRun, HandleOneCase, FinishUp, null, nPerDot);
            threading.run();

            var outputFile = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory + ASETools.VennFilename);
            if (null == outputFile)
            {
                return;
            }

            foreach (var variantType in ASETools.EnumUtil.GetValues<ASETools.VariantType>())
            {
                outputFile.WriteLine("Results for " + ASETools.variantTypeToName[variantType] + "s");
                outputFile.WriteLine("Normalized");
                var histograms = globalState.weighted[variantType].Select(_ => new KeyValuePair<string, ASETools.PreBucketedHistogram>(_.Key.ToString(), _.Value)).ToList();
                ASETools.PreBucketedHistogram.WriteBatchOfHistogramCDFs(outputFile, histograms);

                outputFile.WriteLine("Means: ");
                histograms.ForEach(_ => outputFile.WriteLine(_.Key + "\t" + _.Value.mean()));

                outputFile.WriteLine();

                outputFile.WriteLine("Results for " + ASETools.variantTypeToName[variantType] + "s");
                outputFile.WriteLine("Not normalized");
                histograms = globalState.raw[variantType].Select(_ => new KeyValuePair<string, ASETools.PreBucketedHistogram>(_.Key.ToString(), _.Value)).ToList();
                ASETools.PreBucketedHistogram.WriteBatchOfHistogramCDFs(outputFile, histograms);

                outputFile.WriteLine("Means: ");
                histograms.ForEach(_ => outputFile.WriteLine(_.Key + "\t" + _.Value.mean()));

                outputFile.WriteLine();
            } // variantType

            outputFile.WriteLine("**done**");

            outputFile.Close();
        } // Main

        static void HandleOneCase(ASETools.Case case_, ConcordanceState state)
        {
            var vennResults = ASETools.VennResults.LoadFromFile(case_.perVariantCaller[variantCaller][tumor].venn_filename);
           //
           // If it failed to load it will be null and it printed an error message.  Just let the next statement
           // throw an exception, we don't want to finish with missing data.
           //
           foreach (var variantType in ASETools.EnumUtil.GetValues<ASETools.VariantType>())
            {
                foreach (var alignerSet in ASETools.allAlignerSets)
                {
                    state.raw[variantType][alignerSet].addValue(vennResults.countBySet[variantType][alignerSet]);
                    state.weighted[variantType][alignerSet].addValue((double)vennResults.countBySet[variantType][alignerSet] / vennResults.nVariants[variantType]);
                } // alignerSet
            } // variantType
        } // HandleOneCase

        static void FinishUp(ConcordanceState state)
        {
            lock (globalState)
            {
                globalState.merge(state);
            } // lock
        } // FinishUp
    }
}
