using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;

namespace VCFQualityDistribution
{
    class Program
    {
        static ASETools.CommonData commonData = null;
        static ASETools.Case.ColumnGetter vcfGetter;
        static PerThreadState globalState = new PerThreadState();
        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);
            if (commonData == null)
            {
                return;
            }

            if (commonData.configuration.commandLineArgs.Count() != 1)
            {
                Console.WriteLine("usage: VCFQualityDistribution bwa|bowtie|snap");
                return;
            }

            switch (commonData.configuration.commandLineArgs[0].ToLower())
            {
                case "bwa":
                    vcfGetter = _ => _.BWA_realigned_normal_feeebayes_vcf_filename;
                    break;

                case "bowtie":
                    vcfGetter = _ => _.bowtie_realigned_normal_freebayes_vcf_filename;
                    break;

                case "snap":
                    vcfGetter = _ => _.snap_realigned_normal_freebayes_vcf_filename;
                    break;

                default:
                    Console.WriteLine("Aligner must be one of bwa, bowtie or snap.");
                    return;
            } // switch;

            if (commonData.listOfCases.Any(_ => vcfGetter(_) == ""))
            {
                Console.WriteLine("Missing one or more input VCFs");
                // off for tesing return;
            }

            var casesToRun = commonData.listOfCases.Where(_ => vcfGetter(_) != "").ToList();
            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "cases", casesToRun.Count(), out nPerDot);
            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, PerThreadState>(casesToRun, HandleOneCase, FinishUp, null, nPerDot);
            threading.run(60);
            Console.WriteLine(ASETools.ElapsedTimeInSeconds(commonData.timer));

            var outputFile = ASETools.CreateStreamWriterWithRetry(commonData.configuration.commandLineArgs[0] + "_vcf_quality_distribution.txt");
            outputFile.WriteLine("Max " + globalState.histogram.max() + " mean " + globalState.histogram.mean());
            globalState.histogram.WriteHistogram(outputFile);
            outputFile.Close();

        } // Main

        class PerThreadState
        {
            public ASETools.PreBucketedHistogram histogram = new ASETools.PreBucketedHistogram(0, 300, 1);

            public void merge(PerThreadState peer)
            {
                histogram.merge(peer.histogram);
            }

        } // PerThreadState

        static void HandleOneCase(ASETools.Case case_, PerThreadState state)
        {
            var inputFile = ASETools.CreateStreamReaderWithRetry(vcfGetter(case_));
            if (inputFile == null)
            {
                throw new Exception("Unable to open " + vcfGetter(case_) + " for case id " + case_.case_id);
            }

            string inputLine;
            while (null != (inputLine = inputFile.ReadLine()))
            {
                if (inputLine.Length == 0 || inputLine[0] == '#')
                {
                    continue;
                }

                var fields = inputLine.Split('\t');
                if (fields.Count() < 6)
                {
                    throw new Exception("Malformed VCF: too few fields.  File " + vcfGetter(case_) + " line: " + inputLine);
                }

                state.histogram.addValue(Convert.ToDouble(fields[5]));
            }// while

            inputFile.Close();
        }// HandleOneCase

        static void FinishUp(PerThreadState state)
        {
            lock (globalState)
            {
                globalState.merge(state);
            }
        } // FinishUp
    } // Program
}
