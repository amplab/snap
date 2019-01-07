using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.Diagnostics;


namespace DistanceBetweenMutations
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


            var outputFile = ASETools.CreateStreamWriterWithRetry(commonData.configuration.finalResultsDirectory + ASETools.DistanceBetweenMutationsFilename);

            if (null == outputFile)
            {
                Console.WriteLine("Unable to open output file.");
                return;
            }

            var geneScatterGraphLines = ASETools.GeneScatterGraphLine.LoadAllGeneScatterGraphLines(commonData.configuration.geneScatterGraphsDirectory, false, "*");

            var geneScatterGraphLinesByGene = geneScatterGraphLines.GroupByToDict(_ => _.Hugo_Symbol);

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "genes, ", geneScatterGraphLinesByGene.Count(), out nPerDot);

            int nProcessed = 0;
            foreach (var linesForGene in geneScatterGraphLinesByGene.Select(_ => _.Value))
            {
                var histogram = new ASETools.PreBucketedHistogram(0, 2000, 50);

                var linesPerPatient = linesForGene.GroupByToDict(_ => _.case_id);

                foreach (var perPatientPerGene in linesPerPatient.Select(_ => _.Value).Where(_ => _.Count() > 1))
                {
                    perPatientPerGene.Sort((a, b) => a.Start_Position.CompareTo(b.Start_Position));

                    for (int i = 0; i < perPatientPerGene.Count() - 1; i++)
                    {
                        histogram.addValue(perPatientPerGene[i + 1].Start_Position - perPatientPerGene[i].Start_Position);
                    } // each mutation save the last
                } // each patient with multiple mutations in this gene

                if (histogram.count() >= 10)
                {
                    outputFile.WriteLine(linesForGene[0].Hugo_Symbol);
                    outputFile.WriteLine(ASETools.HistogramResultLine.Header());
                    histogram.ComputeHistogram().ToList().ForEach(_ => outputFile.WriteLine(_));
                    outputFile.WriteLine();
                }

                nProcessed++;
                if (nProcessed % nPerDot == 0)
                {
                    Console.Write(".");
                }
            } // each gene

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine();
            Console.WriteLine("Run time " + ASETools.ElapsedTimeInSeconds(commonData.timer));
        } // Main
    }
}
