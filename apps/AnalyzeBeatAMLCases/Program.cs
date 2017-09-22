using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;

namespace AnalyzeBeatAMLCases
{
    class Program
    {
        static void Main(string[] args)
        {
            var clinicalLines = ASETools.ClinicalSummaryLine.readFile(@"\sequence\BeatAML\Wave 2 FINAL data lock (Jun 2017)\Functional & clinical\clinical-summary-20170703_143804.txt");

            var km = ASETools.KaplanMeier(clinicalLines);

            var outputFile = ASETools.CreateStreamWriterWithRetry(@"\temp\BeatAML.txt");
            outputFile.WriteLine("Days\tFraction Alive");
            km.ForEach(x => outputFile.WriteLine(x.daysSinceInclusion + "\t" + x.fractionAlive));

            outputFile.Close();
        }
    }
}
