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
        static void processOne(StreamWriter outputFile, string header, List<ASETools.ClinicalSummaryLine> clinicalLines)
        {
            int n;
            var km = ASETools.KaplanMeier(clinicalLines, out n);

            outputFile.WriteLine(header + " (n= " + n + ")");
            outputFile.WriteLine("Days\tFraction Alive");
            km.ForEach(x => outputFile.WriteLine(x.daysSinceInclusion + "\t" + x.fractionAlive));
        }
        static void Main(string[] args)
        {
            var clinicalLines = ASETools.ClinicalSummaryLine.readFile(@"\sequence\BeatAML\Wave 2 FINAL data lock (Jun 2017)\Functional & clinical\clinical-summary-20170703_143804.txt");
            var mappings = ASETools.BeatAMLSampleMapping.readFromFile(@"\sequence\beatAML\Wave 2 FINAL data lock (Jun 2017)\BeatAML_sample_mapping_file_8_17_2017.txt");
            var amlOnly = clinicalLines.Where(x => x.dxAtInclusion.ToUpper().Contains("ACUTE MYELOID LEUKAEMIA")).ToList();
            var myelodisplastic = clinicalLines.Where(x => x.dxAtInclusion.ToUpper().Contains("MYELODYSPLASTIC")).ToList();

            int nWithDataAndTP53 = 0;
            foreach (var caseEntry in mappings)
            {

                var patientId = caseEntry.Key;
                if (amlOnly.Where(x => x.patientId == patientId && x.p53Mutant()).Count() == 0)
                {
                    continue;
                }

                var sampleMappings = caseEntry.Value;   // A list of sample mapping lines.
                if (sampleMappings.Where(x => x.Normal_BAM != "").Count() > 0 && sampleMappings.Where(x => x.AML_BAM != "").Count() > 0 && sampleMappings.Where(x => x.RNA_BAM != "").Count() > 0)
                {
                    nWithDataAndTP53++;

                    Console.WriteLine(sampleMappings.Where(x => x.Normal_BAM != "").ToList()[0].Normal_BAM);
                    Console.WriteLine(sampleMappings.Where(x => x.AML_BAM != "").ToList()[0].AML_BAM);
                    Console.WriteLine(sampleMappings.Where(x => x.RNA_BAM != "").ToList()[0].RNA_BAM);
                }
            }

            Console.WriteLine("There are " + nWithDataAndTP53 + " AML cases with TP53 mutations and all three BAM files.");

            var outputFile = ASETools.CreateStreamWriterWithRetry(@"\temp\BeatAML.txt");

            processOne(outputFile, "AML", amlOnly);
            processOne(outputFile, "AML w/TP53 mutation", amlOnly.Where(x => x.p53Mutant()).ToList());
            processOne(outputFile, "AML w/TP53 mutation and hypomethylating", amlOnly.Where(x => x.p53Mutant() &&
                        (x.priorTreatmentRegimens.ToLower().Contains("azacitidine") || x.priorTreatmentRegimens.ToLower().Contains("decitibine") ||
                        x.currentRegimen.ToLower().Contains("azacitidine") || x.currentRegimen.ToLower().Contains("decitibine"))).ToList());
            processOne(outputFile, "AML w/TP53 mutation and no hypomethylating", amlOnly.Where(x => x.p53Mutant() &&
            !(x.priorTreatmentRegimens.ToLower().Contains("azacitidine") || x.priorTreatmentRegimens.ToLower().Contains("decitibine") ||
            x.currentRegimen.ToLower().Contains("azacitidine") || x.currentRegimen.ToLower().Contains("decitibine"))).ToList());
            processOne(outputFile, "Myelodysplastic (all)", myelodisplastic);

            outputFile.Close();
        }
    }
}
