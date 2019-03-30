using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.Diagnostics;

namespace GenerateTentativeASVsWithoutCNVs
{
    class Program
    {

        static ASETools.CommonData commonData;

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);

            if (commonData == null)
            {
                return;
            }

            if (commonData.configuration.commandLineArgs.Count() < 1 || commonData.configuration.commandLineArgs.Any(_ => !commonData.cases.ContainsKey(_)))
            {
                Console.WriteLine("usage: GenerateTentativeASVsWithoutCNVs <one or more case IDs>");
                return;
            }

            var casesToProcess = commonData.configuration.commandLineArgs.Select(_ => commonData.cases[_]).
                Where(_ => _.tentative_annotated_selected_variants_filename != "" && (_.tumor_copy_number_file_id == "" || _.tumor_copy_number_filename != "")).ToList();
            int nPerDot;

            ASETools.PrintMessageAndNumberBar("Processing", "cases", casesToProcess.Count(), out nPerDot);
            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(casesToProcess, HandleOneCase, null, null, nPerDot);
            threading.run();
            Console.WriteLine(ASETools.ElapsedTimeInSeconds(commonData.timer));
        }

        static void HandleOneCase(ASETools.Case case_, int state)
        {
            List<ASETools.CopyNumberVariation> tumorCopyNumber = null;
            if (case_.tumor_copy_number_file_id != "")
            {
                tumorCopyNumber = ASETools.CopyNumberVariation.ReadFile(case_.tumor_copy_number_filename);
                if (null == tumorCopyNumber)
                {
                    Console.WriteLine("unable to read copy number variation for case " + case_.case_id + " from file " + case_.tumor_copy_number_filename);
                    return;
                }
            }

            var asv = ASETools.AnnotatedVariant.readFile(case_.tentative_annotated_selected_variants_filename);
            if (null == asv)
            {
                Console.WriteLine("Unable to read tentative annotated selected variants for case " + case_.case_id + " from file " + case_.tentative_annotated_selected_variants_filename);
                return;
            }

            ASETools.AnnotatedVariant.writeFile(ASETools.GetDirectoryFromPathname(case_.tentative_annotated_selected_variants_filename) + @"\" + case_.case_id + ASETools.tentativeASVsWithoutCNVsExtension,
            asv.Where(_ => tumorCopyNumber.Where(cnv => cnv.OverlapsLocus(ASETools.chromosomeNameToNonChrForm(_.contig), _.locus, _.locus + 1)).Where(r => Math.Abs(r.Segment_Mean) > 1.0).ToList().Count() == 0).ToList());
        } // HandleOneCase
    }
}
