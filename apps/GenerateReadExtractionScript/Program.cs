using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;

namespace GenerateScriptFromVariants
{
    class Program
    {
        static void Main(string[] args)
        {
            var configuration = ASETools.ASEConfirguation.loadFromFile(args);

            if (configuration == null)
            {
                Console.WriteLine("Giving up because we were unable to load configuration");
                return;
            }

            if (configuration.commandLineArgs.Count() != 3 || configuration.commandLineArgs[1] != "-d" && configuration.commandLineArgs[1] != "-r")
            {
                Console.WriteLine("usage: GenerateReadExtractionScript case_id <-d|-r> outputScriptFilename");
                return;
            }

            var case_id = args[0];
            bool forDNA = configuration.commandLineArgs[1] == "-d";

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

            if (!cases.ContainsKey(configuration.commandLineArgs[0]))
            {
                Console.WriteLine("Unable to find case id " + case_id);
                return;
            }

            var case_ = cases[case_id];

            var selectedVariantsFilename = case_.selected_variants_filename;
            var extractedMAFLinesFilename = case_.extracted_maf_lines_filename;
            var inputBamFilename = forDNA ? case_.tumor_dna_filename : case_.tumor_rna_filename;
            var outputScriptFilename = configuration.commandLineArgs[2];

            if (selectedVariantsFilename == "" || extractedMAFLinesFilename == "" || inputBamFilename == "")
            {
                Console.WriteLine("GenerateScriptsFromVariants: at least one input is missing for case " + case_id);
                return;
            }

            int nVariants = 0;

            var outputScript = ASETools.CreateStreamWriterWithRetry(outputScriptFilename);

            var selectedVariantsInputFile = ASETools.CreateStreamReaderWithRetry(selectedVariantsFilename);
 
            string line = selectedVariantsInputFile.ReadLine();

            if (null == line)
            {
                Console.WriteLine("Empty selected variants file " + selectedVariantsFilename);
                return;
            }

            string headerLinePrefix = "SelectGermlineVariants v";
            if (line.Count() < headerLinePrefix.Count() + 2 || line.Substring(0,headerLinePrefix.Count()) != headerLinePrefix)
            {
                Console.WriteLine("Selected variants file " + selectedVariantsFilename + " contains invalid header line");
                return;
            }

            if (line.Substring(headerLinePrefix.Count(), 2) != "1.")
            {
                Console.WriteLine("Selected variants file " + selectedVariantsFilename + " is of the wrong version.");
                return;
            }


            bool seenDone = false;
            bool failed = false;

            while (null != (line = selectedVariantsInputFile.ReadLine()))
            {
                if (seenDone)
                {
                    Console.WriteLine("Selected variants file " + selectedVariantsFilename + " extends beyond **done**");
                    failed = true;
                    break;
                }

                if ("**done**" == line)
                {
                    seenDone = true;
                    continue;
                }

                var fields = line.Split('\t');
                if (fields.Count() < 3) {
                    Console.WriteLine("Bad line in selected variants file " + selectedVariantsFilename + ": " + line);
                    failed = true;
                    break;

                }

                string chromosomeName = fields[0];

                int snvPosition = 0;
                try {
                    snvPosition = Convert.ToInt32(fields[1]);
                } catch (FormatException) {
                    Console.WriteLine("Unparsable pos in selected variants file " + selectedVariantsFilename + ": " + line);
                    failed = true;
                    break;
                }

                outputScript.WriteLine("samtools view " + inputBamFilename + " " + chromosomeName + ":" + Math.Max(1, snvPosition - 200) + "-" + (snvPosition + 10) +
                    @" > " + (forDNA ? case_.tumor_dna_file_id : case_.tumor_rna_file_id) + "-" + chromosomeName + "-" + snvPosition);

                nVariants++;
            }

            if (!seenDone && !failed)
            {
                Console.WriteLine("Selected variants file " + selectedVariantsFilename + " is truncated.");
                failed = true;
            }

            if (failed) {
                //
                // Turn the output script into an empty file, which will in turn cause GenerateConsolodatedExtractedVariants not to produce a file, which will lead to us eventually noticing this error message.
                //

                outputScript.Close();
                outputScript = new StreamWriter(outputScriptFilename);
                return;
            }

            //
            // Now do reads around mutations in the MAF.
            //
            var mafLines = ASETools.MAFLine.ReadFile(case_.extracted_maf_lines_filename, case_.maf_file_id, false);

            foreach (var mafLine in mafLines)
            {
                outputScript.WriteLine("samtools view " + inputBamFilename + " " + mafLine.Chromosome + ":" + Math.Max(1, mafLine.Start_Position - 200) + "-" + (mafLine.End_Positon + 10) +
                    @" > " + (forDNA ? case_.tumor_dna_file_id : case_.tumor_rna_file_id) + "-" + mafLine.Chromosome + "-" + Math.Max(1, mafLine.Start_Position - 200) + "-" + (mafLine.End_Positon + 10));
            }

            if (0 == nVariants)
            {
                outputScript.WriteLine("This script intentionally left blank.");
            }

            outputScript.Close();
        }
    }
}
