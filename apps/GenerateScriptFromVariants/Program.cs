using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ExpressionLib;
using System.IO;

namespace GenerateScriptFromVariants
{
    class Program
    {
        static void Main(string[] args)
        {
            if (args.Count() != 4)
            {
                Console.WriteLine("usage: GenerateScriptsFromVariants selectedVariantsFilename inputBamFilename outputScriptFilename chromosomePrefix");
                return;
            }

            var selectedVariantsFilename = args[0];
            var inputBamFilename = args[1];
            var outputScriptFilename = args[2];
            var chromosomePrefix = args[3];

            var outputScript = new StreamWriter(outputScriptFilename);

            StreamReader inputFile = null;
            try
            {
                inputFile = new StreamReader(selectedVariantsFilename);
            } catch (FileNotFoundException) {
                Console.WriteLine("Error opening input file " + selectedVariantsFilename + ", aborting.");
                return;
            }

            string line = inputFile.ReadLine();

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

            while (null != (line = inputFile.ReadLine()))
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

                string chromosomeName;
                if (fields[0].Count() > 2 && fields[0].Substring(0,2).ToLower() == "gl" || fields[0].Count() > 6 && fields[0].ToLower().Contains("gl")) {
                    chromosomeName = fields[0]; // So as not to generate "chrgl000..."
                } else {
                    if (fields[0].Count() > 3 && fields[0].Substring(0,3).ToLower() == "chr") {
                        chromosomeName = args[3] + fields[0].Substring(3);
                    } else {
                        chromosomeName = args[3] + fields[0];
                    }
                }

                int snvPosition = 0;
                try {
                    snvPosition = Convert.ToInt32(fields[1]);
                } catch (FormatException) {
                    Console.WriteLine("Unparsable pos in selected variants file " + selectedVariantsFilename + ": " + line);
                    failed = true;
                    break;
                }

                outputScript.WriteLine("samtools view " + args[1] + " " + chromosomeName + ":" + Math.Max(1, snvPosition - 200) + "-" + (snvPosition + 10) +
                    @" > " + ExpressionTools.GetAnalysisIdFromPathname(inputBamFilename) + "-" + chromosomeName + "-" + snvPosition);
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
            }

            outputScript.Close();
        }
    }
}
