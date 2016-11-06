using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ExpressionLib;

namespace ApplyBonferroniCorrection
{
    class Program
    {
        const string baseDirectory = @"f:\temp\expression\";
        const string addedExtension = "_bonferroni.txt";

        static void ProcessFile(string filename)
        {
            if (filename.Count() < 5 || filename.Substring(filename.Count() - 4, 4).ToLower() != ".txt")
            {
                Console.WriteLine("ProcessFile: filename " + filename + " does not end in .txt");
                return;
            }

            string outputFilename = filename.Substring(0, filename.Count() - 4) + addedExtension;
            string[] inputLines;

            try
            {
                inputLines = File.ReadAllLines(filename);
            }
            catch
            {
                Console.WriteLine("Error opening " + filename + ", ignoring.");
                return;
            }

            if (inputLines.Count() < 2) {
                Console.WriteLine("Input file " + filename + " appears to be truncated.");
                return;
            }

            var fields = inputLines[1].Split('\t');

            int nFields = fields.Count();
            if (nFields < 6)
            {
                Console.WriteLine("Input file " + filename + " doesn't appear to have enough fields.");
            }

            double[] validCount = new double[nFields - 5];
            for (int i = 0; i < nFields - 5; i++)
            {
                validCount[i] = 0;
            }

            for (int whichLine = 1; whichLine < inputLines.Count(); whichLine++)
            {
                fields = inputLines[whichLine].Split('\t');
                if (fields.Count() != nFields)
                {
                    Console.WriteLine("Input file " + filename + " has inconsistent field count on line " + whichLine + ": " + inputLines[whichLine]);
                    return;
                }
                for (int whichField = 1; whichField < nFields - 4; whichField++)    // First field is gene name, last four are nTumors, nZero, nOne, nMore
                {
                    if (fields[whichField] != "*")
                    {
                        validCount[whichField - 1]++;
                    }
                }
            }

            var outputFile = new StreamWriter(outputFilename);
            outputFile.WriteLine(inputLines[0] + "\tMin p");    // The header

            for (int whichLine = 1; whichLine < inputLines.Count(); whichLine++)
            {
                fields = inputLines[whichLine].Split('\t');

                outputFile.Write(ExpressionTools.ConvertToExcelString(fields[0]) + "\t");    // Gene name, with = and quotes around it to avoid having Excel turn gene names into dates.

                bool foundAny = false;
                double minP = 10000000;
                for (int whichField = 1; whichField < nFields - 4; whichField++)    // First field is gene name, last four are nTumors, nZero, nOne, nMore
                {
                    if (fields[whichField] == "*")
                    {
                        outputFile.Write("*\t");
                    }
                    else
                    {
                        try
                        {
                            double value = Convert.ToDouble(fields[whichField]) * validCount[whichField - 1];
                            minP = Math.Min(minP, value);
                            foundAny = true;
                            outputFile.Write(value  + "\t");
                        } catch (FormatException) {
                            Console.WriteLine("Unparsable field " + fields[whichField] + " in line " + whichLine + " of file " + filename);
                            outputFile.Close();
                            File.Delete(outputFilename);
                            return;
                        }
                    }
                } // for each value field.

                outputFile.WriteLine(fields[nFields - 4] + "\t" + fields[nFields - 3] + "\t" + fields[nFields - 2] + "\t" + fields[nFields - 1] + "\t" + (foundAny ? Convert.ToString(minP) : "*"));   // The counts that are Bonferroni corrected
            }


            outputFile.Close();
        }

        static void Main(string[] args)
        {
            //
            // Get the list of input files.  Note that these lists will also include the output files, so we need to 
            // skip them.
            //
            var inputFiles = Directory.GetFiles(baseDirectory, "ExpressionDistributionByMutationCount*.txt").ToList();
            inputFiles.AddRange(Directory.GetFiles(baseDirectory, "AlleleSpecificExpressionDistributionByMutationCount*.txt"));

            foreach (var inputFile in inputFiles)
            {
                if (inputFile.ToLower().IndexOf(addedExtension) == -1)
                {
                    ProcessFile(inputFile);
                }
            }
        }
    }
}
