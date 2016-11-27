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
        static string baseDirectory = @"f:\temp\expression\";
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

            var headerFields = inputLines[0].Split('\t');

            int nFields = headerFields.Count();
            if (nFields < 6)
            {
                Console.WriteLine("Input file " + filename + " doesn't appear to have enough fields.");
            }

            double[] validCount = new double[nFields];
            bool[] fieldsToConvert = new bool[nFields];
            bool[] zeroValueFields = new bool[nFields];
            bool[] oneValueFields = new bool[nFields];
            int overallValidCount = 0;
            for (int i = 0; i < nFields; i++)
            {
                validCount[i] = 0;
                fieldsToConvert[i] = headerFields[i].IndexOf("vs.") != -1;

                //
                // Note the extra space in these header fields.  Sorry.
                //
                zeroValueFields[i] = headerFields[i].IndexOf("0 mutation  mean") != -1 || headerFields[i].IndexOf("0 mutation mean") != -1 ||
                    headerFields[i].IndexOf("0 mutation  exclusive mean") != -1 || headerFields[i].IndexOf("0 mutation exclusive mean") != -1 ||
                    headerFields[i].IndexOf("0 mutation  mu") != -1 || headerFields[i].IndexOf("0 mutation mu") != -1 ||
                    headerFields[i].IndexOf("0 mutation  exclusive mu") != -1 || headerFields[i].IndexOf("0 mutation exclusive mu") != -1;

                oneValueFields[i] = (headerFields[i].IndexOf("1 mutation  mean") != -1 || headerFields[i].IndexOf("1 mutation mean") != -1 ||
                    headerFields[i].IndexOf("1 mutation  exclusive mean") != -1 || headerFields[i].IndexOf("1 mutation exclusive mean") != -1 ||
                    headerFields[i].IndexOf("1 mutation  mu") != -1 || headerFields[i].IndexOf("1 mutation mu") != -1 ||
                    headerFields[i].IndexOf("1 mutation  exclusive mu") != -1 || headerFields[i].IndexOf("1 mutation exclusive mu") != -1) &&   // Exclude >1 cases, of which the 1 case is a substring
                    (headerFields[i].IndexOf(">1 mutation  mean") == -1 || headerFields[i].IndexOf(">1 mutation mean") == -1 ||
                    headerFields[i].IndexOf(">1 mutation  exclusive mean") == -1 || headerFields[i].IndexOf(">1 mutation exclusive mean") == -1 ||
                    headerFields[i].IndexOf(">1 mutation  mu") == -1 || headerFields[i].IndexOf(">1 mutation mu") == -1 ||
                    headerFields[i].IndexOf(">1 mutation  exclusive mu") == -1 || headerFields[i].IndexOf(">1 mutation exclusive mu") == -1);
            }

            string[] fields;
            for (int whichLine = 1; whichLine < inputLines.Count(); whichLine++)
            {
                fields = inputLines[whichLine].Split('\t');
                if (fields.Count() != nFields)
                {
                    Console.WriteLine("Input file " + filename + " has inconsistent field count on line " + whichLine + ": " + inputLines[whichLine]);
                    return;
                }
                for (int whichField = 0; whichField < nFields; whichField++)    // First field is gene name, last four are nTumors, nZero, nOne, nMore
                {
                    if (fieldsToConvert[whichField] && fields[whichField] != "*")
                    {
                        validCount[whichField - 1]++;
                        overallValidCount++;
                    }
                }
            }

            Console.WriteLine("File " + filename + " has overall valid count (i.e. Bonferroni correction factor) of " + overallValidCount);

            //
            // We rearrange the order of the columns somewhat.  The input is geneName followed by sets of columns with uncorrected p values, raw value means or standard 
            // deviations or counts, followed by four columns that have a count of tumors by mutation count for this gene.  We then add in min p and min p at columns.
            // The final output format is gene name, min p, min p at, the four tumor count columns from the end of the input and then the rest of the columns from
            // the input in order.
            //
            var outputFile = ExpressionTools.CreateStreamWriterWithRetry(outputFilename);
            outputFile.Write(headerFields[0] + "\tMin p\tMin p at\t0 vs. 1 ratio at MinP\tSignificant@.01\tBest 0 vs. 1 ratio for significant results\tBest 0 vs. 1 ratio at");
            for (int i = nFields - 4; i < nFields; i++)
            {
                outputFile.Write("\t" + headerFields[i]);
            }
            for (int i = 1; i < nFields - 4; i++)
            {
                outputFile.Write("\t" + headerFields[i]);
            }
            outputFile.WriteLine();

            for (int whichLine = 1; whichLine < inputLines.Count(); whichLine++)
            {
                fields = inputLines[whichLine].Split('\t');

                outputFile.Write(ExpressionTools.ConvertToExcelString(fields[0]));    // Gene name, with = and quotes around it to avoid having Excel turn gene names into dates.

                string restOfOutputLine = "";   // Since we want to put minP at the beginning but don't know its value until the end, build up the output in this string.

                for (int whichField = nFields - 4; whichField < nFields; whichField++)
                {
                    restOfOutputLine += "\t" + fields[whichField];
                }
         
                bool foundAny = false;
                double minP = 10000000;
                double bestZeroVsOne = 1;
                int bestZeroVsOneAt = -1;
                bool significant = false;
                double zeroVsOne = 0;
                double zeroValue = -1;
                double oneValue = -1;
                bool justSetMinP = false;
                int minPAt = -1;

                for (int whichField = 1; whichField < nFields - 4; whichField++)    // First field is gene name, last four are nTumors, nZero, nOne, nMore
                {
                    if (fields[whichField] == "*" || (!fieldsToConvert[whichField] && !zeroValueFields[whichField] && !oneValueFields[whichField]) )
                    {
                        restOfOutputLine += "\t" + fields[whichField];
                    }
                    else
                    {
                        try
                        {
                            double value = Convert.ToDouble(fields[whichField]);
                            if (fieldsToConvert[whichField])
                            {
                                value *= overallValidCount;
                            }

                            if (fieldsToConvert[whichField] && value < minP)
                            {
                                minP = value;
                                minPAt = whichField;
                                justSetMinP = true;
                                zeroValue = -1;
                                oneValue = -1;
                                zeroVsOne = -1;
                                foundAny = true;
                            }

                            if (fieldsToConvert[whichField] && value <= .01)
                            {
                                significant = true;
                                zeroValue = -1;
                                oneValue = -1;
                            }

                            if (zeroValueFields[whichField])
                            {
                                zeroValue = value;
                            }

                            if (oneValueFields[whichField])
                            {
                                oneValue = value;
                            }

                            if (-1 != zeroValue && -1 != oneValue)
                            {
                                if (justSetMinP)
                                {
                                    if (0 == oneValue)
                                    {
                                        zeroVsOne = 100000; // Something big.
                                    }
                                    else
                                    {
                                        zeroVsOne = zeroValue / oneValue;
                                    }
                                }

                                if (significant && 0 != oneValue)
                                {
                                    double candidate = zeroValue / oneValue;
                                    if (candidate > 1)
                                    {
                                        candidate = 1 / candidate;
                                    }

                                    if (candidate < bestZeroVsOne)
                                    {
                                        bestZeroVsOne = candidate;
                                        bestZeroVsOneAt = whichField;
                                    }
                                }

                                justSetMinP = false;
                                zeroValue = -1;
                                oneValue = -1;
                                significant = false;
                            }

                            restOfOutputLine += "\t" + Convert.ToString(value);
                        } catch (FormatException) {
                            Console.WriteLine("Unparsable field " + fields[whichField] + " in line " + whichLine + " of file " + filename);
                            outputFile.Close();
                            File.Delete(outputFilename);
                            return;
                        }
                    }
                } // for each value field.

                if (foundAny)
                {
                    outputFile.Write("\t" + minP + "\t" + headerFields[minPAt] + "\t" + zeroVsOne + "\t" + ((minP < .01) ? "yes" : "no"));
                }
                else
                {
                    outputFile.Write("\t*\t*\t*\t*");
                }

                if (bestZeroVsOneAt != -1)
                {
                    outputFile.Write("\t" + bestZeroVsOne + "\t" + headerFields[bestZeroVsOneAt]);
                }
                else
                {
                    outputFile.Write("\t*\t*");
                }

                outputFile.WriteLine(restOfOutputLine);
            }


            outputFile.Close();
        }

        static void Main(string[] args)
        {
            if (args.Count() == 1)
            {
                baseDirectory = args[0];
            }
            else if (args.Count() != 0)
            {
                Console.WriteLine("usage: ApplyBonferroniCorrection {base directory}");
                return;
            }
            //
            // Get the list of input files.  Note that these lists will also include the output files, so we need to 
            // skip them.
            //
            var inputFiles = Directory.GetFiles(baseDirectory, "ExpressionDistributionByMutationCount*.txt").ToList();
            inputFiles.AddRange(Directory.GetFiles(baseDirectory, "AlleleSpecificExpressionDistributionByMutationCount*.txt"));

            foreach (var inputFile in inputFiles)
            {
                if (inputFile.ToLower().IndexOf(addedExtension) == -1 && inputFile.ToLower().IndexOf("_by_range_graphs") == -1)
                {
                    ProcessFile(inputFile);
                }
            }
        }
    }
}
