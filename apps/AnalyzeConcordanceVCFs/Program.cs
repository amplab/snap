using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;

namespace AnalyzeConcordanceVCFs
{
    class Program
    {

        static void Main(string[] args)
        {
            var inputFileNames = new List<string>();
            var exclusionFileNames = new List<string>();

            var excludedLoci = new Dictionary<string, HashSet<int>>();  // chrom->hashset of excluded loci on that chromosome.

            int whichArg = 0;
            while (whichArg < args.Count())
            {
                if (args[whichArg] == "-x")
                {
                    if (whichArg + 1 >= args.Count())
                    {
                        Console.WriteLine("Exclusion file name must follow -x");
                        return;
                    }

                    exclusionFileNames.Add(args[whichArg + 1]);
                    whichArg += 2;
                } else
                {
                    inputFileNames.Add(args[whichArg]);
                    whichArg++;
                }
            } // for all the args

            var nInputFiles = inputFileNames.Count();
            if (nInputFiles == 0)
            {
                Console.WriteLine("AnalyzeConcordanceVCFs {-x exclusionList} inputFile(s)");
                return;
            }

            var inputFiles = inputFileNames.Select(_ => ASETools.CreateStreamReaderWithRetryCompressedBasedOnFilename(_)).ToList();
            if (inputFiles.Any(_ => _ == null))
            {
                return;
            }

            foreach (var exclusionFileName in exclusionFileNames)
            {
                var exclusionFile = ASETools.CreateStreamReaderWithRetry(exclusionFileName);
                if (null == exclusionFile)
                {
                    return;
                }

                string line;
                while (null != (line = exclusionFile.ReadLine()))
                {
                    var fields = line.Split(':');
                    if (fields.Count() != 2)
                    {
                        Console.WriteLine("Unparsable line in exclusion file " + exclusionFileName + ": " + line);
                        return;
                    }

                    if (!excludedLoci.ContainsKey(fields[0]))
                    {
                        excludedLoci.Add(fields[0], new HashSet<int>());
                    }

                    excludedLoci[fields[0]].Add(Convert.ToInt32(fields[1].Replace(",", "")));   // Replace is used to get rid of the commas in the human-readable coordinates
                } // while we have an input line

                exclusionFile.Close();
            } // foreach exclusion filename

            var mostRecentLines = inputFiles.Select(_ => ASETools.ConcordanceVCFLine.GetNextVCFLine(_)).ToList();
            var counts = new Dictionary<ASETools.CallType, int>[nInputFiles];

            for (int i = 0; i < nInputFiles; i++)
            {
                counts[i] = new Dictionary<ASETools.CallType, int>();
                foreach (var callType in ASETools.EnumUtil.GetValues<ASETools.CallType>())
                {
                    counts[i].Add(callType, 0);
                } // callType
            } // inputFile

            for (int i = 0; i < nInputFiles; i++)
            {
                while (mostRecentLines[i] != null && !mostRecentLines[i].confidentRegion)
                {
                    mostRecentLines[i] = ASETools.ConcordanceVCFLine.GetNextVCFLine(inputFiles[i]);
                } // while we have a non-confident line

                if (mostRecentLines[i] != null)
                {
                    counts[i][mostRecentLines[i].callType]++;
                }
            } // for each input file name

            while (mostRecentLines.Any(_ => _ != null))
            {
                var refLine = mostRecentLines.Where(_ => _ != null).Min();
                var lineIndicesAtRef = Enumerable.Range(0, nInputFiles).Where(_ => mostRecentLines[_] != null && mostRecentLines[_] == refLine);

                var excluded = excludedLoci.ContainsKey(refLine.chrom) && excludedLoci[refLine.chrom].Contains(refLine.pos);

                if (!excluded)
                {
 
                    if (mostRecentLines.Any(_ => _ == null || _ != refLine))
                    {
                        Console.Write("Locus " + refLine.locus() + " (" + refLine.callType + "," + refLine.insertLength() + ") is missing from these input files: ");
                        foreach (int missingInputFileIndex in Enumerable.Range(0, nInputFiles).Where(_ => mostRecentLines[_] == null || mostRecentLines[_] != refLine))
                        {
                            Console.Write(inputFileNames[missingInputFileIndex] + ", ");
                        }
                        Console.WriteLine();
                    }
                    else if (mostRecentLines.Where(_ => _ == refLine).Any(_ => _.callType != refLine.callType))
                    {
                        Console.Write("Not all calls for " + refLine.locus() + "(" + refLine.variantType + ") have the same call type (file, callType, insert len): ");
                        foreach (int lineIndex in lineIndicesAtRef)
                        {
                            Console.Write("(" + lineIndex + ":" + mostRecentLines[lineIndex].callType + "," + mostRecentLines[lineIndex].insertLength() + ") ");
                        }
                        Console.WriteLine();
                    }
                    else if (mostRecentLines.Where(_ => _ == refLine).Any(_ => _.callType != ASETools.CallType.TP))
                    {
                        Console.WriteLine("All calls for " + refLine.locus() + "(" + refLine.variantType + "," + refLine.insertLength() + ") are incorrect: " + refLine.callType);
                    }
                }


                foreach(int lineIndex in lineIndicesAtRef) 
                {
                    do
                    {
                        mostRecentLines[lineIndex] = ASETools.ConcordanceVCFLine.GetNextVCFLine(inputFiles[lineIndex]);
                    } while (mostRecentLines[lineIndex] != null && !mostRecentLines[lineIndex].confidentRegion);

                    if (mostRecentLines[lineIndex] != null && !excluded)
                    {
                        counts[lineIndex][mostRecentLines[lineIndex].callType]++;
                    }
                }
            } // while we have an input line

            Console.WriteLine();
            for (int i = 0; i < nInputFiles; i++)
            {
                Console.WriteLine(inputFileNames[i] + " " + counts[i][ASETools.CallType.TP] + " TP, " + counts[i][ASETools.CallType.FP] + " FP and, " + counts[i][ASETools.CallType.FN] + " FN");
            }
        } // Main
    } // Program
} // namespace
