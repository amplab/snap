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
            int nInputFiles = args.Count();
            var inputFileNames = args;
           
            if (nInputFiles == 0)
            {
                Console.WriteLine("AnalyzeConcordanceVCFs inputFile(s)");
                return;
            }

            var inputFiles = new StreamReader[nInputFiles];
            var mostRecentLines = new ASETools.ConcordanceVCFLine[nInputFiles];
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
                inputFiles[i] = ASETools.CreateStreamReaderWithRetryCompressedBasedOnFilename(inputFileNames[i]);
                if (inputFiles[i] == null)
                {
                    return;
                }

                mostRecentLines[i] = ASETools.ConcordanceVCFLine.GetNextVCFLine(inputFiles[i]);

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
                var lineIndicesAtRef = Enumerable.Range(0, nInputFiles - 1).Where(_ => mostRecentLines[_] != null && mostRecentLines[_] == refLine);

                if (mostRecentLines.Any(_ => _ == null || _ != refLine))
                {
                    Console.Write("Locus " + refLine.locus() + " is missing from these input files: ");
                    foreach (int missingInputFileIndex in Enumerable.Range(0,nInputFiles - 1).Where(_ => mostRecentLines[_] == null || mostRecentLines[_] != refLine))
                    {
                        Console.Write(inputFileNames[missingInputFileIndex] + ", ");
                    }
                    Console.WriteLine();
                }

                if (mostRecentLines.Where(_ => _ == refLine).Any(_ => _.callType != refLine.callType))
                {
                    Console.Write("Not all calls for " + refLine.locus() + " have the same call type (file, callType): ");
                    foreach (int lineIndex in lineIndicesAtRef)
                    {
                        Console.Write("(" + lineIndex + ":" + mostRecentLines[lineIndex].callType + ") ");
                    }
                    Console.WriteLine();
                }


                foreach(int lineIndex in lineIndicesAtRef) 
                {
                    do
                    {
                        mostRecentLines[lineIndex] = ASETools.ConcordanceVCFLine.GetNextVCFLine(inputFiles[lineIndex]);
                    } while (mostRecentLines[lineIndex] != null && !mostRecentLines[lineIndex].confidentRegion);

                    if (mostRecentLines[lineIndex] != null)
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
