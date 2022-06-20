using System;
using System.CodeDom;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Net;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using ASELib;

//
// Deal with namesorted, split SAM files that aren't exactly perfectly lined up.
// This can happen because of supplementary alignments.
//
// Find any reads that aren't in the whole set and write them to "extra" files.
//

namespace NamesortSplitFixup
{
    internal class Program
    {

        class ComparableSAMLine : ASETools.SAMLine, IComparable<ComparableSAMLine>
        {
            public ComparableSAMLine(string rawLine) : base(rawLine) { }

            public int CompareTo(ComparableSAMLine peer)
            {
                if (qname.CompareTo(peer.qname) == 0)
                {
                    if (seqIfMappedForward == peer.seqIfMappedForward)
                    {
                        return 0;
                    }

                    //
                    // Code for dealing with hard clipping of Ns by novoalign stolen from CompareAlignments.
                    //
                    var read0 = this;
                    var read1 = peer;
                    //
                    // Now check that the extra bases are all Ns.  Exactly one will be clipped either at the beginning or the end (but not both).
                    // We can't tell which end from seqIfMappedForward, since the read may be RC or not.  So just check both ways and if either works
                    // the reads match.
                    //
                    var minLength = Math.Min(read0.seqIfMappedForward.Length, read1.seqIfMappedForward.Length);

                    //
                    // Recall that Enumerable.All returns true when run on the empty set.
                    // clipped at the end
                    //
                    if (read0.seqIfMappedForward.Substring(0, minLength) == read1.seqIfMappedForward.Substring(0, minLength) &&
                            Enumerable.Range(minLength, Math.Min(0, read0.seqIfMappedForward.Length - minLength)).All(_ => read0.seqIfMappedForward[_] == 'N') &&
                            Enumerable.Range(minLength, Math.Min(0, read1.seqIfMappedForward.Length - minLength)).All(_ => read1.seqIfMappedForward[_] == 'N'))
                    {
                        return 0;
                    }

                    //
                    // clipped at the begining.
                    //
                    if (read0.seqIfMappedForward.Substring(read0.seqIfMappedForward.Length - minLength) != read1.seqIfMappedForward.Substring(read1.seqIfMappedForward.Length - minLength))
                    {
                        //
                        // The non-clipped parts don't match.
                        //
                        return read0.seqIfMappedForward.Substring(read0.seqIfMappedForward.Length - minLength).CompareTo(read1.seqIfMappedForward.Substring(read1.seqIfMappedForward.Length - minLength));
                    }

                    //
                    // Recall that Enumerable.All returns true when run on the empty set.
                    //
                    if (Enumerable.Range(0, read0.seqIfMappedForward.Length - minLength).All(_ => read0.seqIfMappedForward[_] == 'N') &&
                           Enumerable.Range(0, read1.seqIfMappedForward.Length - minLength).All(_ => read1.seqIfMappedForward[_] == 'N'))
                    {
                        return 0;
                    }

                    return read0.seqIfMappedForward.Substring(0, read0.seqIfMappedForward.Length - minLength).CompareTo(read1.seqIfMappedForward.Substring(0, read1.seqIfMappedForward.Length - minLength));
                } // if the qnames match

                return qname.CompareTo(peer.qname);
            } // CompareTo
        } // ComparableSAMLine


        static bool getNextLine(StreamReader input, out string line, out ComparableSAMLine samLine)
        {
            samLine = null;
            line = input.ReadLine();
            while (line != null && (line.StartsWith("\t") || (samLine = new ComparableSAMLine(line)).isSupplementaryAlignment() || samLine.isSecondaryAlignment()))
            {
                line = input.ReadLine();
            }

            if (line == null)
            {
                samLine = null;
                return false;
            }

            return true;
        }

        static string tailBinary = @"c:\bolosky\bin\tail.exe";

        //
        // A class for getting lines from the end of a text file indexed by lines from the end.
        //
        class EndLinesOfFile
        {
            public EndLinesOfFile(string inputFilename_)
            {
                inputFilename = inputFilename_;
            }

            public string getLine(int countFromEnd) // 0-based, so getLine(0) returns the last line of the file
            {
                if (countFromEnd < 0)
                {
                    throw new Exception("countFromEnd can't be negative");
                }


                if (countFromEnd >= linesRequested)
                {
                    //
                    // Need more lines.  Start with 100, then grow exponentially to 100K, then grow lineraly.
                    //

                    if (linesRequested < 100)
                    {
                        linesRequested = 100;
                    } else if (linesRequested > 100000)
                    {
                        linesRequested += 100000;
                    } else
                    {
                        linesRequested *= 2;
                    }

                    lines = ASETools.RunProcessAndGetOutput(tailBinary, "-" + linesRequested + " " + inputFilename);    // This may not get linesRequested if the file isn't that big.
                    if (lines.Count() > 0 && lines[0] == "")
                    {
                        lines.RemoveAt(0);  // Sometimes tail prints a blank line to start.
                    }
                }

                if (countFromEnd >= lines.Count())
                {
                    return null;    // BOF (beginning of file)
                }

                return lines[lines.Count() - 1 - countFromEnd]; // Recall that lines[0] is linesRequested + 1 before the EOF
            }

            public readonly string inputFilename;
            List<string> lines = null;
            int linesRequested = 0;
        }

        static void Main(string[] args)
        {
            if (args.Count() < 2)
            {
                Console.WriteLine("usage: FixSplitNamesorted input1 input2 ... inputN");
                Console.WriteLine(@"'input' is the base name of a namesorted, split SAM file.So, for instance, \\machine\d$\dir\hg001.snap.namesorted");
                Console.WriteLine("This program will then read the inputs with .nnnnnnnn appended to the names and write into outputs");
                Console.WriteLine("<file>.extra any reads that aren't in all of the inputs");
                return;
            }

            var nInputFiles = args.Length;
            var inputFilenameBases = args;

            var inputsDone = new bool[nInputFiles];
            var outputFiles = new StreamWriter[nInputFiles];

            for (int i = 0; i < nInputFiles; i++ )
            {
                inputsDone[i] = false;

                outputFiles[i] = ASETools.CreateStreamWriterWithRetry(inputFilenameBases[i] + ".extra");
                if (outputFiles[i] == null)
                {
                    Console.WriteLine("Failed to open output file " + inputFilenameBases[i] + ".extra");
                    return;
                }
            } // for each input filename

            long nExtraLines = 0;

            int splitSerialNumber = /*BJB0*/1;
            while (inputsDone.All(_ => !_))
            {
                //
                // First do the beginning of the files.
                //
                var inputs = new List<StreamReader>();
                var latestLines = new List<string>();
                var latestSAMLines = new List<ComparableSAMLine>();

                for (int i = 0; i < nInputFiles; i++)
                {

                    inputs.Add(ASETools.CreateStreamReaderWithRetry(inputFilenameBases[i] + "." + String.Format("{0:00000000}", splitSerialNumber)));

                    if (inputs[i] == null)
                    {
                        inputsDone[i] = true;
                        continue;
                    }

                    string line;
                    ComparableSAMLine samLine;

                    if (!getNextLine(inputs[i], out line, out samLine))
                    {
                        inputsDone[i] = true;
                        break;
                    }

                    latestLines.Add(line);
                    latestSAMLines.Add(samLine);
                }

                if (inputsDone.Any(_ => _))
                {
                    break;  // If we're missing even one line, then everything left is extra.
                }

                while (latestSAMLines.Any(_ => _ != null && _.CompareTo(latestSAMLines[0]) != 0))
                {
                    //
                    // They don't all match.  Emit the smallest one to the appropriate extra file and
                    // replace it with the next line.
                    //
                    int smallest = 0;
                    for (int i = 1; i < nInputFiles; i++)
                    {
                        if (latestSAMLines[i].CompareTo(latestSAMLines[smallest]) < 0)
                        {
                            smallest = i;
                        }
                    } // for each input file

                    outputFiles[smallest].WriteLine(latestLines[smallest]);
                    nExtraLines++;

                    string line;
                    ComparableSAMLine samLine;

                    if (!getNextLine(inputs[smallest], out line, out samLine))
                    {
                        throw new Exception("Write me");    // Hit EOF in front of file processing.
                    }


                    latestLines[smallest] = line;
                    latestSAMLines[smallest] = samLine;
                }

                if (latestSAMLines.Any(_ => _ == null))
                {
                    throw new Exception("Didn't find match in input file chunk " + splitSerialNumber);
                }

                //
                // Now do the end of the files.
                //
                inputs.ForEach(_ => _.Close());
                inputs = null;

                var ends = new List<EndLinesOfFile>();
                var nextLineToConsider = new List<int>();
                for (int i = 0; i < nInputFiles; i++)
                {
                    ends.Add(new EndLinesOfFile(inputFilenameBases[i] + "." + String.Format("{0:00000000}", splitSerialNumber)));
                    latestLines[i] = ends[i].getLine(0);
                    latestSAMLines[i] = new ComparableSAMLine(latestLines[i]);
                    nextLineToConsider.Add(0);
                }

                while (latestSAMLines.Any(_ => _ != null && _.CompareTo(latestSAMLines[0]) != 0))
                {
                    //
                    // They don't all match.  Emit the largest one to the appropriate extra file and
                    // replace it with the next line.
                    //
                    int largest = 0;
                    for (int i = 1; i < nInputFiles; i++)
                    {
                        if (latestSAMLines[i].CompareTo(latestSAMLines[largest]) > 0)
                        {
                            largest = i;
                        }
                    } // for each input file

                    outputFiles[largest].WriteLine(latestLines[largest]);
                    nExtraLines++;

                    nextLineToConsider[largest]++;

                    if (null == (latestLines[largest] = ends[largest].getLine(nextLineToConsider[largest])))
                    {
                        throw new Exception("Write me");    // Hit EOF in front of file processing.
                    }

                    latestSAMLines[largest] = new ComparableSAMLine(latestLines[largest]);
                }
            }

            //
            // Go along emitting extras until all lines match.
            //

        } // Main
    } // Program
} // namespace
