using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.Diagnostics;
using System.Threading;
using System.IO;

//
// Take a set of long reads drawn from a sample and mapped to the reference and a set of short reads
// mapped to the long reads and remap the short reads to the original reference.
//

namespace RemapReadsToRef
{
    internal class Program
    {
        static int getRemappedPos(int pos, ASETools.SAMLine mappedLongRead)
        {
            if (pos == 0)
            {
                return 0;
            }

            int currentPos = mappedLongRead.pos;
            int remainingOffset = pos;

            //
            // Walk the cigar string of mappedLongRead and use that to adjust.
            //
            int offsetInCigar = 0;
            while (offsetInCigar < mappedLongRead.cigar.Length)
            {

            }

        } // getRemappedPos

        static string rewriteSAMLine(ASETools.SAMLine samLine, string readnamePrefix, ASETools.SAMLine mappedLongRead)
        {
            var fields = samLine.line.Split('\t');
            if (fields.Count() < 11)
            {
                throw new Exception("rewriteSAMLine: too few fields in SAM line: " + samLine);
            }

            // QNAME
            var outputLine = readnamePrefix + samLine.qname + "\t";

            if (mappedLongRead.isUnmapped())
            {
                //
                // The long read to which it mapped is itself unmapped.  Mark the read (and its mate) unmapped.
                //
                var newFlags = ((samLine.flag | ASETools.SAMLine.Unmapped) & ~ASETools.SAMLine.AllSegmentsProperlyAligned);
                if ((newFlags & ASETools.SAMLine.MultipleSegments) != 0)
                {
                    newFlags |= ASETools.SAMLine.NextUnmapped;
                }

                // FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN
                outputLine += newFlags + "\t*\t0\t0\t*\t*\t0\t0";// Do not add the final tab; it comes later
            } else
            {
                // FLAG
                outputLine += samLine.flag + "\t";

                // RNAME
                if (samLine.rname == "*")
                {
                    outputLine += "*\t";
                } else
                {
                    outputLine += mappedLongRead.rname + "\t";
                }

                // POS MAPQ CIGAR (in some cases CIGAR should be rewritten, but we'll just not bother)
                outputLine += getRemappedPos(samLine.pos, mappedLongRead) + "\t" + samLine.mapq + "\t" + samLine.cigar + "\t";

                // RNEXT
                if (samLine.rnext == "=" || samLine.rnext == "*")
                {
                    outputLine += samLine.rnext + "\t";
                }
                else
                {
                    outputLine += mappedLongRead.rname + "\t";
                }

                // PNEXT 
                outputLine += getRemappedPos(samLine.pnext, mappedLongRead) + "\t";

                // TLEN (this is a pain; for now, we'll just copy the tlen from the original read and hope the long read didn't map over any indels)
                outputLine += samLine.tlen; // Do not add the tab; it comes later
            }

            // SEQ, QUAL and the optional fields remain.  Copy them verbatim.

            for (int i = 9; i < fields.Length; i++)
            {
                outputLine += "\t" + fields[i];
            }

            return outputLine;            
        }

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            if (args.Length != 3)
            {
                Console.WriteLine("usage: RemapReadsToRef longReads.sam shortReadFilenameTemplate outputFile.sam");
                Console.WriteLine("shortReadFilenameTemplate will have _n.fastq added for n from 0 to whatever fails first.");
                return;
            }
            var commandLine = args[0] + " " + args[1] + " " + args[2];

            var longreadInputFilename = args[0];
            var shortReadFilenameTemplate = args[1];
            var outputFilename = args[2];

            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
            if (outputFile == null)
            {
                Console.WriteLine("Unable to open " + outputFilename + " for write.");
                return;
            }

            var longreadInputFile = ASETools.CreateStreamReaderWithRetryCompressedBasedOnFilename(longreadInputFilename);
            if (null == longreadInputFile)
            {
                Console.WriteLine("unable to open " + longreadInputFile + " for read");
                return;
            }

            var longReads = ASETools.SAMLine.ReadFromFile(longreadInputFile);
            var longReadsByName = new Dictionary<string, ASETools.SAMLine>();
            longReads.ForEach(_ => longReadsByName.Add(_.qname, _));    // These are unique because in addition to source contig and offset, they have a serial number appended to their name (i.e., _n#)

            //
            // Copy the SAM header from the long read file to the output, adding our @PG line
            //
            longreadInputFile.Close();
            longreadInputFile = ASETools.CreateStreamReaderWithRetryCompressedBasedOnFilename(longreadInputFilename);
            if (null == longreadInputFile)
            {
                Console.WriteLine("unable to open " + longreadInputFile + " for read (for second time)");
                return;
            }

            string line = "";
            while (null != (line = longreadInputFile.ReadLine()))
            {
                if (!line.StartsWith("@"))
                {
                    break;
                }

                outputFile.WriteLine(line);
            }
            longreadInputFile.Close();
            outputFile.WriteLine("@PG\tID:RemapReadsToRef\tPN:RemapReadsToRef\tCL:RemapReadsToRef.exe " + commandLine);

            int inputFileNumber = 0;

            foreach (var shortReadFilaname in Directory.EnumerateFiles(ASETools.GetDirectoryFromPathname(shortReadFilenameTemplate), ASETools.GetFileNameFromPathname(shortReadFilenameTemplate)))
            {
                var shortReadsFile = ASETools.CreateStreamReaderWithRetryCompressedBasedOnFilename(shortReadFilaname);
                var shortReads = ASETools.SAMLine.ReadFromFile(shortReadsFile);

                shortReadsFile.Close();

                foreach (var shortRead in shortReads)
                {
                    
                }
            }

        } // Main
    } // Program
} // namespace
