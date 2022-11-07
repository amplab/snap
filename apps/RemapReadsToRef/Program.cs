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
        static int getRemappedPos(int posRelativeToLongRead, ASETools.SAMLine mappedLongRead)
        {
            if (posRelativeToLongRead == 0)
            {
                return 0;
            }

            var cigar = ASETools.ParseCIGARString(mappedLongRead.cigar);

            if (cigar.Count() == 0)
            {
                return 0;   // This is a "*" cigar string, which typically means unmapped.
            }
            
            int currentPos = 0;
            int remainingOffset = posRelativeToLongRead;

            for (int indexInCigar = 0; indexInCigar < cigar.Count(); indexInCigar++)
            {
                switch (cigar[indexInCigar].type)
                {
                    case ASETools.CIGARType.Match:
                    case ASETools.CIGARType.Unequal:
                    case ASETools.CIGARType.Equal:

                        var amount = Math.Min(cigar[indexInCigar].count, remainingOffset);
                        currentPos += amount;
                        remainingOffset -= amount;
                        break;

                }
            }
        } // getRemappedPos

        static string convertContigToRef(string mappedLongReadContig)
        {
            //
            // The long reads have names like contig_offset_n#.  So, just take up to the first _
            //
            if (mappedLongReadContig == "*")
            {
                return "*";
            }

            if (mappedLongReadContig == "=")
            {
                return "=";
            }

            return mappedLongReadContig.Substring(0, mappedLongReadContig.IndexOf("_"));
        }

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
                    outputLine += convertContigToRef(mappedLongRead.rname) + "\t";
                }

                // POS MAPQ CIGAR (in some cases CIGAR should be rewritten, but we'll just not bother)
                outputLine += getRemappedPos(samLine.pos, mappedLongRead) + "\t" + samLine.mapq + "\t" + samLine.cigar + "\t";

                // RNEXT
                outputLine += convertContigToRef(samLine.rnext);

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
