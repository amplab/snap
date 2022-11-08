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

            posRelativeToLongRead--; // Because it's 1-based.

            var cigar = ASETools.ParseCIGARString(mappedLongRead.cigar);

            if (cigar.Count() == 0)
            {
                return 0;   // This is a "*" cigar string, which typically means unmapped.
            }
            
            int currentPos = 0;
            int remainingOffset = posRelativeToLongRead;

            for (int indexInCigar = 0; indexInCigar < cigar.Count(); indexInCigar++)
            {
                if (remainingOffset < 0)
                {
                    throw new Exception("greRemappedPos: remaining offset < 0: " + posRelativeToLongRead + " " + mappedLongRead.cigar);
                }

                if (remainingOffset == 0)
                {
                    break;
                }

                var amount = Math.Min(cigar[indexInCigar].count, remainingOffset);

                switch (cigar[indexInCigar].type)
                {
                    case ASETools.CIGARType.Match:
                    case ASETools.CIGARType.Unequal:
                    case ASETools.CIGARType.Equal:

                        currentPos += amount;
                        remainingOffset -= amount;
                        break;

                    case ASELib.ASETools.CIGARType.Insertion:
                    case ASETools.CIGARType.SoftClip:
                        //
                        // Bases in an insertion or soft clip do not correspond to bases in the referece, so we eat the
                        // bases in the read without advancing in the reference.
                        //
                        remainingOffset -= amount;
                        break;

                    case ASETools.CIGARType.Deletion:
                        //
                        // A deletion is bases in the reference that don't correspond to references in the read, 
                        // so we move along in the reference without using any bases in the read.
                        //
                        currentPos += cigar[indexInCigar].count; // NOT amount: we might have a deletion bigger than the remaining bases in the read (i.e., 10M10D3M)
                        break;

                    case ASETools.CIGARType.HardClip:
                        //
                        // In a hard clip the bases are deleted from the read and don't correspond to the reference, so
                        // we just ignore the cigar entry entirely.
                        //
                        break;

                    default:
                        throw new Exception("Unimplemented CIGAR string type: " + mappedLongRead.cigar);

                } // switch (cigar[indexInCigar].type)
            } // for each cigar element

            return currentPos + mappedLongRead.pos;
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
                if (samLine.rnext == "*" || samLine.rnext == "=")
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

            var longReads = ASETools.SAMLine.ReadFromFile(longreadInputFile, true, true);
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

            var nShortReadSets = Directory.EnumerateFiles(ASETools.GetDirectoryFromPathname(shortReadFilenameTemplate), ASETools.GetFileNameFromPathname(shortReadFilenameTemplate)).Count();
            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "long reads", nShortReadSets, out nPerDot);

            int nProcessed = 0;
            long nShortReads = 0;

            bool rcNoticePrinted = false;

            foreach (var shortReadFilaname in Directory.EnumerateFiles(ASETools.GetDirectoryFromPathname(shortReadFilenameTemplate), ASETools.GetFileNameFromPathname(shortReadFilenameTemplate)))
            {
                //
                // This assumes that there will be a short 
                // Do the dot before actually processing it, so that we don't have to handle printing the dot when all the short reads are unmapped.
                //
                nProcessed++;
                if (nProcessed % nPerDot == 0)
                {
                    Console.Write(".");
                }

                var shortReadsFile = ASETools.CreateStreamReaderWithRetryCompressedBasedOnFilename(shortReadFilaname);
                var shortReads = ASETools.SAMLine.ReadFromFile(shortReadsFile);

                nShortReads += shortReads.Count();

                if (shortReads.All(_ => _.isUnmapped()))
                {
                    //
                    // Probably a long read of all Ns or something.  Don't bother writing them out.
                    //
                    continue;
                }

                //
                // The RNAME of any mapped read should be the same as the QNAME of the long read to which it's mapped.
                //
                var aMappedRead = shortReads.Where(_ => !_.isUnmapped()).First();

                var longRead = longReadsByName[aMappedRead.rname];

                if (!rcNoticePrinted && longRead.isRC())
                {
                    Console.WriteLine("Long read " + longRead.qname + " is mapped RC");
                    rcNoticePrinted = true;
                }

                shortReadsFile.Close();

                foreach (var shortRead in shortReads)
                {
                    outputFile.WriteLine(rewriteSAMLine(shortRead, longRead.qname + "_", longRead));
                }
            } // for each short read file

            outputFile.Close();

            Console.WriteLine();
            Console.WriteLine("Processed " + nShortReads + " short reads mapped to " + longReadsByName.Count() + " long reads in " + ASETools.ElapsedTime(timer));

        } // Main
    } // Program
} // namespace
