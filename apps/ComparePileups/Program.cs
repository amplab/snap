using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;

namespace ComparePileups
{
    class Program
    {
        class BAMAndSAMPair
        {
            public string bamFilename;
            public string nameSortedSamFilename;

            public FileStream nameSortedSamFile;

            public List<ASETools.SAMLine> readsMappingLocus;

            public ASETools.SAMLine getMatchingRead(ASETools.SAMLine query)
            {
                var matchingReads = readsMappingLocus.Where(_ => _.qname == query.qname && _.seqIfMappedForward == query.seqIfMappedForward).ToList();
                if (matchingReads.Count() == 0)
                {
                    return null;
                } else if (matchingReads.Count() == 1)
                {
                    return matchingReads[0];
                } else
                {
                    Console.WriteLine("!!!input file " + bamFilename + " has more than one read matching " + query.qname);
                    return matchingReads[0];
                }
            } // getMatchingRead
        } // BAMAndSAMPair
        static void Main(string[] args)
        {
            if (args.Count() < 3 || args.Count() > 5)
            {
                Console.WriteLine("usage: ComparePileups locus alignment0.bam alignment1.bam {alignment0.name-sorted.sam {alignment1.name-sorted.sam}}");
                Console.WriteLine("If you don't provide the name-sorted SAM files, you won't find out where reads that mapped elsewhere went.");
                return;
            } // If we don't have the right number of args

            var locus = args[0];
            var inputFiles = new BAMAndSAMPair[2];

            for (int i = 0; i < 2; i++)
            {
                inputFiles[i] = new BAMAndSAMPair();
                inputFiles[i].bamFilename = args[i + 1];
            }
  
            if (args.Count() > 3)
            {
                inputFiles[0].nameSortedSamFilename = args[3];
                if (args.Count() > 4)
                {
                    inputFiles[1].nameSortedSamFilename = args[4];
                }
            }

            var locusFields = locus.Split(':');
            if (locusFields.Count() != 2)
            {
                Console.WriteLine("Locus must be of the form contig:pos but is instead " + locus);
                return;
            }
            var contig = locusFields[0];
            var pos = Convert.ToInt32(locusFields[1].Replace(",",""));
            if (pos < 1)
            {
                Console.WriteLine("Locus pos must be > 0");
                return;
            }

            for (int i = 0; i < 2; i++)
            {
                if (inputFiles[i].nameSortedSamFilename != null)
                {
                    inputFiles[i].nameSortedSamFile = File.OpenRead(inputFiles[i].nameSortedSamFilename);
                    if (inputFiles[i].nameSortedSamFile == null)
                    {
                        Console.WriteLine("Unable to open " + inputFiles[i].nameSortedSamFilename);
                        return;
                    }
                }
            } // foreach input file


            string samtoolsLocus = contig + ":" + pos + "-" + (pos + 1);    // For whatever reason, if you just list the locus you get a bunch of nonsense in addition to the real data.

            for (int i = 0; i < 2; i++)
            {
                //
                // Run samtools to get the reads mapping locus; we need to filter them since samtools can be somewhat generous in what it reports.
                //
                inputFiles[i].readsMappingLocus = ASETools.RunProcessAndGetOutput(@"c:\bolosky\bin\samtools.exe", "view " + inputFiles[i].bamFilename + " " + samtoolsLocus).
                                                  Select(_ => new ASETools.SAMLine(_)).
                                                  Where(_ => !_.isUnmapped() && _.rname == contig && _.pos <= pos && _.pos + _.basesOfReferenceMapped > pos).ToList();
            }
            
            //
            // Find the reads mapped by both files.  This does *not* mean they're mapped to the same pos; just that they both cover the locus in question.
            //
            var readsMappedByBothFiles = inputFiles[0].readsMappingLocus.Where(read0 => inputFiles[1].readsMappingLocus.Any(read1 => read1.qname == read0.qname && read0.seqIfMappedForward == read1.seqIfMappedForward)).ToList();   // Check seqIfMappedForward to make sure it's not a mate
            var readsOnlyIn = new List<ASETools.SAMLine>[2];

            Console.WriteLine(readsMappedByBothFiles.Count() + " reads map the locus in both files.");

            //
            // Report reads for each file that aren't mapped in the other
            //
            for (int i = 0; i < 2; i++)
            {
                readsOnlyIn[i] = inputFiles[i].readsMappingLocus.Where(readI => readsMappedByBothFiles.All(readOther => readI.qname != readOther.qname || readI.seqIfMappedForward != readOther.seqIfMappedForward)).ToList();

                if (readsOnlyIn[i].Count() != 0)
                {
                    Console.WriteLine();
                    Console.WriteLine("Reads mapping " + locus + " only in " + inputFiles[i].bamFilename + ":");

                    foreach (var readOnlyInI in readsOnlyIn[i])
                    {
                        if (inputFiles[1 - i].nameSortedSamFile != null)
                        {
                            var readsInOtherFile = ASETools.FindSAMLinesInNameSortedSAMFile(inputFiles[1 - i].nameSortedSamFile, readOnlyInI.qname).Select(_ => new ASETools.SAMLine(_)).Where(_ => _.seqIfMappedForward == readOnlyInI.seqIfMappedForward).ToList();
                            if (readsInOtherFile.Count() == 0)
                            {
                                Console.WriteLine("Read " + readOnlyInI.qname + " is missing from " + inputFiles[1 - i].nameSortedSamFilename);
                                continue;
                            } 
                            else if (readsInOtherFile.Count() > 1)
                            {
                                Console.WriteLine("Read " + readOnlyInI.qname + " occurs more than once in " + inputFiles[1 - i].nameSortedSamFilename + ", choosing one.");
                            }

                            var readInOtherFile = readsInOtherFile[0];
 
                            if (readInOtherFile.isUnmapped())
                            {
                                Console.WriteLine(readOnlyInI.qname + " is unmapped in the other file");
                            } else 
                            {
                                Console.WriteLine(readOnlyInI.qname + " which is mapped at MAPQ " + readOnlyInI.mapq + ", NM " + readOnlyInI.NM() + 
                                    " in the other file is mapped to " + readInOtherFile.rname + ":" + readInOtherFile.pos + " at mapq " + readInOtherFile.mapq + " and NM " + readInOtherFile.NM());
                            }
                        } 
                        else
                        {
                            // no name sorted sam file for the other input, just report the read name
                            Console.WriteLine(readOnlyInI.qname);
                        }
                    }
                } else
                {
                    Console.WriteLine();
                    Console.WriteLine("File " + inputFiles[i].bamFilename + " has no reads mapping " + locus + " that aren't also mapping it in the other input file.");
                }
            } // each input file

            // Now look for differences.  Put each read into exactly one of these, taking condidions in order.
            var differInPOS = new List<ASETools.SAMLine>();
            var differInCIGAR = new List<ASETools.SAMLine>();
            var differInMAPQ = new List<ASETools.SAMLine>();

            var identicalReads = new List<ASETools.SAMLine>();

            foreach (var readIn0 in readsMappedByBothFiles)
            {
                var mateRead = inputFiles[1].readsMappingLocus.Where(_ => _.qname == readIn0.qname && _.seqIfMappedForward == readIn0.seqIfMappedForward).ToList()[0];

                var differentPOS = mateRead.pos != readIn0.pos;
                var differentCIGAR = mateRead.cigar != readIn0.cigar;
                var differentMAPQ = Math.Min(mateRead.mapq,60) != Math.Min(readIn0.mapq,60);    // SNAP uses MAPQ 70, others 60 as their max.  Smooth over that difference.

                if (differentPOS)
                {
                    differInPOS.Add(readIn0);
                }
                else if (differentCIGAR)
                {
                    differInCIGAR.Add(readIn0);
                }
                else if (differentMAPQ)
                {
                    differInMAPQ.Add(readIn0);
                }
                else
                {
                    identicalReads.Add(readIn0);
                }
            }

            Console.WriteLine();
            Console.WriteLine("Of the reads that map the location from both input files, " + differInPOS.Count() + " differ in mapped location, " + 
                differInCIGAR.Count() + " have different CIGAR strings, " + differInMAPQ.Count() + " differ in MAPQ, and " + identicalReads.Count() + " are identical");
            Console.WriteLine();

            if (differInPOS.Count() != 0)
            {
                Console.WriteLine();
                Console.WriteLine("Reads that differ in mapped location: ");
                foreach (var read0 in differInPOS)
                {
                    var read1 = inputFiles[1].getMatchingRead(read0);
                    Console.WriteLine(read0.qname + " file 0 maps it to " + read0.rname + ":" + read0.pos + " and file 1 maps it to " + read1.rname + ":" + read1.pos);
                }
                Console.WriteLine();
            }

            if (differInCIGAR.Count() != 0)
            {
                Console.WriteLine();
                Console.WriteLine("Reads that differ in CIGAR string:");

                foreach (var read0 in differInCIGAR)
                {
                    var read1 = inputFiles[1].getMatchingRead(read0);
                    Console.WriteLine(read0.qname + " mapped to " + read0.rname + ":" + read0.pos + " has CIGAR string " + read0.cigar + " in input 0 and " + read1.cigar + " in input 1.");
                }

                Console.WriteLine();
            }

            if (differInMAPQ.Count() != 0)
            {
                Console.WriteLine();
                Console.WriteLine("Reads that differ in MAPQ:");

                foreach (var read0 in differInMAPQ)
                {
                    var read1 = inputFiles[1].getMatchingRead(read0);
                    Console.WriteLine(read0.qname + " mapped to " + read0.rname + ":" + read0.pos + " has MAPQ " + read0.mapq + " in input 0 and " + read1.mapq + " in input 1.");
                }

                Console.WriteLine();
            }

        } // Main
    } // Program
} // Namespace
