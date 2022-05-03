using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using ASELib;
using System.Net;
using System.Xml.Linq;
using System.IO;
using System.Runtime.Remoting;
using System.Threading;
using System.Net.Sockets;
using System.Data;
using System.Diagnostics.Eventing.Reader;

namespace CompareAlignments
{
    internal class Program
    {
        const int fuzziness = 50;   // How close is close enough (except that to make it equivalence classes we just divide by fuzziness, so fuzziness - 2 and fuzziness + 2 don't match)

        //
        // This class is used when we're doing wgsim to indicate whether the line is correctly aligned.
        // If not wgsim then its "correctlyAligned" bools are always just false.
        //
        class SAMLineWithCorrectAlignmentIndicator : ASETools.SAMLine
        {
            public SAMLineWithCorrectAlignmentIndicator(string rawLine_, bool wgsim): base(rawLine_)
            {
                if (!wgsim)
                {
                    return;
                }

                //
                // The "correct" alignment is determined from the read name.
                // The read name format is chr_forwardPos_reversePos.  The actual alignment you'll see for reversePos is reversePos - readLen + 1.
                // 
                // Alas, contigs contain "_" as their name, so we have to carefully parse here.
                //

                if (isUnmapped())
                {
                    return;
                }

                if (!qname.StartsWith(rname + "_"))
                {
                    //
                    // Different contigs.
                    //
                    return;
                }

                var fields = qname.Substring(rname.Length).Split('_');    // fields[0] will be "" since we know the first char of the substring is "_"

                if (fields[1][0] < '0' || fields[1][0] > '9' || fields[2][0] < '0' || fields[2][0] > '9') // Shortcutting the (very slow) exception path below
                {
                    return;
                }

                int pos0, pos1;

                try
                {
                    pos0 = Convert.ToInt32(fields[1]);
                    pos1 = Convert.ToInt32(fields[2]) - seq.Length + 1;
                } 
                catch
                {
                    //
                    // This is the case where the rname is a substring of the real contig, which means that fields[1] isn't really the pos
                    return;
                }

                correctlyAligned = pos0 / fuzziness == pos / fuzziness || pos1 / fuzziness == pos / fuzziness;
                corretlyAlignedExactly = pos0 == pos || pos1 == pos;
            }

            public bool correctlyAligned = false;       // Within fuzziness
            public bool corretlyAlignedExactly = false; // Exact
        }


        interface ReadGetter
        {
            SAMLineWithCorrectAlignmentIndicator GetRead();
        }

        class StreamReaderReadGetter : ReadGetter
        {
            public StreamReaderReadGetter(StreamReader reader_, bool wgsim_)
            {
                reader = reader_;
                wgsim = wgsim_;

                var thread = new Thread(() => threadWorker(this));
                thread.Start();
            }

            static void threadWorker(StreamReaderReadGetter getter)
            {
                getter.threadWorker();
            }

            void threadWorker()
            {
                string line;
                while ((line = reader.ReadLine()) != null)
                {
                    if (line.StartsWith("@"))
                    {
                        continue;
                    }

                    var samLine = new SAMLineWithCorrectAlignmentIndicator(line, wgsim);
                    if (samLine.isSupplementaryAlignment() || samLine.isSecondaryAlignment())
                    {
                        continue;   // ignore these
                    }
                    queue.Enqueue(samLine);
                }

                queue.TerminateWriter();
            }

            public SAMLineWithCorrectAlignmentIndicator GetRead()
            {
                SAMLineWithCorrectAlignmentIndicator line;
                if (!queue.Dequeue(out line))
                {
                    return null;
                }

                return line;
            } // GetRead

            ASETools.ThrottledParallelQueue<SAMLineWithCorrectAlignmentIndicator> queue = new ASETools.ThrottledParallelQueue<SAMLineWithCorrectAlignmentIndicator>(20000, 1);

            StreamReader reader;
            bool wgsim;
        } // StreamReaderReadGetter

        class BAMReadGetter : ReadGetter
        {
            public BAMReadGetter(string inputFilename_, bool wgsim_)
            {
                inputFilename = inputFilename_;
                wgsim = wgsim_;

                var thread = new Thread(() => threadWorker(this));
                thread.Start();
            }

            static void threadWorker(BAMReadGetter getter)
            {
                getter.threadWorker();
            }

            void threadWorker()
            {
                string samtoolsBinary = @"c:\bolosky\bin\samtools.exe";
                if (!File.Exists(samtoolsBinary))
                {
                    Console.WriteLine("Can't find samtools binary at " + samtoolsBinary);
                    throw new Exception("Can't find samtools binary at " + samtoolsBinary);
                }

                ASETools.RunProcess(samtoolsBinary, "view " + inputFilename, _ => GotReadLine(_));
                queue.TerminateWriter();
            }

            void GotReadLine(string rawLine)
            {
                //
                // Don't need to worry about header lines, since samtools doesn't print them.
                //
                var line = new SAMLineWithCorrectAlignmentIndicator(rawLine, wgsim);

                if (line.isSecondaryAlignment() || line.isSupplementaryAlignment()) // Drop secondary and supplementary alignments.
                {
                    return;
                }
                queue.Enqueue(line);
            }

            public SAMLineWithCorrectAlignmentIndicator GetRead()
            {
                SAMLineWithCorrectAlignmentIndicator line;
                if (!queue.Dequeue(out line))
                {
                    return null;
                }

                return line;
            } // GetRead


            ASETools.ThrottledParallelQueue<SAMLineWithCorrectAlignmentIndicator> queue = new ASETools.ThrottledParallelQueue<SAMLineWithCorrectAlignmentIndicator>(20000, 1);
            string inputFilename;
            bool wgsim;
        }


        delegate bool MatchFunction(int a, int b, bool wgsim);  // This must be equivalence classes.  i.e., matchFunction(a, b) && matchFunction(b, c) => matchFunction(a, c)
        static string partitionString(int nClasses, MatchFunction matchFunction, bool wgsim)
        {
            var result = "";

            var used = new bool[nClasses];
            for (int i = 0; i < nClasses; i++)
            {
                used[i] = false;
            }

            bool needBar = false;
            for (int i = 0; i < nClasses; i++)
            {
                if (used[i])
                {
                    continue;
                }

                used[i] = true;

                if (needBar)
                {
                    result += "|";
                } else
                {
                    needBar = true;
                }

                result += (char)('a' + i);

                for (int j = i+1; j < nClasses; j++)
                {
                    if (!used[j] && matchFunction(j, i, wgsim)) // !used is unnecessary and here just to avoid calling matchFunction when we don't have to
                    {
                        result += (char)('a' + j);
                        used[j] = true;
                    }
                } // for each potential partner
            } // for each input file

            return result;
        } // paritionString

        delegate bool SAMComparitor(SAMLineWithCorrectAlignmentIndicator read0, SAMLineWithCorrectAlignmentIndicator read1);    // a null read means use the "correct" alignment from wgsim

        static bool doesAlignmentMatchExactly(SAMLineWithCorrectAlignmentIndicator read0, SAMLineWithCorrectAlignmentIndicator read1)
        {
            if (read0 == null && read1 == null)
            {
                return true;
            }

            if (read0 == null || read1 == null)
            {
                var read = (read0 == null) ? read1 : read0;

                //
                // The null read is the "correct" alignment, which is determined from the read name in the other read and 
                // was attached to the read when it was ingested.
                //

                return read.corretlyAlignedExactly;
            }

            if (read0.isUnmapped() != read1.isUnmapped())
            {
                return false;
            }

            if (read0.isUnmapped())
            {
                return true;
            }

            return read0.pos == read1.pos && read0.rname == read1.rname;
        }
        static bool doesAlignmentMatch(SAMLineWithCorrectAlignmentIndicator read0, SAMLineWithCorrectAlignmentIndicator read1)
        {

            if (read0 == null && read1 == null)
            {
                return true;
            }

            if (read0 == null || read1 == null)
            {
                var read = (read0 == null) ? read1 : read0;

                //
                // The null read is the "correct" alignment, which is determined from the read name in the other read and 
                // was attached to the read when it was ingested.
                //

                return read.correctlyAligned;
            }

            if (read0.isUnmapped() != read1.isUnmapped())
            {
                return false;
            }

            if (read0.isUnmapped())
            {
                return true;
            }

            return read0.pos / fuzziness == read1.pos / fuzziness && read0.rname == read1.rname;
        }

        static bool doesAlignmentWithMAPQMatch(SAMLineWithCorrectAlignmentIndicator read0, SAMLineWithCorrectAlignmentIndicator read1)
        {
            return read0.isUnmappedOrMapq0() && read1.isUnmappedOrMapq0() || !read0.isUnmappedOrMapq0() && !read1.isUnmappedOrMapq0() && doesAlignmentMatch(read0, read1);
        }

        static bool doesAlignmentMAPQAndCigarMatch(SAMLineWithCorrectAlignmentIndicator read0, SAMLineWithCorrectAlignmentIndicator read1)
        {
            return doesAlignmentWithMAPQMatch(read0, read1) && (read0.mapq == read1.mapq || read0.cigar == read1.cigar);
        }

        static bool compareReads(Dictionary<string, List<SAMLineWithCorrectAlignmentIndicator>> [] unmatchedReads, SAMLineWithCorrectAlignmentIndicator example, int exampleFrom, int a, int b, SAMComparitor comparitor)
        {
            SAMLineWithCorrectAlignmentIndicator read0, read1;
            if (a == exampleFrom)
            {
                read0 = example;
            } else
            {
                read0 = unmatchedReads[a][example.qname].Where(_ => _.seqIfMappedForward == example.seqIfMappedForward).FirstOrDefault();
            }

            if (b == exampleFrom)
            {
                read1 = example;
            } else
            {
                read1 = unmatchedReads[b][example.qname].Where(_ => _.seqIfMappedForward == example.seqIfMappedForward).FirstOrDefault();
            }

            return comparitor(read0, read1);
        }

        static void WriteMatchClasses(List<StreamWriter> outputStreams, string header, Dictionary<string, long> matches, int nClasses)
        {
            ASETools.WriteLineToMultipleStreams(outputStreams, header);

            var total = matches.Select(_ => _.Value).Sum();

            var list = matches.ToList();
            list.Sort((a, b) => b.Value.CompareTo(a.Value));    // Backward comparison to pu the big ones on top

            var concordance = new long[nClasses, nClasses];

            foreach (var match in list)
            {
                ASETools.WriteLineToMultipleStreams(outputStreams, match.Key + "\t" + ASETools.NumberWithCommas(match.Value) + "\t" + ASETools.Percentage(match.Value, total, 2));

                foreach (var group in match.Key.Split('|')) 
                {
                    for (int i = 0; i < group.Count(); i++)
                    {
                        for (int j = 0; j < group.Count(); j++)
                        {
                            concordance[(int)(group[i] - 'a'), (int)(group[j] - 'a')] += match.Value;
                        } // j 
                    } // i
                } // group
            } // match

            ASETools.WriteLineToMultipleStreams(outputStreams);
            ASETools.WriteLineToMultipleStreams(outputStreams, "Concordance (raw):");
            for (int i = 0; i < nClasses; i++)
            {
                ASETools.WriteToMultipleStreams(outputStreams, "\t" + (char)('a' + i));
            }
            ASETools.WriteLineToMultipleStreams(outputStreams);

            for (int i = 0; i < nClasses; i++)
            {
                ASETools.WriteToMultipleStreams(outputStreams, (char)('a' + i) + "");
                for (int j = 0; j < nClasses; j++)
                {
                    ASETools.WriteToMultipleStreams(outputStreams, "\t" + concordance[i, j]);
                } // j
                ASETools.WriteLineToMultipleStreams(outputStreams);
            } // i

            ASETools.WriteLineToMultipleStreams(outputStreams);
            ASETools.WriteLineToMultipleStreams(outputStreams, "Concordance (percentage):");
            for (int i = 0; i < nClasses; i++)
            {
                ASETools.WriteToMultipleStreams(outputStreams, "\t" + (char)('a' + i));
            }
            ASETools.WriteLineToMultipleStreams(outputStreams);

            for (int i = 0; i < nClasses; i++)
            {
                ASETools.WriteToMultipleStreams(outputStreams, (char)('a' + i) + "");
                for (int j = 0; j < nClasses; j++)
                {
                    ASETools.WriteToMultipleStreams(outputStreams, "\t" + ASETools.Percentage(concordance[i, j], total, 2));
                } // j
                ASETools.WriteLineToMultipleStreams(outputStreams);
            } // i


        } // WriteMatchClasses

        class MatchGroup
        {
            public MatchGroup(SAMComparitor comparitor_, Dictionary<string,long> matched_, ASETools.ThrottledParallelQueue<List<SAMLineWithCorrectAlignmentIndicator>> inputQueue_, bool wgsim_, int nClasses_)
            {
                comparitor = comparitor_;
                matched = matched_;
                inputQueue = inputQueue_;
                wgsim = wgsim_;
                nClasses = nClasses_;

                thread = new Thread(() => Worker(this));
                thread.Start();
            }

            public void waitForDone()
            {
                thread.Join();
            }

            static void Worker(MatchGroup group)
            {
                group.Worker();
            }

            static bool compareReads(List<SAMLineWithCorrectAlignmentIndicator> reads, int a, int b, SAMComparitor comparitor)
            {
                //
                // If we have wgsim, then a and/or be can be reads.Count() (i.e., one too big).  That indicates
                // the "correct" alignment.  In that case, we pass in null for the read to comparitor
                //

                return comparitor((a == reads.Count()) ? null : reads[a], (b == reads.Count()) ? null : reads[b]);
            }

            void Worker()
            {
                List<SAMLineWithCorrectAlignmentIndicator> reads;
                while (inputQueue.Dequeue(out reads))
                {
                    var partition = partitionString(nClasses, (a, b, wgsim) => compareReads(reads, a, b, comparitor), wgsim);
                    if (!matched.ContainsKey(partition))
                    {
                        matched.Add(partition, 1);
                    }
                    else
                    {
                        matched[partition]++;
                    }
                } // while we have input read sets
            } // Worker

            Thread thread;
            SAMComparitor comparitor;
            ASETools.ThrottledParallelQueue<List<SAMLineWithCorrectAlignmentIndicator>> inputQueue;
            Dictionary<string, long> matched;
            bool wgsim;
            int nClasses;
        } // MatchGroup
        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            if (args.Length < 3 || args[0].ToLower().EndsWith(".sam") || args[0].ToLower().EndsWith(".bam"))
            {
                Console.WriteLine("usage: CompareAligments {-wgsim} outputFile inputFile1 inputFile2 <...inputFileN>");
                Console.WriteLine("output file cannot end with .sam or .bam (to keep you from accidentally ovewriting an input)");
                Console.WriteLine("-wgsim means that the reads came from wgsim and have the correct alignment as their name");
                return;
            }

            int nextArg = 0;
            bool wgsim = args[nextArg] == "-wgsim";
            if (wgsim)
            {
                nextArg++;
            }

            var outputFile = ASETools.CreateStreamWriterWithRetry(args[nextArg]);
            if (outputFile == null)
            {
                Console.WriteLine("Unable to open output file " + args[nextArg]);
                return;
            }
            nextArg++;

            var outputStreams = new List<StreamWriter>();
            outputStreams.Add(new StreamWriter(Console.OpenStandardOutput()));
            outputStreams.Add(outputFile);

            var nInputFiles = args.Length - nextArg;
            var unmatchedReads = new Dictionary<string, List<SAMLineWithCorrectAlignmentIndicator>>[nInputFiles];
            var queueDone = new bool[nInputFiles];

            var readerQueues = new ReadGetter[nInputFiles];

            for (int i = 0; i < nInputFiles; i++)
            {
                var inputFilename = args[i+nextArg]; 
                if (inputFilename.ToLower().EndsWith(".bam"))
                {
                    readerQueues[i] = new BAMReadGetter(inputFilename, wgsim);
                }
                else
                {
                    var inputFile = ASETools.CreateStreamReaderWithRetry(inputFilename);
                    if (null == inputFile)
                    {
                        Console.WriteLine("Unable to open " + inputFilename);
                        return;
                    }

                    readerQueues[i] = new StreamReaderReadGetter(inputFile, wgsim);
                }

                unmatchedReads[i] = new Dictionary<string, List<SAMLineWithCorrectAlignmentIndicator>>();
                queueDone[i] = false;
            } // start reader threads

            var nClasses = nInputFiles + (wgsim ? 1 : 0);    // if wgsim, then the correct alignment is an extra class

            //
            // These map partition strings to counts.  They're sparse because the space of possible partitons grows hyperexponentially in nInputFiles, so most won't be used
            // if we have many aligners.
            //
            var matchedAlignmentExactly = new Dictionary<string, long>();
            var matchedAlignment = new Dictionary<string, long>();
            var matchedAlignmentWithMAPQ = new Dictionary<string, long>();
            var matchedAlignmentMAPQAndCigar = new Dictionary<string, long>();

            var matchGroups = new List<MatchGroup>();
            var queues = new List<ASETools.ThrottledParallelQueue<List<SAMLineWithCorrectAlignmentIndicator>>>();
 
            { // scope to make the definition of queue be local
                var queue = new ASETools.ThrottledParallelQueue<List<SAMLineWithCorrectAlignmentIndicator>>(5000, 1);
                queues.Add(queue);
                matchGroups.Add(new MatchGroup(doesAlignmentMatchExactly, matchedAlignmentExactly, queue, wgsim, nClasses));

                queue = new ASETools.ThrottledParallelQueue<List<SAMLineWithCorrectAlignmentIndicator>>(5000, 1);
                queues.Add(queue);
                matchGroups.Add(new MatchGroup(doesAlignmentMatch, matchedAlignment, queue, wgsim, nClasses));

                queue = new ASETools.ThrottledParallelQueue<List<SAMLineWithCorrectAlignmentIndicator>>(5000, 1);
                queues.Add(queue);
                matchGroups.Add(new MatchGroup(doesAlignmentWithMAPQMatch, matchedAlignmentWithMAPQ, queue, false, nInputFiles));    // wgsim doesn't make sense for MAPQ, since the "correct" answer doesn't come with MAPQ

                queue = new ASETools.ThrottledParallelQueue<List<SAMLineWithCorrectAlignmentIndicator>>(5000, 1);
                queues.Add(queue);
                matchGroups.Add(new MatchGroup(doesAlignmentMAPQAndCigarMatch, matchedAlignmentMAPQAndCigar, queue, false, nInputFiles)); // wgsim doesn't make sense for CIGAR, since the "correct" answer doesn't have one (or MAPQ)
            }


            int nextQueue = 0;
            long nReads = 0;
            long nMatches = 0;

            while (!queueDone.All(_ => _))
            {
                if (!queueDone[nextQueue])
                {
                    SAMLineWithCorrectAlignmentIndicator samLine;
                    if (null == (samLine = readerQueues[nextQueue].GetRead()))
                    {
                        queueDone[nextQueue] = true;
                        nextQueue = (nextQueue + 1) % nInputFiles;
                        continue;
                    } // if the queue was newly done

                    nReads++;
                    ///*BJB*/ if (nReads >= 5000000) break;
                    if (nReads % 1000000 == 0)
                    {
                        Console.WriteLine(ASETools.ElapsedTime(timer) + ": " + ASETools.NumberWithCommas(nReads)
                                          + " (" + ASETools.Percentage(nMatches * nInputFiles, nReads, 2) + " matched), " +
                                          ASETools.NumberWithCommas(nReads - nMatches * nInputFiles) + " unmatched");
                    }

                    //
                    // See if this line is now read in all input files
                    //
                    if (Enumerable.Range(0,nInputFiles).All(_ => (_ == nextQueue) || unmatchedReads[_].ContainsKey(samLine.qname) && unmatchedReads[_][samLine.qname].Any(x => x.seqIfMappedForward == samLine.seqIfMappedForward))) {
                        nMatches++;

                        var reads = new List<SAMLineWithCorrectAlignmentIndicator>();
                        for (int i = 0; i < nInputFiles; i++)
                        {
                            if (i == nextQueue)
                            {
                                reads.Add(samLine);
                            } else
                            {
                                reads.Add(unmatchedReads[i][samLine.qname].Where(_ => _.seqIfMappedForward == samLine.seqIfMappedForward).ToList()[0]);
                            }
                        }

                        foreach(var queue in queues)
                        {
                            queue.Enqueue(reads);
                        }

                        for (int inputToRemove = 0; inputToRemove < nInputFiles; inputToRemove++)
                        {
                            if (inputToRemove != nextQueue)
                            {
                                for (int i = 0; i < unmatchedReads[inputToRemove][samLine.qname].Count(); i++)
                                {
                                    if (unmatchedReads[inputToRemove][samLine.qname][i].seqIfMappedForward == samLine.seqIfMappedForward)
                                    {
                                        unmatchedReads[inputToRemove][samLine.qname].RemoveRange(i, 1);
                                        break;  // This is intentional; we could have two identical reads in which case we only want to remove one
                                    }
                                }
                            }

                            if (unmatchedReads[inputToRemove].ContainsKey(samLine.qname) && unmatchedReads[inputToRemove][samLine.qname].Count() == 0)
                            {
                                unmatchedReads[inputToRemove].Remove(samLine.qname);
                            }
                        }
                    } // matched everywhere
                    else
                    {
                        if (!unmatchedReads[nextQueue].ContainsKey(samLine.qname))
                        {
                            unmatchedReads[nextQueue].Add(samLine.qname, new List<SAMLineWithCorrectAlignmentIndicator>());
                        }

                        unmatchedReads[nextQueue][samLine.qname].Add(samLine);
                    } // Not yet matched everywhere
                }

                nextQueue = (nextQueue + 1) % nInputFiles;

            } // while we have SAM lines to process

            foreach (var queue in queues)
            {
                queue.TerminateWriter();
            }

            foreach (var matchGroup in matchGroups)
            {
                matchGroup.waitForDone();
            }

            //
            // Now write the output
            //

            for (nextQueue = 0; nextQueue < nInputFiles; nextQueue++)
            {
                outputFile.WriteLine("There are " + unmatchedReads[nextQueue].Count() + " unmatched reads from " + args[nextQueue + 1] + ".  Here are the first few: ");

                int nPrinted = 0;
                foreach (var unmatchedRead in unmatchedReads[nextQueue].Select(_ => _.Value))
                {
                    outputFile.WriteLine(unmatchedRead[0].qname);
                    if (nPrinted++ >= 10)
                    {
                        break;
                    }
                } // foreach 
            } // next queue

            ASETools.WriteLineToMultipleStreams(outputStreams, "Processed " + ASETools.NumberWithCommas(nReads) + " total from " + nInputFiles + " input files, generating " + ASETools.NumberWithCommas(nMatches) + " matches");
            for (int i = 0; i < nInputFiles; i++)
            {
                ASETools.WriteLineToMultipleStreams(outputStreams, (char)('a' + i) + ": " + args[i+nextArg] + " (" + ASETools.NumberWithCommas(unmatchedReads[i].Count()) + " unmatched)");
            }

            if (wgsim)
            {
                ASETools.WriteLineToMultipleStreams(outputStreams, (char)('a' + nInputFiles) + ": correct answer from wgsim (only exact and normal matches)");
            }

            ASETools.WriteLineToMultipleStreams(outputStreams, "");
            WriteMatchClasses(outputStreams, "Matched alignment exactly:", matchedAlignmentExactly, nClasses);

            ASETools.WriteLineToMultipleStreams(outputStreams);
            ASETools.WriteLineToMultipleStreams(outputStreams, "");
            WriteMatchClasses(outputStreams, "Matched alignment:", matchedAlignment, nClasses);

            ASETools.WriteLineToMultipleStreams(outputStreams);
            ASETools.WriteLineToMultipleStreams(outputStreams, "");
            WriteMatchClasses(outputStreams, "Matched alignment w/MAPQ:", matchedAlignmentWithMAPQ, nInputFiles);   // NB: no wgsim here, so nInputFiles

            ASETools.WriteLineToMultipleStreams(outputStreams);
            ASETools.WriteLineToMultipleStreams(outputStreams, "");
            WriteMatchClasses(outputStreams, "Matched alignment w/MAPQ and cigar:", matchedAlignmentMAPQAndCigar, nInputFiles); // NB: no wgsim here, so nInputFiles

            outputStreams.ForEach(_ => _.Close());

            Console.WriteLine("Finished in " + ASETools.ElapsedTimeInSeconds(timer));

        } // Main
    } // Program
} // namespace
