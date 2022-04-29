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

namespace CompareAlignments
{
    internal class Program
    {
        interface ReadGetter
        {
            ASETools.SAMLine GetRead();
        }

        class StreamReaderReadGetter : ReadGetter
        {
            public StreamReaderReadGetter(StreamReader reader_)
            {
                reader = reader_;
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

                    var samLine = new ASETools.SAMLine(line);
                    if (samLine.isSupplementaryAlignment() || samLine.isSecondaryAlignment())
                    {
                        continue;   // ignore these
                    }
                    queue.Enqueue(samLine);
                }

                queue.TerminateWriter();
            }

            public ASETools.SAMLine GetRead()
            {
                ASETools.SAMLine line;
                if (!queue.Dequeue(out line))
                {
                    return null;
                }

                return line;
            } // GetRead

            ASETools.ThrottledParallelQueue<ASETools.SAMLine> queue = new ASETools.ThrottledParallelQueue<ASETools.SAMLine>(20000, 1);

            StreamReader reader;
        } // StreamReaderReadGetter

        class BAMReadGetter : ReadGetter
        {
            public BAMReadGetter(string inputFilename_)
            {
                inputFilename = inputFilename_;
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
                var line = new ASETools.SAMLine(rawLine);

                if (line.isSecondaryAlignment() || line.isSupplementaryAlignment()) // Drop secondary and supplementary alignments.
                {
                    return;
                }
                queue.Enqueue(line);
            }

            public ASETools.SAMLine GetRead()
            {
                ASETools.SAMLine line;
                if (!queue.Dequeue(out line))
                {
                    return null;
                }

                return line;
            } // GetRead


            ASETools.ThrottledParallelQueue<ASETools.SAMLine> queue = new ASETools.ThrottledParallelQueue<ASETools.SAMLine>(20000, 1);
            string inputFilename;
        }


        delegate bool MatchFunction(int a, int b);  // This must be equivalence classes.  i.e., matchFunction(a, b) && matchFunction(b, c) => matchFunction(a, c)
        static string partitionString(int nInputFiles, MatchFunction matchFunction)
        {
            var result = "";

            var used = new bool[nInputFiles];
            for (int i = 0; i < nInputFiles; i++)
            {
                used[i] = false;
            }

            bool needBar = false;
            for (int i = 0; i < nInputFiles; i++)
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

                for (int j = i+1; j < nInputFiles; j++)
                {
                    if (!used[j] && matchFunction(j, i)) // !used is unnecessary and here just to avoid calling matchFunction when we don't have to
                    {
                        result += (char)('a' + j);
                        used[j] = true;
                    }
                } // for each potential partner
            } // for each input file

            return result;
        } // paritionString

        delegate bool SAMComparitor(ASETools.SAMLine read0, ASETools.SAMLine read1);

        static bool doesAlignmentMatchExactly(ASETools.SAMLine read0, ASETools.SAMLine read1)
        {
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
        static bool doesAlignmentMatch(ASETools.SAMLine read0, ASETools.SAMLine read1)
        {
            const int fuzziness = 50;   // How close is close enough (except that to make it equivalence classes we just divide by fuzziness, so fuzziness - 2 and fuzziness + 2 don't match)

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

        static bool doesAlignmentWithMAPQMatch(ASETools.SAMLine read0, ASETools.SAMLine read1)
        {
            return read0.isUnmappedOrMapq0() && read1.isUnmappedOrMapq0() || !read0.isUnmappedOrMapq0() && !read1.isUnmappedOrMapq0() && doesAlignmentMatch(read0, read1);
        }

        static bool doesAlignmentMAPQAndCigarMatch(ASETools.SAMLine read0, ASETools.SAMLine read1)
        {
            return doesAlignmentWithMAPQMatch(read0, read1) && (read0.mapq == read1.mapq || read0.cigar == read1.cigar);
        }

        static bool compareReads(Dictionary<string, List<ASETools.SAMLine>> [] unmatchedReads, ASETools.SAMLine example, int exampleFrom, int a, int b, SAMComparitor comparitor)
        {
            ASETools.SAMLine read0, read1;
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

        static void WriteMatchClasses(List<StreamWriter> outputStreams, string header, Dictionary<string, long> matches)
        {
            ASETools.WriteLineToMultipleStreams(outputStreams, header);

            var total = matches.Select(_ => _.Value).Sum();

            var list = matches.ToList();
            list.Sort((a, b) => b.Value.CompareTo(a.Value));    // Backward comparison to pu the big ones on top

            foreach (var match in list)
            {
                ASETools.WriteLineToMultipleStreams(outputStreams, match.Key + "\t" + ASETools.NumberWithCommas(match.Value) + "\t" + ASETools.Percentage(match.Value, total, 2));
            }
        } // WriteMatchClasses

        public static int nInputFiles;

        class MatchGroup
        {
            public MatchGroup(SAMComparitor comparitor_, Dictionary<string,long> matched_, ASETools.ThrottledParallelQueue<List<ASETools.SAMLine>> inputQueue_)
            {
                comparitor = comparitor_;
                matched = matched_;
                inputQueue = inputQueue_;

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

            static bool compareReads(List<ASETools.SAMLine> reads, int a, int b, SAMComparitor comparitor)
            {
                return comparitor(reads[a], reads[b]);
            }

            void Worker()
            {
                List<ASETools.SAMLine> reads;
                while (inputQueue.Dequeue(out reads))
                {
                    var partition = partitionString(nInputFiles, (a, b) => compareReads(reads, a, b, comparitor));
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
            ASETools.ThrottledParallelQueue<List<ASETools.SAMLine>> inputQueue;
            Dictionary<string, long> matched;
        } // MatchGroup
        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            if (args.Length < 3 || args[0].ToLower().EndsWith(".sam") || args[0].ToLower().EndsWith(".bam"))
            {
                Console.WriteLine("usage: CompareAligments outputFile inputFile1 inputFile2 <...inputFileN>");
                Console.WriteLine("output file cannot end with .sam or .bam (to keep you from accidentally ovewriting an input)");
                return;
            }

            var outputFile = ASETools.CreateStreamWriterWithRetry(args[0]);
            if (outputFile == null)
            {
                Console.WriteLine("Unable to open output file " + args[0]);
                return;
            }

            var outputStreams = new List<StreamWriter>();
            outputStreams.Add(new StreamWriter(Console.OpenStandardOutput()));
            outputStreams.Add(outputFile);

            nInputFiles = args.Length - 1;
            var unmatchedReads = new Dictionary<string, List<ASETools.SAMLine>>[nInputFiles];
            var queueDone = new bool[nInputFiles];

            var readerQueues = new ReadGetter[nInputFiles];

            for (int i = 0; i < nInputFiles; i++)
            {
                var inputFilename = args[i+1]; // +1 skips over output file
                if (inputFilename.ToLower().EndsWith(".bam"))
                {
                    readerQueues[i] = new BAMReadGetter(inputFilename);
                }
                else
                {
                    var inputFile = ASETools.CreateStreamReaderWithRetry(inputFilename);
                    if (null == inputFile)
                    {
                        Console.WriteLine("Unable to open " + inputFilename);
                        return;
                    }

                    readerQueues[i] = new StreamReaderReadGetter(inputFile);
                }

                unmatchedReads[i] = new Dictionary<string, List<ASETools.SAMLine>>();
                queueDone[i] = false;
            } // start reader threads

            //
            // These map partition strings to counts.  They're sparse because the space of possible partitons grows hyperexponentially in nInputFiles, so most won't be used
            // if we have many aligners.
            //
            var matchedAlignmentExactly = new Dictionary<string, long>();
            var matchedAlignment = new Dictionary<string, long>();
            var matchedAlignmentWithMAPQ = new Dictionary<string, long>();
            var matchedAlignmentMAPQAndCigar = new Dictionary<string, long>();

            var matchGroups = new List<MatchGroup>();
            var queues = new List<ASETools.ThrottledParallelQueue<List<ASETools.SAMLine>>>();
 
            { // scope to make the definition of queue be local
                var queue = new ASETools.ThrottledParallelQueue<List<ASETools.SAMLine>>(5000, 1);
                queues.Add(queue);
                matchGroups.Add(new MatchGroup(doesAlignmentMatchExactly, matchedAlignmentExactly, queue));

                queue = new ASETools.ThrottledParallelQueue<List<ASETools.SAMLine>>(5000, 1);
                queues.Add(queue);
                matchGroups.Add(new MatchGroup(doesAlignmentMatch, matchedAlignment, queue));

                queue = new ASETools.ThrottledParallelQueue<List<ASETools.SAMLine>>(5000, 1);
                queues.Add(queue);
                matchGroups.Add(new MatchGroup(doesAlignmentWithMAPQMatch, matchedAlignmentWithMAPQ, queue));

                queue = new ASETools.ThrottledParallelQueue<List<ASETools.SAMLine>>(5000, 1);
                queues.Add(queue);
                matchGroups.Add(new MatchGroup(doesAlignmentMAPQAndCigarMatch, matchedAlignmentMAPQAndCigar, queue));
            }


            int nextQueue = 0;
            long nReads = 0;
            long nMatches = 0;

            while (!queueDone.All(_ => _))
            {
                if (!queueDone[nextQueue])
                {
                    ASETools.SAMLine samLine;
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

                        var reads = new List<ASETools.SAMLine>();
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
                            unmatchedReads[nextQueue].Add(samLine.qname, new List<ASETools.SAMLine>());
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
                ASETools.WriteLineToMultipleStreams(outputStreams, (char)('a' + i) + ": " + args[i+1] + " (" + ASETools.NumberWithCommas(unmatchedReads[i].Count()) + " unmatched)");
            }

            ASETools.WriteLineToMultipleStreams(outputStreams, "");
            WriteMatchClasses(outputStreams, "Matched alignment exactly:", matchedAlignmentExactly);

            ASETools.WriteLineToMultipleStreams(outputStreams, "");
            WriteMatchClasses(outputStreams, "Matched alignment:", matchedAlignment);

            ASETools.WriteLineToMultipleStreams(outputStreams, "");
            WriteMatchClasses(outputStreams, "Matched alignment w/MAPQ:", matchedAlignmentWithMAPQ);

            ASETools.WriteLineToMultipleStreams(outputStreams, "");
            WriteMatchClasses(outputStreams, "Matched alignment w/MAPQ and cigar:", matchedAlignmentMAPQAndCigar);

            outputStreams.ForEach(_ => _.Close());

            Console.WriteLine("Finished in " + ASETools.ElapsedTimeInSeconds(timer));

        } // Main
    } // Program
} // namespace
