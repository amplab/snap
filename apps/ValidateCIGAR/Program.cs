using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.Diagnostics;
using System.Threading;

namespace ValidateCIGAR
{
    class Program
    {
        static ASETools.ThrottledParallelQueue<List<string>> dataQueue;
        static long nReads = 0;

        static void ReaderThread()
        {
            List<string> batch;

            while (dataQueue.Dequeue(out batch))
            {
                foreach (var readLine in batch)
                {
                    var read = new ASETools.SAMLine(readLine);

                    if (read.cigar != "*" && !read.cigar.Contains("M") && !read.cigar.Contains("I") && !read.cigar.Contains("D") && !read.cigar.Contains("N") && !read.cigar.Contains("="))
                    {
                        lock (dataQueue)
                        {
                            Console.WriteLine();
                            Console.WriteLine("Read's cigar string doesn't contain any of MIDN=: " + readLine);
                            Console.WriteLine();
                        }
                    } // if it's a bogus CIGAR
                } // foreach readline in the batch

                lock (dataQueue)
                {
                    if (nReads / 10000000 != (nReads + batch.Count()) / 10000000)
                    {
                        Console.Write(".");
                    }

                    nReads += batch.Count();
                } // lock
            } // while the queue isn't terminated
        } // Reader Thread
        

        static void Main(string[] args)
        {
            dataQueue = new ASETools.ThrottledParallelQueue<List<string>>(100, 1);

            var readers = new List<Thread>();
            for (int i = 0; i < 15; i++)
            {
                var thread = new Thread(() => ReaderThread());
                thread.Start();
                readers.Add(thread);
            } // thread start loop


            if (args.Count() != 1 || args[0] == "-?")
            {
                Console.WriteLine("usage: ValidateCIGAR inputBAMFile");
                return;
            }

            var startInfo = new ProcessStartInfo(@"c:\bolosky\bin\samtools.exe", "view " + args[0]);

            startInfo.RedirectStandardOutput = true;
            startInfo.UseShellExecute = false;

            var process = Process.Start(startInfo);

            Console.Write("Processing reads, one dot/10M: ");

            int nInBatch = 0;
            var batch = new List<string>();
            string inputLine;

            while (null != (inputLine = process.StandardOutput.ReadLine()))
            {
                batch.Add(inputLine);
                nInBatch++;
                if (nInBatch > 1000)
                {
                    lock (dataQueue)
                    {
                        dataQueue.Enqueue(batch);
                        batch = new List<string>();
                        nInBatch = 0;
                    } // lock
                } // if we have a full batch
            } // while we have input lines

            if (nInBatch != 0)
            {
                dataQueue.Enqueue(batch);
            }

            dataQueue.TerminateWriter();
            readers.ForEach(_ => _.Join()); // Wait for all the readers to terminate

            Console.WriteLine();
            Console.WriteLine("Processed " + nReads + " reads.");

        } // Main
    } // Program
} // namespace
