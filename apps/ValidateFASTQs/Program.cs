using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.Threading;
using System.IO;

namespace ValidateFASTQs
{
    class Program
    {
        static ASETools.CommonData commonData = null;
        class SecondEnd { }

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);
            if (null == commonData)
            {
                return;
            }

            if (commonData.configuration.commandLineArgs.Count() == 0)
            {
                Console.WriteLine("usage: ValidateVCFs {caseId...caseId...| -a}");
            }

            bool doAllCases = commonData.configuration.commandLineArgs[0] == "-a";

            if (doAllCases)
            {
                Console.WriteLine("Writeme");
                return;
            }

            List<ASETools.Case> casesToDo;
            if (doAllCases)
            {
                casesToDo = commonData.listOfCases.Where(_ => _.normal_dna_fastq_filename != "").ToList();
            }
            else
            {
                casesToDo = new List<ASETools.Case>();
                foreach (var caseId in commonData.configuration.commandLineArgs)
                {

                    if (!commonData.cases.ContainsKey(caseId))
                    {
                        Console.WriteLine("Couldn't find case ID " + caseId);
                        return;
                    }
                    casesToDo.Add(commonData.cases[caseId]);
                }
            }

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "cases", casesToDo.Count(), out nPerDot);

            var workerThreads = new ASETools.WorkerThreadHelper<ASETools.Case, int>(casesToDo, HandleOneCase, null, null, nPerDot);
            workerThreads.run();

            Console.Write(ASETools.ElapsedTimeInSeconds(commonData.timer));

        } // Main

        static void HandleOneCase(ASETools.Case case_, int state)
        {
            if (case_.normal_dna_fastq_filename == "")
            {
                Console.WriteLine("Case " + case_.case_id + " doesn't have a normal DNA FASTQ to check.");
                return;
            }

            if (case_.normal_dna_fastq_second_end == "")
            {
                HandleEnd(case_.normal_dna_fastq_filename, null);
            } else
            {
                ASETools.ThrottledParallelQueue<string[]>[] queues = new ASETools.ThrottledParallelQueue<string[]>[2];
                var threads = new Thread[2];

                for (int i = 0; i < 2; i++)
                {
                    queues[i] = new ASETools.ThrottledParallelQueue<string[]>(10000, 1);
                    var filename = i == 0 ? case_.normal_dna_fastq_filename : case_.normal_dna_fastq_second_end;
                    var queue = queues[i];
                    threads[i] = new Thread(() => HandleEnd(filename, queue));
                    threads[i].Start();
                }

                string[] read0, read1;
                while (queues[0].Dequeue(out read0))
                {
                    if (!queues[1].Dequeue(out read1))
                    {
                        FailOneCase("Premature EOF on " + case_.normal_dna_fastq_second_end, queues);
                        return;
                    }

                    if (!read0[0].EndsWith("/1"))
                    {
                        FailOneCase(case_.normal_dna_fastq_filename + ": Read 1 ID doesn't end with /1", queues);
                        return;
                    }

                    if (!read1[0].EndsWith("/2"))
                    {
                        FailOneCase(case_.normal_dna_fastq_second_end + ": Read 2 ID doesn't end with /2", queues);
                        return;
                    }

                    if (read0[0].Substring(0, read0[0].Length - 2) != read1[0].Substring(0, read1[0].Length - 2))
                    {
                        FailOneCase("Read ID mismatch for case " + case_.case_id + ": " + read0[0] + " and " + read1[0], queues);
                        return;
                    }
                } // While queues[0] isn't empty

                if (queues[1].Dequeue(out read1))
                {
                    FailOneCase("Premature EOF on " + case_.normal_dna_fastq_filename, queues);
                }

                foreach (var thread in threads)
                {
                    thread.Join();
                }
            } // if it's paired-end

        } // HandleOneCase

        static void FailOneCase(string errorMessage, ASETools.ThrottledParallelQueue<string[]>[] queues)
        {
            Console.WriteLine(errorMessage);
            foreach (var queue in queues)
            {
                string[] value;
                while (queue.Dequeue(out value))
                {
                    // Nothing.
                }
            }
        } // Drain Queues



        static void HandleEnd(string inputFilename, ASETools.ThrottledParallelQueue<string[]> outputQueue)
        {
            var inputFile = ASETools.CreateStreamReaderWithRetry(inputFilename);
            if (null == inputFile)
            {
                Console.WriteLine("Unable to open " + inputFilename);
                if (outputQueue != null)
                {
                    outputQueue.TerminateWriter();
                    return;
                }
            }

            string inputLine;
            while (null != (inputLine = inputFile.ReadLine()))
            {
                string[] inputLines = new string[4];
                inputLines[0] = inputLine;

                for (int i = 1; i < 4; i++)
                {
                    inputLines[i] = inputFile.ReadLine();
                    if (inputLines[i] == null)
                    {
                        FailHandleEnd("EOF in mid-read on " + inputFilename, inputFile, outputQueue);
                        return;
                    }
                }

                if (inputLines[0].Length == 0 || inputLines[0][0] != '@')
                {
                    FailHandleEnd(inputFilename + " has a read that has an ID line that doesn't start with @: " + inputLines[0], inputFile, outputQueue);
                    return;
                }

                if (inputLines[1].Length == 0)
                {
                    Console.WriteLine(inputFilename + " has a read with no bases", inputFile, outputQueue);
                    return;
                }

                if (inputLines[2] != "+")
                {
                    //Technically +ReadID is permitted, but bedtools never generates it.
                    FailHandleEnd(inputFilename + " has read with third line not +: " + inputLines[2], inputFile, outputQueue);
                    return;
                }

                if (inputLines[3].Length != inputLines[1].Length)
                {
                    FailHandleEnd(inputFilename + " has mismatched bases and quality: " + inputLines[1] + " and " + inputLines[3], inputFile, outputQueue);
                    return;
                }

                if (null != outputQueue)
                {
                    outputQueue.Enqueue(inputLines);
                }

            }

            outputQueue.TerminateWriter();
            inputFile.Close();
        } // HandleEnd

        static void FailHandleEnd(string errorMessage, StreamReader inputFile, ASETools.ThrottledParallelQueue<string[]> outputQueue)
        {
            Console.WriteLine(errorMessage);
            inputFile.Close();
            if (null != outputQueue)
            {
                outputQueue.TerminateWriter();
            }
        } // FailHandleEnd
    }
}
