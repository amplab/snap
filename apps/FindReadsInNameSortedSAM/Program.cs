using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;


namespace FindReadsInNameSortedSAM
{
    class Program
    {
        static void Main(string[] args)
        {
            if (args.Count() < 3)
            {
                Console.WriteLine("usage: input_sam_file output_sam_file {readName(s)}");
                return;
            }

            var inputStream = ASETools.CreateStreamReaderWithRetry(args[0]);
            if (null == inputStream)
            {
                Console.WriteLine("Unable to open " + args[0] + " for input.");
                return;
            }
 
            var outputFile = ASETools.CreateStreamWriterWithRetry(args[1]);
            if (null == outputFile)
            {
                Console.WriteLine("Unable to open " + args[1] + " for output");
                return;
            }

            //
            // Copy the header lines, which all start with "@" and come before any other lines.
            //
            string line;
            while (null != (line = inputStream.ReadLine()) && line.StartsWith("@"))
            {
                outputFile.WriteLine(line);
            }

            inputStream.Close();

            FileStream inputFile;
            try
            {
                inputFile = File.OpenRead(args[0]);
            }
            catch (Exception e)
            {
                Console.WriteLine("Failed to open input file " + args[0] + ", exception text: " + e.Message);
                return;
            }


            for (int whichRead = 2; whichRead < args.Count(); whichRead++)
            {
                var qnameToMatch = args[whichRead];

                long min = 0;
                long max = inputFile.Length;

                while (true)
                {
                    long probe = (min + max) / 2;

                    line = ASETools.ReadNextLineStartingAt(inputFile, probe);

                    bool tooSmall;

                    if (line == null)
                    {
                        //
                        // EOF
                        //
                        tooSmall = false;
                    } 
                    else if (line.StartsWith("@"))
                    {
                        //
                        // header line is less than any read
                        //
                        tooSmall = true;
                    } 
                    else
                    {
                        var samLine = new ASETools.SAMLine(line);

                        if (samLine.qname != qnameToMatch)
                        {
                            tooSmall = samLine.CompareTo(qnameToMatch) < 0;
                        }
                        else
                        {
                            //
                            // We don't know if this is the first occurrence of this line.  Start backing up until we find
                            // one that's not it (or hit beginning of file).
                            //
                            long currentOffset = Math.Max(probe - 1024, 0);

                            while (currentOffset > 0)
                            {
                                line = ASETools.ReadNextLineStartingAt(inputFile, currentOffset);
                                if (line == null)
                                {
                                    Console.WriteLine("Unexpected EOF???");
                                    return;
                                }

                                if (line.StartsWith("@"))
                                {
                                    break;
                                }

                                samLine = new ASETools.SAMLine(line);
                                if (samLine.qname != qnameToMatch)
                                {
                                    break;
                                }

                                currentOffset = Math.Max(currentOffset - 1024, 0);
                            } // while we're not at the begining of the file and haven't found a line that doesn't match the qnameToMatch

                            //
                            // Now read until we find a line that does match, since we might have backed up too much.
                            //
                            while (line.StartsWith("@") || (samLine = new ASETools.SAMLine(line)).qname != qnameToMatch)
                            {
                                currentOffset += line.Length + 1;   // ??? is this right for crlf text?
                                line = ASETools.ReadNextLineStartingAt(inputFile, currentOffset);
                            }

                            //
                            // Now read and copy out all matching lines.
                            //
                            while (samLine.qname == qnameToMatch )
                            {
                                outputFile.WriteLine(line);
                                Console.WriteLine(line);

                                currentOffset += line.Length + 1;   // ??? is this right for crlf text?
                                line = ASETools.ReadNextLineStartingAt(inputFile, currentOffset);

                                if (line == null)
                                {
                                    break;
                                }

                                samLine = new ASETools.SAMLine(line);
                            } // while the qnames match

                            break; // out of the min/max loop (i.e., we're done with this read name)
                        } // matching qname
                    } // read a good line

                    if (min >= max)
                    {
                        break;
                    }

                    if (max - min == 1)
                    {
                        //
                        // This handles the case where they stop converging because of roundoff in integer division.
                        //
                        if (tooSmall)
                        {
                            min = max;
                        } 
                        else
                        {
                            max = min;
                        }
                    }
                    else if (tooSmall)
                    {
                        min = probe;
                    } else
                    {
                        max = probe;
                    }
                } // while min <= max
            } // for each qname to find

            outputFile.Close();
            inputFile.Close();

        } // Main
    } // Program
} // namespace
