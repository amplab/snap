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
                var lines  = ASETools.FindSAMLinesInNameSortedSAMFile(inputFile, args[whichRead]);
                foreach (var samLine in lines)
                {
                    outputFile.WriteLine(samLine);
                    Console.WriteLine(samLine);
                }
            } // for each qname to find

            outputFile.Close();
            inputFile.Close();

        } // Main
    } // Program
} // namespace
