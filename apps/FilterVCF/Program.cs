using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;
using System.Diagnostics;

namespace FilterVCF
{
    class Program
    {
        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            if (args.Count() != 3)
            {
                Console.WriteLine("usage: filterVCF minQuality inputFile outputFile");
                return;
            }

            var minQuality = Convert.ToDouble(args[0]);

            var inputFile = ASETools.CreateStreamReaderWithRetry(args[1]);
            if (inputFile == null)
            {
                return;
            }

            var outputFile = ASETools.CreateStreamWriterWithRetry(args[2]);
            if (outputFile == null)
            {
                return;
            }

            string inputLine;
            int nHeader = 0;
            int nEmitted = 0;
            int nSkipped = 0;
            while (null != (inputLine = inputFile.ReadLine()))
            {
                if (inputLine.Length == 0 || inputLine[0] == '#')
                {
                    nHeader++;
                    outputFile.WriteLine(inputLine);
                    continue;
                }

                var fields = inputLine.Split('\t');
                if (fields.Count() < 6)
                {
                    throw new Exception("Too few fields in input line: " + inputLine);
                }

                if (Convert.ToDouble(fields[5]) < minQuality) {
                    nSkipped++;
                    continue;
                }

                nEmitted++;
                outputFile.WriteLine(inputLine);
            } // while we have input

            inputFile.Close();
            outputFile.Close();

            Console.WriteLine(nHeader + " header lines, " + nEmitted + " good enough lines and " + nSkipped + " too low quality in " + ASETools.ElapsedTimeInSeconds(timer));
        } // Main
    }
}
