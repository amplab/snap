using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.ComponentModel;

namespace ExtractReadFromSAMFile
{
    class Program
    {
        static void Main(string[] args)
        {
            if (args.Count() != 3)
            {
                Console.WriteLine("usage: ExtractReadFromSAMFile input.sam output.sam readID");
                return;
            }

            var inputFile = ASETools.CreateStreamReaderWithRetry(args[0]);
            var outputFile = ASETools.CreateStreamWriterWithRetry(args[1]);

            var idAndTab = args[2] + "\t";

            int nReadsFound = 0;

            string line;
            while (null != (line = inputFile.ReadLine()))
            {
                if (line.StartsWith("@"))
                {
                    outputFile.WriteLine(line);
                    continue;
                }

                if (line.StartsWith(idAndTab))
                {
                    outputFile.WriteLine(line);
                    Console.WriteLine(line);

                    var samLine = new ASETools.SAMLine(line);
                    if (!samLine.isPaired() || nReadsFound == 1)
                    {
                        outputFile.Close();
                        return;
                    }

                    nReadsFound++;
                } // line starts with id\t
            } // while we have an input line

            outputFile.Close();
            File.Delete(args[1]);

            Console.WriteLine("Never found read");
        } // Main
    } // Program
} // namespace
