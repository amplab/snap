using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;

namespace DumpSubfile
{
    class Program
    {
        static void Main(string[] args)
        {
            if (args.Count() != 2 && args.Count() != 3)
            {
                Console.WriteLine("usage: DumpSubfile compositeFileName subfileName {outputFile}");
                Console.WriteLine("If you omit outputfile, then the subfile will be written to stdout.");
                return;
            }

            var reader = new ASETools.ConsolodatedFileReader();

            if (!reader.open(args[0]))
            {
                Console.WriteLine("Unable to open containing file " + args[0]);
                return;
            }

            var subfileReader = reader.getSubfile(args[1]);
            if (null == subfileReader)
            {
                Console.WriteLine("Unable to open subfile " + args[1]);
                return;
            }

            StreamWriter writer = null;
            if (args.Count() > 2)
            {
                writer = ASETools.CreateStreamWriterWithRetry(args[2]);

                if (null == writer)
                {
                    Console.WriteLine("Failed to open output file " + args[2]);
                    return;
                }
            }

            string line;
            while (null != (line = subfileReader.ReadLine()))
            {
                if (null == writer)
                {
                    Console.WriteLine(line);
                } else
                {
                    writer.WriteLine(line);
                }
            }

            if (null != writer)
            {
                writer.Close();
            }


        }
    }
}
