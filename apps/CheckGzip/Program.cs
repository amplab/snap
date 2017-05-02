using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;

namespace CheckGzip
{
    class Program
    {
        static void Main(string[] args)
        {
            if (args.Count() == 0 || args.Count() > 3)
            {
                Console.WriteLine("usage: CheckGZip filename {firstLine|* {lastLine}}");
                return;
            }

            string firstLine = "*";
            string lastLine = "*";

            if (args.Count() > 1)
            {
                firstLine = args[1];

                if (args.Count() > 2)
                {
                    lastLine = args[2];
                }
            }

            try
            {
                var reader = ASETools.CreateCompressedStreamReaderWithRetry(args[0]);

                if (reader == null)
                {
                    Console.WriteLine(args[0]);
                    return;
                }

                var line = reader.ReadLine();

                if (firstLine == "*" && lastLine == "*")
                {
                    return;
                }

                if (line != firstLine && firstLine != "*")
                {
                    Console.WriteLine(args[0]);
                    return;
                }

                if (lastLine == "*")
                {
                    return;
                }

                while (true)
                {
                    var nextLine = reader.ReadLine();
                    if (nextLine == null)
                    {
                        break;
                    }

                    line = nextLine;
                }

                if (line != lastLine)
                {
                    Console.WriteLine(args[0]);
                }
            }
            catch
            {
                Console.WriteLine(args[0]);
            }
        }
    }
}
