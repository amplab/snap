using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;

namespace CountMappedBases
{
    class Program
    {
        static int Main(string[] args)
        {
            if (args.Count() != 2)
            {
                Console.WriteLine("usage: CountMappedBases inputAllcountFilename outputFilename");
                return -1;
            }

            var reader = new ASETools.AllcountReader(args[0]);

            if (!reader.openFile())
            {
                Console.WriteLine("unable to open or corrupt input file :" + args[0]);
                return -1;
            }

            long totalMappedBaseCount = 0;

            if (!reader.ReadAllcountFile((contigName, location, mappedCount) => totalMappedBaseCount += mappedCount))
            {
                Console.WriteLine("Error reading file " + args[0]);
                return -1;
            }

            var writer = ASETools.CreateStreamWriterWithRetry(args[1]);
            if (writer == null)
            {
                Console.WriteLine("Unable to open output file " + args[1]);
                return -1;
            }

            writer.WriteLine(totalMappedBaseCount + "\t" + args[0]);
            writer.WriteLine("**done**");
            writer.Close();

            return 0;
        }
    }
}
