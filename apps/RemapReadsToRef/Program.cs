using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;

//
// Take a set of long reads drawn from a sample and mapped to the reference and a set of short reads
// mapped to the long reads and remap the short reads to the original reference.
//

namespace RemapReadsToRef
{
    internal class Program
    {
        static void Main(string[] args)
        {
            if (args.Length != 3)
            {
                Console.WriteLine("usage: RemapReadsToRef longReads.sam shortReadFilenameTemplate outputFile.sam");
                Console.WriteLine("shortReadFilenameTemplate will have _n.fastq added for n from 0 to whatever fails first.");
                return;
            }

            var longreadInputFilename = args[0];
            var shortReadFilenameTemplate = args[1];
            var outputFilename = args[2];

            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
            if (outputFile == null)
            {
                Console.WriteLine("Unable to open " + outputFilename + " for write.");
                return;
            }

            var longreadInputFile = ASETools.CreateStreamReaderWithRetry(longreadInputFilename);
            if (null == longreadInputFile)
            {
                Console.WriteLine("unable to open " + longreadInputFile + " for read");
                return;
            }


            var longReads = ASETools.SAMLine.ReadFromFile(longreadInputFile);
            var longReadsByName = new Dictionary<string, ASETools.SAMLine>();
            longReads.ForEach(_ => longReadsByName.Add(_.qname, _));    // Wonder if these are all unique....

        } // Main
    } // Program
} // namespace
