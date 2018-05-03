using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;
using System.Diagnostics;

namespace CountMappedBases
{
    class Program
    {
        class FilenamePair
        {
            public FilenamePair(string inputFilename_, string outputFilename_)
            {
                inputFilename = inputFilename_;
                outputFilename = outputFilename_;
            }

            public readonly string inputFilename;
            public readonly string outputFilename;
        }

        static void HandleOneInputFile(FilenamePair filenamePair, int state)
        {
            var reader = new ASETools.AllcountReader(filenamePair.inputFilename);

            if (!reader.openFile())
            {
                Console.WriteLine("unable to open or corrupt input file :" + filenamePair.inputFilename);
                return;
            }

            long totalMappedBaseCount = 0;
            long basesCovered = 0;

            if (!reader.ReadAllcountFile((contigName, location, mappedCount) => { totalMappedBaseCount += mappedCount; basesCovered++; }))
            {
                Console.WriteLine("Error reading file " + filenamePair.inputFilename);
                return;
            }

            var writer = ASETools.CreateStreamWriterWithRetry(filenamePair.outputFilename);
            if (writer == null)
            {
                Console.WriteLine("Unable to open output file " + filenamePair.outputFilename);
                return;
            }

            writer.WriteLine(totalMappedBaseCount + "\t" + filenamePair.outputFilename + "\t" + basesCovered);
            writer.WriteLine("**done**");
            writer.Close();
        }

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            var configuration = ASETools.Configuration.loadFromFile(args);
            if (null == configuration)
            {
                Console.WriteLine("Unable to load configuration");
                return;
            }

            if (configuration.commandLineArgs.Count() < 2 || configuration.commandLineArgs.Count() % 2 != 0)
            {
                Console.WriteLine("usage: CountMappedBases {inputAllcountFilename outputFilename}");
                return;
            }

            var filesToProcess = new List<FilenamePair>();

            for (int i = 0; i < configuration.commandLineArgs.Count(); i += 2)
            {
                filesToProcess.Add(new FilenamePair(configuration.commandLineArgs[i], configuration.commandLineArgs[i + 1]));
            }

            var threading = new ASETools.WorkerThreadHelper<FilenamePair, int>(filesToProcess, HandleOneInputFile, null, null, 1);

            Console.Write("Processing " + filesToProcess.Count() + " input files, 1 dot/file: ");
            threading.run();

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(timer));
        }
    }
}
