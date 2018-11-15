using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.Diagnostics;
using System.IO;
using System.Threading;

namespace RemoveEmptyDirectories
{
    class Program
    {

        static long nEmptyDirectoriesRemoved = 0;
        static ASETools.Configuration configuration;

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            configuration = ASETools.Configuration.loadFromFile(args);

            if (null == configuration)
            {
                Console.WriteLine("Giving up because we were unable to load configuration.");
                return;
            }

            int nPerDot;

            ASETools.PrintMessageAndNumberBar("Processing", "data directories", configuration.dataDirectories.Count(), out nPerDot);

            var threading = new ASETools.WorkerThreadHelper<string, int>(configuration.dataDirectories, HandleOneDataDirectory, null, null, nPerDot);
            threading.run();

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(timer));

            Console.WriteLine("Removed " + nEmptyDirectoriesRemoved + " empty directories");
        } // Main

        static void HandleOneDataDirectory(string dataDirectory, int state)
        {
            string[] dirsToCheck = { dataDirectory, dataDirectory + @"..\" + configuration.derivedFilesDirectory };

            int nEmpty = 0;

            foreach (var dirToCheck in dirsToCheck)
            {
                try
                {
                    foreach (var subdir in Directory.EnumerateDirectories(dirToCheck))
                    {
                        bool anyFound = false;
                        foreach (var containedFile in Directory.EnumerateFiles(subdir))
                        {
                            anyFound = true;
                            break;
                        }

                        if (!anyFound)
                        {
                            foreach (var subsubDir in Directory.EnumerateDirectories(subdir))
                            {
                                anyFound = true;
                            }
                        }

                        if (!anyFound)
                        {
                            Directory.Delete(subdir);   // This version won't delete directories with contents, so it can't screw up in the awful way.
                            nEmpty++;
                        }
                    } // each subdir
                } catch
                {
                    Console.Write("Exception processing data directory " + dirToCheck);
                }
            } // downloaded or derived

            Interlocked.Add(ref nEmptyDirectoriesRemoved, nEmpty);
        } // HandleOneDataDirectory
    }
}
