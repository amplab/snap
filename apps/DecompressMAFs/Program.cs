using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.Diagnostics;
using System.IO;

namespace DecompressMAFs
{
    class Program
    {
        static ASETools.CommonData commonData;

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);
            if (null == commonData)
            {
                return;
            }

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "diseases", commonData.diseases.Count(), out nPerDot);

            var threading = new ASETools.WorkerThreadHelper<string, int>(commonData.diseases, HandleOneDisease, null, null, nPerDot);
            threading.run();

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(commonData.timer));
        }

        static void HandleOneDisease(string disease, int state)
        {
            var cases = commonData.listOfCases.Where(_ => _.disease() == disease && _.maf_filename != "").ToList();

            if (cases.Count() == 0)
            {
                Console.WriteLine("No cases for disease " + disease + " have a MAF file.");
                return;
            }

            var bin = commonData.configuration.binariesDirectory;
            ASETools.RunAndWaitForProcess(bin + "cmd.exe", "/c " + bin + "cat.exe " + cases[0].maf_filename + " | " + bin + "zcat.exe -d > " + ASETools.GetDirectoryFromPathname(cases[0].maf_filename) + @"\" + disease + ".maf");
        } // HandleOneDisease
    }
}
