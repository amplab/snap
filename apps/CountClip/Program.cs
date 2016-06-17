using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace CountClip
{
    class Program
    {
        static void Main(string[] args)
        {
            var reader = new StreamReader(args[0]);

            int nSoft = 0;
            int nHard = 0;
            int nBoth = 0;
            int nUnmapped = 0;
            int totalReads = 0;

            string read;
            while (null != (read = reader.ReadLine())) 
            {
                if (read[0] == '@') continue;
                var cigar = read.Split('\t')[5];
                int flags = Convert.ToInt32(read.Split('\t')[1]);

                bool foundS = false;
                bool foundH = false;
                bool unmapped = false;

                for (int i = 0; i < cigar.Count(); i++)
                {
                    if ((flags & 4) == 4) unmapped = true;
                    else if (cigar[i] == 'S') foundS = true;
                    else if (cigar[i] == 'H') foundH = true;
                }

                totalReads++;
                if (foundH && foundS) nBoth++;
                else if (foundS) nSoft++;
                else if (foundH) nHard++;
                else if (unmapped) nUnmapped++;

                if (!(foundS || foundH || unmapped))
                {
                    Console.WriteLine(read);
                }
            }

            Console.WriteLine("" + totalReads + " reads " + nSoft + " soft " + nHard + " hard " + nBoth + " both " + nUnmapped + " unmapped " + (totalReads - nSoft - nHard - nBoth - nUnmapped) + " unclipped.");
        }
    }
}
