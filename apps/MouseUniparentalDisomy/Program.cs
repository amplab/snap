using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;

namespace MouseUniparentalDisomy
{
    class Program
    {

        enum LocusTypes { NoData, Neither, B, Both, C};
        static void Main(string[] args)
        {

            Dictionary<string, Dictionary<string, LocusTypes>> results = new Dictionary<string, Dictionary<string, LocusTypes>>();  // cells -> locus -> type

            var inputFile = ASETools.CreateStreamReaderWithRetry(@"f:\temp\AustMms1stRdTest_w2s.txt");

            var headerLine = inputFile.ReadLine();
            var cells = headerLine.Split('\t').ToList();
            cells.RemoveAt(0);  // The first column is blank.

            cells.ForEach(_ => results.Add(_, new Dictionary<string, LocusTypes>()));

            string line;
            while (null != (line = (inputFile.ReadLine())))
            {
                var results = line.Split('\t');

                cells.ForEach(cell => results[cell].Add())
            }

            inputFile.Close();

        }
    }
}
