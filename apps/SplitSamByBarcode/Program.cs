using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace SplitSamByBarcode
{
    class Program
    {


        static void Main(string[] args)
        {
            if (args.Count() != 3)
            {
                Console.WriteLine("usage: SplitSamByBarcode inputFile.sam outputFilePrefix groupSize");
                return;
            }
            StreamReader inputFile = new StreamReader(args[0]);

            var outputFiles = new Dictionary<string, StreamWriter>();
            var header = new List<string>();

            //
            // Read the header from the input file.  It's all the lines that start with a '@'
            //

            var inputLine = inputFile.ReadLine();
            int groupSize = Convert.ToInt32(args[2]);

            while (inputLine != null && inputLine != "" && inputLine[0] == '@')
            {
                header.Add(inputLine);
                inputLine = inputFile.ReadLine();
            }

            var seenChromosomes = new Dictionary<string, bool>();

            //
            // Now put each line into an output file based on its barcode.
            //
            int groupNumber = -1;
            int nRemainingInCurrentGroup = 0;
            long totalLines = 0;
            StreamWriter currentGroupStream = null;
            var allStreams = new List<StreamWriter>();
            Console.Write("Processing contig: ");
            while (null != inputLine)
            {
                totalLines++;
                var fields = inputLine.Split('\t');
                if (fields.Count() < 12)
                {
                    Console.WriteLine("Malformed input line or no optional fields in SAM file: " + inputLine);
                }
                else
                {
                    if (!seenChromosomes.ContainsKey(fields[2]))
                    {
                        Console.Write(fields[2] + " ");
                        seenChromosomes.Add(fields[2], true);
                    }

                    //
                    // Find the barcode.
                    //
                    string barcode = "BarcodeNotFound";
                    for (int i = 11; i < fields.Count(); i++)
                    {
                        var subfields = fields[i].Split(':');
                        if (subfields.Count() == 3 && subfields[0] == "RX" && subfields[1] == "Z")
                        {
                            barcode = subfields[2];
                            break;
                        }
                    }

                    if (!outputFiles.ContainsKey(barcode))
                    {
                        if (nRemainingInCurrentGroup > 0)
                        {
                            nRemainingInCurrentGroup--;
                            outputFiles.Add(barcode, currentGroupStream);
                        }
                        else
                        {
                            groupNumber++;
                            currentGroupStream = new StreamWriter(args[1] + "group-" + groupNumber + ".sam");
                            allStreams.Add(currentGroupStream);
                            outputFiles.Add(barcode, currentGroupStream);
                            foreach (var headerLine in header)
                            {
                                outputFiles[barcode].WriteLine(headerLine);
                            }
                            nRemainingInCurrentGroup = groupSize - 1;
                        }
                    }
                    outputFiles[barcode].WriteLine(inputLine);
                }

                inputLine = inputFile.ReadLine();
            }

            Console.WriteLine();
            Console.WriteLine("Processed " + outputFiles.Count() + " barcodes in " + (groupNumber + 1) + " output files among " + totalLines + " of input file.");

            foreach (var stream in allStreams) {
                stream.Close();
            }
        }
    }
}
