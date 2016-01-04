using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace ColumnizeSam
{
    class Program
    {
        static void Main(string[] args)
        {
            if (args.Count() < 2 || args.Count() > 3)
            {
                Console.WriteLine("usage: ColumnizeSam inputFile offset {chr}");
                return;
            }

            int target = Convert.ToInt32(args[1]);
            string chr = null;
            if (args.Count() == 3)
            {
                chr = args[2];
            }

            int lowestEmitted = -1;
            StreamReader inputFile = new StreamReader(args[0]);
            string line;
            while (null != (line = inputFile.ReadLine())) {
                var fields = line.Split('\t');
                if (fields.Count() < 10) {
                    Console.WriteLine("Ill-formed SAM line " + line);
                    continue;
                }

                int pos = Convert.ToInt32(fields[3]);
                int len = fields[9].Count();
                int flags = Convert.ToInt32(fields[1]);
                string seq = fields[9];
                string cigar = fields[5];

                if (pos + len <= target || pos > target || chr != null && fields[2] != chr || (flags & 0x4) != 0) { // The flags bit is checking for unmapped
                    continue;
                }

                string mate = "";
                if ((flags & 0x1) == 1)
                {
                    int tlen = Convert.ToInt32(fields[8]);
                    if (tlen < 0)
                    {
                        mate = "<";
                    }
                    else if (tlen > 0)
                    {
                        mate = ">";
                    }
                    else
                    {
                        mate = "=";
                    }
                }

                //
                // Apply the cigar to seq as best we can.
                //
                int indexIntoSeq = 0;
                string processedCigar = cigar;
                while (processedCigar != "" && processedCigar != "*")
                {
                    string digits = "";
                    while (processedCigar[0] >= '0' && processedCigar[0] <= '9')
                    {
                        digits += processedCigar[0];
                        processedCigar = processedCigar.Substring(1);
                    }

                    int count = Convert.ToInt32(digits);
                    switch (processedCigar[0])
                    {
                        case 'S' :
                            seq = seq.Substring(0, indexIntoSeq) + seq.Substring(indexIntoSeq + count);
                            break;

                        case 'H':
                            break;

                        case 'M':
                        case '=':
                        case 'X':
                        case 'I':   // Unfortunate, but we need to be more sophisticated to handle this.
                            indexIntoSeq += count;
                            break;

                        case 'D' :
                            string Ds = "";
                            for (int i = 0; i < count; i++) {
                                Ds = Ds + "d";
                            }

                            seq = seq.Substring(0,indexIntoSeq) + Ds + seq.Substring(indexIntoSeq);
                            indexIntoSeq += count;
                            break;

                    }

                    processedCigar = processedCigar.Substring(1);   // Strip off the letter
                }



                if (-1 == lowestEmitted)
                {
                    lowestEmitted = pos;
                    Console.Write("" + target + " ");
                    for (int i = lowestEmitted; i < target; i++)
                    {
                        Console.Write(" ");
                    }
                    Console.WriteLine("*");
                }
                string spaces = " ";
                for (int i = lowestEmitted; i < pos; i++) {
                    spaces += " ";
                }
                string direction;
                if ((flags & 0x10) == 0)
                {
                    direction = "+";
                }
                else
                {
                    direction = "-";
                }

                Console.Write("" + pos + spaces);
                for (int i = 0; i < seq.Count(); i++)
                {
                    switch (seq[i])
                    {
                        case 'A':
                            Console.BackgroundColor = ConsoleColor.Red;
                            break;
                        case 'T':
                            Console.BackgroundColor = ConsoleColor.Green;
                            break;
                        case 'C':
                            Console.BackgroundColor = ConsoleColor.Yellow;
                            break;
                        case 'G':
                            Console.BackgroundColor = ConsoleColor.Blue;
                            break;

                        default:
                            Console.BackgroundColor = ConsoleColor.Black;
                            break;
                    }
                    Console.Write(seq[i]);
                }
                Console.BackgroundColor = ConsoleColor.Black;
                Console.WriteLine("\t" + fields[5] + "\t" + direction + "\t" + mate);
            }
        }
    }
}
