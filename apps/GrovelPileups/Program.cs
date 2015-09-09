using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace GrovelPileups
{
    class Program
    {
        static void Main(string[] args)
        {
            foreach (var file in Directory.GetFiles(@"f:\sequence\reads\tcga\cdkn2a_pileup\")) {
                if (file.Count() < 7 || file.Substring(file.Count() - 7) != ".pileup") continue;
                string[] lines = File.ReadAllLines(file);

                for (int i = 1; i < lines.Count(); i++)
                {
                    string[] fields = lines[i].Split('\t');

                    if (fields.Count() >= 5)
                    {
                        int t = 0, c = 0, g = 0, a = 0;

                        for (int j = 0; j < fields[4].Count(); j++ )
                        {
                            char putativeBase = fields[4][j];
                            if (putativeBase == 't' || putativeBase == 'T') t++;
                            else if (putativeBase == 'c' || putativeBase == 'C') c++;
                            else if (putativeBase == 'g' || putativeBase == 'G') g++;
                            else if (putativeBase == 'a' || putativeBase == 'A') a++;
                        }

                        if (t + c + g + a > 15) {
                            int nPossible = 0;
                            if (t * 10 >= (t + c + g + a) * 3) nPossible++;
                            if (c * 10 >= (t + c + g + a) * 3) nPossible++;
                            if (g * 10 >= (t + c + g + a) * 3) nPossible++;
                            if (a * 10 >= (t + c + g + a) * 3) nPossible++;

                            if (nPossible > 1) {
                                Console.WriteLine("Possibly heterozygous locus in file " + file + " line " + lines[i]);
                            }
                        }
                    }
                }
            }
        }
    }
}
