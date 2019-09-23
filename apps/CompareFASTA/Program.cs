using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;

namespace CompareFASTA
{
    class Program
    {
        static void Main(string[] args)
        {
            if (args.Count() != 2)
            {
                Console.WriteLine("Usage: ComapreFASTA fasta1 fasta2");
                return;
            }

            ASETools.FASTA one = ASETools.FASTA.loadFromFile(args[0]);
            ASETools.FASTA two = ASETools.FASTA.loadFromFile(args[1]);

            if (one == null || two == null)
            {
                Console.WriteLine("Unable to load at least one FASTA");
                return;
            }

            ASETools.FASTA.compare(one, two);
        }
    }
}
