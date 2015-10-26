using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace ComputeP53Levels
{
    class Program
    {
        static public string GenerateP53Line(string analysisID, string disease_abbr, string allcountFile)
        {
            return "";
        }
        static void Main(string[] args)
        {
            var outputFile = new StreamWriter(@"f:\temp\expression\p53-levels.txt");
            var outputFileNormal = new StreamWriter(@"f:\temp\expression\p53-levels-normal.txt");

            var experimentReader = new StreamReader(@"f:\temp\expression\experiments.txt");
            experimentReader.ReadLine();    // Skip the header

            string experimentLine;
            while (null != (experimentLine = experimentReader.ReadLine()))
            {
                string[] fields = experimentLine.Split('\t');


            }

            outputFile.Close();
            outputFileNormal.Close();
        }
    }
}
