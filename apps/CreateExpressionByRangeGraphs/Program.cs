using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ExpressionLib;

namespace CreateExpressionByRangeGraphs
{
    class Program
    {

        static void ProcessFile(string filename)
        {
            if (filename.IndexOf(".txt") != filename.Count() - 4)
            {
                Console.WriteLine("ProcessFile: filename doesn't end in .txt: " + filename);
                return;
            }

            string outputFilename = filename.Substring(0, filename.Count() - 4) + "_by_range_graphs.txt";

            var reader = ExpressionTools.CreateStreamReaderWithRetry(filename);
            var writer = ExpressionTools.CreateStreamWriterWithRetry(outputFilename);

            //
            // We pull out the columns of the input of the form "...n mutation  mean", and then pivot them into rows and genes into columns.  Read the header to find them.
            //

            string header = reader.ReadLine();
            var headerFields = header.Split('\t');

            var inputFields = new List<string[]>();

            string line;
            while ((line = reader.ReadLine()) != null) {
                inputFields.Add(line.Split('\t'));
            }

            reader.Close();

            //
            // Now find the columns to turn into rows.  We sort by 0, 1, >1 mutation to make the final graph building in Excel easier
            //
            List<int>[] columns = new List<int>[6];
            for (int i = 0; i < columns.Count(); i++)
            {
                columns[i] = new List<int>();
            }

            for (int i = 0; i < headerFields.Count(); i++)
            {
                if (headerFields[i].IndexOf("0 mutation mean") != -1 || headerFields[i].IndexOf("0 mutation  mean") != -1)
                {
                    columns[3].Add(i);
                }
                else if (headerFields[i].IndexOf(">1 mutation mean") != -1 || headerFields[i].IndexOf(">1 mutation  mean") != -1)    // Have to do this before 1, because 1's a substring
                {
                    columns[5].Add(i);
                }
                else if (headerFields[i].IndexOf("1 mutation mean") != -1 || headerFields[i].IndexOf("1 mutation  mean") != -1)
                {
                    columns[4].Add(i);
                } 
                else if (headerFields[i].IndexOf("0 mutation exclusive mean") != -1 || headerFields[i].IndexOf("0 mutation  exclusive mean") != -1) 
                {
                    columns[0].Add(i);
                }
                else if (headerFields[i].IndexOf(">1 mutation exclusive mean") != -1 || headerFields[i].IndexOf(">1 mutation  exclusive mean") != -1) // Have to do this before 1, because 1's a substring
                {
                    columns[2].Add(i);
                }
                else if (headerFields[i].IndexOf("1 mutation exclusive mean") != -1 || headerFields[i].IndexOf("1 mutation  exclusive mean") != -1) {
                    columns[1].Add(i);
                }
            }

            var allColumns = new List<int>();
            for (int i = 0; i < columns.Count(); i++)
            {
                allColumns.AddRange(columns[i]);
            }

            //
            // And print out the result.  Print the gene names along the top row.
            //
            for (int i = 0; i < inputFields.Count(); i++)
            {
                writer.Write("\t" + ExpressionTools.ConvertToExcelString(inputFields[i][0]));
            }
            writer.WriteLine();

            //
            // The first column is the 

            foreach (var column in allColumns)
            {
                writer.Write(headerFields[column]);
                for (int i = 0; i < inputFields.Count(); i++)
                {
                    if (inputFields[i][column] == "*") 
                    {
                        writer.Write("\t"); // Replace * with a blank, since that renders as a gap in an excel graph, while * renders as 0.
                    }
                    else
                    {
                        writer.Write("\t" + inputFields[i][column]);
                    }
                }
                writer.WriteLine();
            }

            writer.Close();


        }
        static void Main(string[] args)
        {
            string baseDirectory = @"f:\temp\expression\";

            // No need to do the Bonferroni ones, since they have the same values in them; only p-values are corrected.
            ProcessFile(baseDirectory + "AlleleSpecificExpressionDistributionByMutationCount.txt"); 
            ProcessFile(baseDirectory + "ExpressionDistributionByMutationCount.txt"); 
        }
    }
}
