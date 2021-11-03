using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace FindHole
{
    class Program
    {
        static void Main(string[] args)
        {
            if (args.Count() != 1 && args.Count() != 2)
            {
                Console.WriteLine("usage: FindHole inputFilename <starting offset>");
                Console.WriteLine("Finds any regions of the file with one or more aligned 4K regions of all zeroes");
                return;
            } // usage

            var inputfile = new FileStream(args[0], FileMode.Open);
            if (inputfile == null)
            {
                Console.WriteLine("Unable to open input file " + args[0]);
                return;
            }

            long lastDot = 0;
            long address = 0;
            const int readSize = 4096;
            
            if (args.Count() == 2)
            {
                address = Convert.ToInt64(args[1]);
                inputfile.Seek(address, SeekOrigin.Begin);
                lastDot = address;
            }


            long beginOfAllZeroesRegion = -1;

            var buffer = new byte[readSize];
            int bytesInBuffer;

            while (0 != (bytesInBuffer = inputfile.Read(buffer, 0, readSize)))
            {
                bool anyNonZero = false;
                for (int i = 0; i < bytesInBuffer; i++)
                {
                    if (buffer[i] != '\0')
                    {
                        anyNonZero = true;
                        break;
                    }
                }

                if (anyNonZero)
                {
                    if (beginOfAllZeroesRegion != -1)
                    {
                        Console.WriteLine(beginOfAllZeroesRegion + " - " + address);
                        beginOfAllZeroesRegion = -1;
                    }
                } else
                {
                    if (beginOfAllZeroesRegion == -1)
                    {
                        beginOfAllZeroesRegion = address;
                    }
                } // anyNonZero

                address += bytesInBuffer;

                if (address - lastDot > 1024 * 1024 * 1024)
                {
                    lastDot = address;
                    Console.Write(".");
                }
            }
            
            if (beginOfAllZeroesRegion != -1)
            {
                Console.WriteLine(beginOfAllZeroesRegion + " - " + address + " (EOF)");
            }

            inputfile.Close();
        } // Main
    } // Program
} // FindHole
