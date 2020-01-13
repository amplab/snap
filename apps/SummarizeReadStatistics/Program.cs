using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;

namespace SummarizeReadStatistics
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

            if (commonData.listOfCases.Any(_ => _.read_statictics_filename == ""))
            {
                Console.WriteLine("At least one case doesn't have read statistics.");
                return;
            }

            foreach (var case_ in commonData.listOfCases)
            {
                var readStatictics = ASETools.ReadStatistics.ReadFromFile(case_.read_statictics_filename);


            }
        }
    }
}
