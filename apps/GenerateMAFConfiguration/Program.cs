using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Net;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Json;
using System.Diagnostics;
using System.Web;
using ASELib;
using System.IO;

namespace GenerateMAFConfiguration
{
    class Program
    {
        static void Main(string[] args)
        {
            var stopwatch = new Stopwatch();
            stopwatch.Start();

            var configuration = ASETools.Configuration.loadFromFile(args);

            if (configuration == null)
            {
                return;
            }

            var webClient = ASETools.getWebClient();

            if (null == webClient)
            {
                //
                // Probably a version change at gdc.  Give up.
                //
                Console.WriteLine("Unable to get web client, quitting.");
                return;
            }

            //
            // Filter for MAFs in the right projects.  We narrow more later.
            //
            string mafFilesRequest = "{\"op\": \"and\", \"content\": [" +
                                        "{\"op\":\"in\",\"content\":{\"field\":\"data_format\",\"value\":[\"MAF\"]}}," +
                                        "{\"op\":\"in\",\"content\":{\"field\":\"cases.project.program.name\",\"value\":" + ASETools.generateFilterList(configuration.programNames) + "}}" +
                                     "]}";

            var filesSerializer = new DataContractJsonSerializer(typeof(ASETools.GDCData<ASETools.GDCFile>));

            List<ASETools.GDCFile> selectedMAFs = new List<ASETools.GDCFile>();

            Console.Write("Downloading file metadata...");

            Stopwatch printoutTimer = new Stopwatch();
            printoutTimer.Start();

            int from = 1;
            for (; ; )
            {
                var response = webClient.DownloadString(ASETools.urlPrefix + "files?from=" + from + "&size=100&pretty=true&filters=" + mafFilesRequest + "&fields=data_type,updated_datetime,created_datetime,file_name,md5sum,data_format,access,platform,state,file_id,data_category,file_size,type,experimental_strategy,submitter_id,cases.samples.sample_id,analysis.workflow_link&sort=file_id:asc");
                ASETools.GDCData<ASETools.GDCFile> fileData = (ASETools.GDCData<ASETools.GDCFile>)filesSerializer.ReadObject(new MemoryStream(Encoding.ASCII.GetBytes(response)));

                if (null == fileData || null == fileData.data || null == fileData.data.hits)
                {
                    Console.WriteLine("Unable to parse server response for file data.");
                    return;
                }

                //
                // For whatever reason the JSON decoder is failing to get the pagination information, so we extract it directly.
                //
                var pagination = ASETools.GDCPagination.extractFromString(response);

                if (null == pagination)
                {
                    Console.WriteLine("Couldn't parse pagination from server on cases download, from = " + from);
                    return;
                }

                if (pagination.from != from)
                {
                    Console.WriteLine("Pagination from server shows data starting at the wrong place, " + pagination.from + " != " + from);
                    return;
                }

                if (pagination.count != fileData.data.hits.Count())
                {
                    Console.WriteLine("Pagination from server has count mismatch with returned data, " + pagination.count + " != " + fileData.data.hits.Count());
                    return;
                }

                foreach (var file in fileData.data.hits) 
                {
                    if (file.data_type.ToLower() == "aggregated somatic mutation" && file.file_name.ToLower().IndexOf(configuration.mutationCaller.ToLower()) != -1)
                    {
                        selectedMAFs.Add(file);
                    }
                }


                from += pagination.count;

                if (from >= pagination.total)
                {
                    break;
                }

                if (printoutTimer.ElapsedMilliseconds >= 60000)
                {
                    Console.Write(" " + from + @"/" + pagination.total);
                    printoutTimer.Stop();
                    printoutTimer.Reset();
                    printoutTimer.Start();
                }
            } // for ever (pagination loop with internal exit)

            Console.WriteLine();


            var manifestWriter = ASETools.CreateStreamWriterWithRetry(configuration.mafManifestPathname);

            if (null == manifestWriter)
            {
                Console.WriteLine("Unable to open manifest file for write.");
                return;
            }

            manifestWriter.WriteLine("MAF Manifest v1.0 generated at " + DateTime.Now.ToLocalTime().ToString());

            foreach (var file in selectedMAFs)
            {
                manifestWriter.WriteLine(file.file_id + "\t" + file.md5sum + "\t" + file.file_name);
            }

            manifestWriter.WriteLine("**done**");
            manifestWriter.Close();

            Console.WriteLine("Selected " + selectedMAFs.Count() + " MAFs and saved them to " + configuration.mafManifestPathname + " in " + ASETools.ElapsedTimeInSeconds(stopwatch) + ".");
        }
    }
}
