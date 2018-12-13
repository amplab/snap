using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;
using System.Net;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Json;
using System.Diagnostics;
using System.Web;
using System.Threading;

namespace GenerateCases
{
    //
    // Creates a set of Cases by looking at the existing MAFs and using the Genome Data Commons API to look at what is available.  It then writes them into the cases file.
    //
    class Program
    {
        class DiseaseSelector
        {
            public DiseaseSelector(string[] diseasesToHandle_)
            {
                diseasesToHandle = diseasesToHandle_.ToList();
                for (int i = 0; i < diseasesToHandle.Count(); i++)
                {
                    diseasesToHandle[i] = diseasesToHandle[i].ToLower();
                }
                handleAllDiseases = diseasesToHandle.Count() == 0;
            }

            public bool useThisProject(string projectId)
            {
                if (handleAllDiseases)
                {
                    return true;
                }

                projectId = projectId.ToLower();
                foreach (var disease in diseasesToHandle)
                {
                    if (projectId.Contains(disease))
                    {
                        return true;
                    }
                }

                return false;
            }

            List<string> diseasesToHandle;
            bool handleAllDiseases;
        }

        static void Main(string[] args)
        {
            var stopwatch = new Stopwatch();
            stopwatch.Start();

            var configuration = ASETools.Configuration.loadFromFile(args);

            if (null == configuration)
            {
                return;
            }

            var diseaseSelector = new DiseaseSelector(configuration.commandLineArgs);

            Dictionary<string, ASETools.DownloadedFile> downloadedFiles = new Dictionary<string, ASETools.DownloadedFile>();
            Dictionary<string, List<ASETools.DerivedFile>> derivedFiles = new Dictionary<string, List<ASETools.DerivedFile>>();

            ASETools.ScanFilesystems(configuration, out downloadedFiles, out derivedFiles);

            var mafManifest = ASETools.MAFInfo.LoadMAFManifest(configuration.mafManifestPathname);

            if (null == mafManifest || mafManifest.Count() == 0)
            {
                Console.WriteLine("Unable to load MAF manifest.  You probably need to generate it first.");
                return;
            }

            var mafManifestFileIds = new List<string>();
            foreach (var mafInfoEntry in mafManifest)
            {
                mafManifestFileIds.Add(mafInfoEntry.Value.file_id);
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

            //string filesForCaseRequest = "{\"op\":\"in\",\"content\":{\"field\":\"cases.case_id\",\"value\":[\"fcf64b55-5c4f-4d82-ac6a-4713d01143cb\",\"25ebc29a-7598-4ae4-ba7f-618d448882cc\",\"fe660d7c-2746-4b50-ab93-b2ed99960553\"]}}";
            //Console.WriteLine(filesForCaseRequest);

            //var response = webClient.DownloadString(ASETools.urlPrefix + "files?pretty=true&filters=" + filesForCaseRequest + "&fields=data_type,updated_datetime,created_datetime,file_name,md5sum,data_format,access,platform,state,file_id,data_category,file_size,type,experimental_strategy,submitter_id,cases.samples.sample_type_id,cases.samples.sample_type,cases.samples.sample_id,analysis.workflow_link");

            //Console.WriteLine("response " + response);

            var filesSerializer = new DataContractJsonSerializer(typeof(ASETools.GDCData<ASETools.GDCFile>));
            //var filesData = (ASETools.GDCData<ASETools.GDCFile>)filesSerializer.ReadObject(new MemoryStream(Encoding.ASCII.GetBytes(response)));

            var allCases = new Dictionary<string, ASETools.GDCCase>();   // Indexed by caseId

            //
            // Call the web server to download all of the cases.  They come back paginated, so we need to step through them.
            //

            var casesSerializer = new DataContractJsonSerializer(typeof(ASETools.GDCData<ASETools.GDCCase>));
            int from = 1;
            int nDuplicateCases = 0;

            string tcgaFiltersString = "{\"op\":\"in\", \"content\": {\"field\": \"project.program.name\", \"value\": " + ASETools.generateFilterList(configuration.programNames) + "}}";

            //Console.WriteLine("Encoded string is " + HttpUtility.UrlEncode(tcgaFiltersString));

            var periodicStopwatch = new Stopwatch();
            periodicStopwatch.Start();

            //
            // Run through all the case Ids we saw and load the info on all cases.
            //
            Console.Write("Loading cases...");

            int nProcessed = 0;
            from = 0;

            while (true)
            {
                var serverData = ASETools.DownloadStringWithRetry(webClient, ASETools.urlPrefix + "cases?from=" + from + "&size=100&filters=" + HttpUtility.UrlEncode(tcgaFiltersString) +
                    "&fields=" + ASETools.GDCCase.fields + "&sort=case_id:asc");

                ASETools.GDCData<ASETools.GDCCase> caseData = (ASETools.GDCData<ASETools.GDCCase>)casesSerializer.ReadObject(new MemoryStream(Encoding.ASCII.GetBytes(serverData)));

                var pagination = ASETools.GDCPagination.extractFromString(serverData);

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

                if (pagination.count != caseData.data.hits.Count())
                {
                    Console.WriteLine("Pagination from server has count mismatch with returned data, " + pagination.count + " != " + caseData.data.hits.Count());
                    return;
                }

                foreach (var case_ in caseData.data.hits)
                {

                    if (allCases.ContainsKey(case_.case_id))
                    {
                        nDuplicateCases++;
                        Console.WriteLine("Duplicate caseID " + case_.case_id);
                        Console.WriteLine(case_.debugString());
                        Console.WriteLine(allCases[case_.case_id].debugString());
                    }
                    else
                    {
                        if (diseaseSelector.useThisProject(case_.project.project_id))
                        {
                            allCases.Add(case_.case_id, case_);
                        }
                    }

                    nProcessed++;
                }

                from += pagination.count;

                if (from >= pagination.total)
                {
                    break;
                }

                if (periodicStopwatch.ElapsedMilliseconds >= 60000)
                {
                    Console.Write(nProcessed + " ");
                    periodicStopwatch.Stop();
                    periodicStopwatch.Reset();
                    periodicStopwatch.Start();
                }

                if (allCases.Count() > 100 && false) break;  // BJB - debugging - just load a few.
            }

            Console.WriteLine(allCases.Count() + " (plus " + nDuplicateCases + " duplicates) in " + ASETools.ElapsedTimeInSeconds(stopwatch) + ", " + allCases.Count() * 1000 / stopwatch.ElapsedMilliseconds + @"/s");

            //
            // Now run through all of the cases looking to see which have files for tumor and matched normal DNA and tumor RNA and copy number.
            //
            int nMissingTumorDNA = 0;
            int nMissingNormalDNA = 0;
            int nMissingMatchedDNA = 0; // We have DNA, but not a pair that corresponds to what was used to generate the MAF
            int nMissingTumorRNA = 0;
            int nMissingMAF = 0;
            int nMissingCopyNumber = 0;
            int nComplete = 0;
            int nWithNormalRNA = 0;
            int nWithTumorMethylation = 0;
			int nWithTNMethylation = 0;
			int nWithNormalCopyNumber = 0;

            Console.Write("Loading files for cases...");
            periodicStopwatch.Stop();
            periodicStopwatch.Reset();
            periodicStopwatch.Start();

            int nCasesProcessed = 0;

            var cases = new Dictionary<string, ASETools.Case>();

			// Load in old cases to get old maf lines
			//var pastCases = ASETools.Case.LoadCases(configuration.casesFilePathname);

			const string mafCondition = "{\"op\":\"in\",\"content\":{\"field\":\"data_format\",\"value\":[\"MAF\"]}}";
            const string bamOrTxtCondition = "{\"op\":\"in\",\"content\":{\"field\":\"data_format\",\"value\":[\"BAM\",\"TXT\"]}}";
 
            foreach (var caseEntry in allCases)
            {
                var gdcCase = caseEntry.Value;

                ASETools.GDCFile tumorRNA = null;
                List<ASETools.GDCFile> tumorDNA = new List<ASETools.GDCFile>();
                List<ASETools.GDCFile> normalDNA = new List<ASETools.GDCFile>();
                ASETools.GDCFile normalRNA = null;
				ASETools.GDCFile tumorMethylation = null;
				ASETools.GDCFile normalMethylation = null;
				List<ASETools.GDCFile> maf = new List<ASETools.GDCFile>();
				ASETools.GDCFile tumorCopyNumber = null;
				ASETools.GDCFile normalCopyNumber = null;
				ASETools.GDCFile tumorFPKM = null;
				ASETools.GDCFile normalFPKM = null;

				from = 0;

                for (; ; )
                {
                    if (periodicStopwatch.ElapsedMilliseconds > 60000)
                    {
                        Console.Write(nCasesProcessed + " ");
                        periodicStopwatch.Stop();
                        periodicStopwatch.Reset();
                        periodicStopwatch.Start();
                    }

                    string fileCondition = "{\"op\":\"in\",\"content\":{\"field\":\"cases.case_id\",\"value\":[\"" + gdcCase.case_id + "\"]}}";
                    string notMAFFilter = "{\"op\":,\"and\",\"content\":[" + fileCondition + "," + bamOrTxtCondition + "]}";
                    string mafFFilter = "{\"op\":,\"and\",\"content\":[" + fileCondition + "," + mafCondition + "]}";

                    var response = ASETools.DownloadStringWithRetry(webClient, 
                        ASETools.urlPrefix + "files?from=" + from + "&size=100&pretty=true&filters=" + fileCondition + 
                        "&fields=data_type,updated_datetime,created_datetime,file_name,md5sum,data_format,access,platform,state,file_id,data_category,file_size,type,experimental_strategy,submitter_id," + 
                        "cases.samples.sample_type_id,cases.samples.sample_type,platform,cases.samples.sample_id,analysis.workflow_link&sort=file_id:asc");

                    ASETools.GDCPagination pagination;
                    ASETools.GDCData<ASETools.GDCFile> fileData;

                    try
                    {
                       fileData  = (ASETools.GDCData<ASETools.GDCFile>)filesSerializer.ReadObject(new MemoryStream(Encoding.ASCII.GetBytes(response)));

                        pagination = ASETools.GDCPagination.extractFromString(response);
                    } catch (SerializationException)
                    {
                        Console.WriteLine("Got SerializationException handling server response, it was probably truncated.  Sleeping 10s and retrying.");
                        Thread.Sleep(10000);
                        continue;   // Since we didn't increase from, we'll just reread the same thing again.
                    }

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

                    if (null == fileData)
                    {
                        Console.WriteLine("Unable to parse file data " + response);
                    } 
                    else
                    {
                        foreach (var file in fileData.data.hits)
                        {
                            if (file.state != "live" && file.state != "submitted" && file.state != "released" || configuration.excludedFileIDs.Contains(file.file_id))
                            {
                                continue;
                            }                         

							if (file.data_format == "BAM" && file.data_type == "Aligned Reads")
							{
								if (0 == file.cases.Count() || 0 == file.cases[0].samples.Count())
								{
									continue;
								}

								int type_id = -1;
								try
								{
									type_id = Convert.ToInt32(file.cases[0].samples[0].sample_type_id);
								}
								catch (FormatException)
								{
									//
									// If it doesn't parse, leave the default -1 value in type_id, which will cause us to ignore this.
									//
								}

								//
								// Type IDs from 1 to 9 are tumor samples and 10 to 19 are matched normal.
								//
								if (type_id < 1 || type_id >= 20)
								{
									//
									// Either a format error, or else a non-tumor, non-normal sample (like a cell line or control sample or something).
									// Ignore it.
									//
									continue;
								}

								bool tumor = type_id < 10;


								if (file.experimental_strategy == "RNA-Seq")
								{
									if (tumor)
									{
										tumorRNA = ASETools.GDCFile.selectNewestUpdated(tumorRNA, file);
									}
									else
									{
										normalRNA = ASETools.GDCFile.selectNewestUpdated(normalRNA, file);
									}
								}
								else if (file.experimental_strategy == "WXS" || file.experimental_strategy == "WGS")
								{
									if (tumor)
									{
										tumorDNA.Add(file);
									}
									else
									{
										normalDNA.Add(file);
									}
								}
							}
							else if (file.data_format == "MAF")
							{
								maf.Add(file);
							}
							else if (file.data_format == "TXT")
							{
								if (0 == file.cases.Count() || 0 == file.cases[0].samples.Count())
								{
									continue;
								}

								int type_id = -1;
								try
								{
									type_id = Convert.ToInt32(file.cases[0].samples[0].sample_type_id);
								}
								catch (FormatException)
								{
									//
									// If it doesn't parse, leave the default -1 value in type_id, which will cause us to ignore this.
									//
								}

								//
								// Type IDs from 1 to 9 are tumor samples and 10 to 19 are matched normal.
								//
								if (type_id < 1 || type_id >= 20)
								{
									//
									// Either a format error, or else a non-tumor, non-normal sample (like a cell line or control sample or something).
									// Ignore it.
									//
									continue;
								}

								bool tumor = type_id < 10;

								if (file.data_type == "Methylation Beta Value")
								{
									if (tumor)
									{
										tumorMethylation = ASETools.GDCFile.selectNewestUpdated(tumorMethylation, file);
									}
									else
									{
										normalMethylation = ASETools.GDCFile.selectNewestUpdated(normalMethylation, file);
									}
								}
								else if (file.data_type == "Copy Number Segment")
								{
									if (tumor)
									{
										tumorCopyNumber = ASETools.GDCFile.selectNewestUpdated(tumorCopyNumber, file);
									}
									else
									{
										normalCopyNumber = ASETools.GDCFile.selectNewestUpdated(normalCopyNumber, file);
									}
								}
								else if (file.data_type == "Gene Expression Quantification" && file.file_name.Contains("FPKM.txt")) // we only want normal FPKM
								{
									if (tumor)
									{
										tumorFPKM = ASETools.GDCFile.selectNewestUpdated(tumorFPKM, file);
									}
									else
									{
										normalFPKM = ASETools.GDCFile.selectNewestUpdated(normalFPKM, file);
									}
								}

								
							}
                        } // Foreach file in this batch
                    } // If file data parsed

                    from += pagination.count;


                    if (from >= pagination.total)
                    {
                        break;
                    }

                } // forever (looping over file pagination)

				var mafFileIds = maf.Select(x => x.file_id).ToList();

                if (tumorDNA.Count() == 0)
                {
                    nMissingTumorDNA++;
                }
                else if (normalDNA.Count() == 0)
                {
                    nMissingNormalDNA++;
                }
                else if (null == tumorRNA)
                {
                    nMissingTumorRNA++;
                }
                else if (mafFileIds.Intersect(mafManifestFileIds).Count() == 0)
                {
                    nMissingMAF++;
                }
                else if (tumorCopyNumber == null)
                {
                    nMissingCopyNumber++;
                }
                else
                {
                    ASETools.GDCFile normalDNAFile = null;
                    ASETools.GDCFile tumorDNAFile = null;

                    //
                    // Just take the newest ones for now.
                    //
                    foreach (var file in normalDNA)
                    {
                        normalDNAFile = ASETools.GDCFile.selectNewestUpdated(normalDNAFile, file);
                    }

                    foreach (var file in tumorDNA)
                    {
                        tumorDNAFile = ASETools.GDCFile.selectNewestUpdated(tumorDNAFile, file);
                    }

                    if (normalDNAFile == null) {
                        nMissingMatchedDNA++;
                    } else {
                        nComplete++;

                        var case_ = new ASETools.Case();
                        case_.case_id = gdcCase.case_id;
                        case_.tumor_dna_file_id = tumorDNAFile.file_id;
                        case_.normal_dna_file_id = normalDNAFile.file_id;
                        case_.tumor_rna_file_id = tumorRNA.file_id;

                        case_.tumor_dna_size = tumorDNAFile.file_size;
                        case_.normal_dna_size = normalDNAFile.file_size;
                        case_.tumor_rna_size = tumorRNA.file_size;

                        case_.tumor_rna_file_bam_md5 = tumorRNA.md5sum;
                        case_.normal_dna_file_bam_md5 = normalDNAFile.md5sum;
                        case_.tumor_dna_file_bam_md5 = tumorDNAFile.md5sum;

                        case_.project_id = gdcCase.project.project_id;

                        if (downloadedFiles.ContainsKey(case_.tumor_dna_file_id)) {
                            case_.tumor_dna_filename = downloadedFiles[case_.tumor_dna_file_id].fileInfo.FullName;
                        }

                        if (downloadedFiles.ContainsKey(case_.normal_dna_file_id)) {
                            case_.normal_dna_filename = downloadedFiles[case_.normal_dna_file_id].fileInfo.FullName;
                        }

                        if (downloadedFiles.ContainsKey(case_.tumor_rna_file_id)) {
                            case_.tumor_rna_filename = downloadedFiles[case_.tumor_rna_file_id].fileInfo.FullName;
                        }

                        if (downloadedFiles.ContainsKey(case_.tumor_dna_file_id))
                        {
                            case_.tumor_dna_filename = downloadedFiles[case_.tumor_dna_file_id].fileInfo.FullName;
                        }

                        if (null != normalRNA)
                        {
                            nWithNormalRNA++;

                            case_.normal_rna_file_id = normalRNA.file_id;
                            case_.normal_rna_size = normalRNA.file_size;
                            case_.normal_rna_file_bam_md5 = normalRNA.md5sum;

                            if (downloadedFiles.ContainsKey(case_.normal_rna_file_id)) {
                                case_.normal_rna_filename = downloadedFiles[case_.normal_rna_file_id].fileInfo.FullName;
                            }
                        }

						var hasTumorMethylation = false;

                        if (null != tumorMethylation)
                        {
							hasTumorMethylation = true;
                            nWithTumorMethylation++;

							case_.tumor_methylation_file_id = tumorMethylation.file_id;
                            case_.tumor_methylation_size = tumorMethylation.file_size;
                            case_.tumor_methylation_file_md5 = tumorMethylation.md5sum;

                            if (downloadedFiles.ContainsKey(case_.tumor_methylation_file_id)) {
                                case_.tumor_methylation_filename = downloadedFiles[case_.tumor_methylation_file_id].fileInfo.FullName;
                            }
                        }

						if (null != normalMethylation)
						{
							// keep track of samples that have both tumor and normal methylation
							if (hasTumorMethylation)
							{
								nWithTNMethylation++;
							}

							case_.normal_methylation_file_id = normalMethylation.file_id;
							case_.normal_methylation_size = normalMethylation.file_size;
							case_.normal_methylation_file_md5 = normalMethylation.md5sum;

							if (downloadedFiles.ContainsKey(case_.normal_methylation_file_id))
							{
								case_.normal_methylation_filename = downloadedFiles[case_.normal_methylation_file_id].fileInfo.FullName;
							}
						}

					    case_.tumor_copy_number_file_id = tumorCopyNumber.file_id;
                        case_.tumor_copy_number_size = tumorCopyNumber.file_size;
                        case_.tumor_copy_number_file_md5 = tumorCopyNumber.md5sum;

					    if (downloadedFiles.ContainsKey(case_.tumor_copy_number_file_id))
                        {
                            case_.tumor_copy_number_filename = downloadedFiles[case_.tumor_copy_number_file_id].fileInfo.FullName;
                        }

						if (null != normalCopyNumber)
						{
							nWithNormalCopyNumber++;

							case_.normal_copy_number_file_id = normalCopyNumber.file_id;
							case_.normal_copy_number_size = normalCopyNumber.file_size;
							case_.normal_copy_number_file_md5 = normalCopyNumber.md5sum;

							if (downloadedFiles.ContainsKey(case_.normal_copy_number_file_id))
							{
								case_.normal_copy_number_filename = downloadedFiles[case_.normal_copy_number_file_id].fileInfo.FullName;
							}
						}

						if (null != tumorFPKM)
						{
							case_.tumor_fpkm_file_id = tumorFPKM.file_id;
							case_.tumor_fpkm_size = tumorFPKM.file_size;
							case_.tumor_fpkm_file_md5 = tumorFPKM.md5sum;

							if (downloadedFiles.ContainsKey(case_.tumor_fpkm_file_id))
							{
								case_.tumor_fpkm_filename = downloadedFiles[case_.tumor_fpkm_file_id].fileInfo.FullName;
							}
						}

						if (null != normalFPKM)
						{
							case_.normal_fpkm_file_id = normalFPKM.file_id;
							case_.normal_fpkm_size = normalFPKM.file_size;
							case_.normal_fpkm_file_md5 = normalFPKM.md5sum;

							if (downloadedFiles.ContainsKey(case_.normal_fpkm_file_id))
							{
								case_.tumor_fpkm_filename = downloadedFiles[case_.normal_fpkm_file_id].fileInfo.FullName;
							}
						}

						case_.maf_file_id = mafFileIds.Intersect(mafManifestFileIds).ToList()[0];

                        if (downloadedFiles.ContainsKey(case_.maf_file_id))
                        {
                            case_.maf_filename = downloadedFiles[case_.maf_file_id].fileInfo.FullName;
                        }

                        case_.sample_ids = gdcCase.sample_ids.ToList();

                        cases.Add(case_.case_id, case_);
                    }
                }

                nCasesProcessed++;

            } // foreach case


            Console.WriteLine(ASETools.ElapsedTimeInSeconds(periodicStopwatch));

            Console.WriteLine("Of " + allCases.Count() + ", " + nMissingTumorDNA + " don't have tumor DNA, " + nMissingNormalDNA + " don't have matched normal DNA, " + nMissingTumorRNA + " don't have tumor RNA, " 
                + nMissingMAF + " don't have MAFs," + nMissingCopyNumber + " don't have copy number, and " + nComplete + " are complete.  Of the complete, " + nWithNormalRNA + " also have normal RNA, " 
				+ nWithTumorMethylation + " have tumor methylation, and " + nWithTNMethylation + " have both normal and tumor methylation, and " 
				+ nWithNormalCopyNumber + " have normal copy number.");

			ASETools.Case.SaveToFile(cases, configuration.casesFilePathname);

            Console.WriteLine("Overall run took " + ASETools.ElapsedTimeInSeconds(stopwatch));
        } // Main

    }
}
