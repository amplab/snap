using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;

namespace ManageBeatAMLDownloads
{
    class Program
    {

        class BeatAMLFolder
        {
            public readonly int versionComment;
            public readonly string eTag;
            public readonly string id;
            public readonly int createdByPrincipalId;
            public readonly int versionNumber;
            public readonly int projectId;
            public readonly int modifiedByPrincipalId;
            public readonly long modifiedOn;
            public readonly string concreteType;
            public readonly long createdOn;
            public readonly int benefactorId;
            public readonly string alias;
            public readonly string name;
            public readonly string nodeType;
            public readonly int version;
            public readonly string parentId;

            public static List<BeatAMLFolder> readFile(string filename)
            {
                var inputFile = ASETools.CreateStreamReaderWithRetry(filename);

                if (null == inputFile)
                {
                    Console.WriteLine("BeatAMLFolder.ReadFile(): unable to open " + filename + ", failing.");
                    return null;
                }

                string[] wantedFields = {
                    "folder.versionComment",
                    "folder.eTag",
                    "folder.id",
                    "folder.createdByPrincipalId",
                    "folder.versionNumber",
                    "folder.projectId",
                    "folder.modifiedByPrincipalId",
                    "folder.modifiedOn",
                    "folder.concreteType",
                    "folder.createdOn",
                    "folder.benefactorId",
                    "folder.alias",
                    "folder.name",
                    "folder.nodeType",
                    "folder.versionLabel",
                    "folder.parentId",
                };

                var headerizedFile = new ASETools.HeaderizedFile<BeatAMLFile>(inputFile, false, false, "", wantedFields.ToList());

                List<BeatAMLFile> retVal;

                if (!headerizedFile.ParseFile(parser, out retVal))
                {
                    Console.WriteLine("BeaAMLFile.ReadFile(): Failed to parse input.");
                    return null;
                }

                return retVal;
            }


        }

        class BeatAMLFile
        {
            public readonly long fileSize;
            public readonly string platform;
            public readonly string diagnosis;
            public readonly long modifiedOn;
            public readonly string format;
            public readonly string coreProjectId;
            public readonly string diagnosisClass;
            public readonly int modifiedByPrincipalId;
            public readonly string fileId;
            public readonly long createdOn;
            public readonly long modifiedDate;
            public readonly string speciminType;
            public readonly string name;
            public readonly string coreSampleId;
            public readonly string labId;
            public readonly int versionNumber;
            public readonly string parentId;
            public readonly string dataWave;
            public readonly string fileType;
            public readonly string eTag;
            public readonly string dataType;
            public readonly string mutation;

            public static List<BeatAMLFile> readFile(string filename)
            {
                var inputFile = ASETools.CreateStreamReaderWithRetry(filename);

                if (null == inputFile) {
                    Console.WriteLine("BeatAMLFile.ReadFile(): unable to open " + filename + ", failing.");
                    return null;
                }

                string[] wantedFields = {
                    "file.fileSize",
                    "file.platform",
                    "file.diagnosis",
                    "file.modifiedOn",
                    "file.fileFormat",
                    "file.coreProjectId",
                    "file.diagnosisClass",
                    "file.modifiedByPrincipalId",
                    "file.id",
                    "file.createdOn",
                    "file.modifiedDate",
                    "file.specimenType",
                    "file.name",
                    "file.coreSampleId",
                    "file.labId",
                    "file.versionNumber",
                    "file.parentId",
                    "file.datawave",
                    "file.fileType",
                    "file.eTag",
                    "file.dataType",
                    "file.mutation",
                };

                var headerizedFile = new ASETools.HeaderizedFile<BeatAMLFile>(inputFile, false, false, "", wantedFields.ToList());

                List<BeatAMLFile> retVal;

                if (!headerizedFile.ParseFile(parser, out retVal))
                {
                    Console.WriteLine("BeaAMLFile.ReadFile(): Failed to parse input.");
                    return null;
                }

                return retVal;
            }

            BeatAMLFile(long fileSize_, string platform_, string diagnosis_, long modifiedOn_, string format_, string coreProjectId_, string diagnosisClass_, int modifiedByPrincipalId_, string fileId_,
                long createdOn_, long modifiedDate_, string speciminType_, string name_, string coreSampleId_, string labId_, int versionNumber_, string parentId_, string dataWave_, string fileType_,
                string eTag_, string dataType_, string mutation_)
            {
                fileSize = fileSize_;
                platform = platform_;
                diagnosis = diagnosis_;
                modifiedOn = modifiedOn_;
                format = format_;
                coreProjectId = coreProjectId_;
                diagnosisClass = diagnosisClass_;
                modifiedByPrincipalId = modifiedByPrincipalId_;
                fileId = fileId_;
                createdOn = createdOn_;
                modifiedDate = modifiedDate_;
                speciminType = speciminType_;
                name = name_;
                coreSampleId = coreSampleId_;
                labId = labId_;
                versionNumber = versionNumber_;
                parentId = parentId_;
                dataWave = dataWave_;
                fileType = fileType_;
                eTag = eTag_;
                dataType = dataType_;
                mutation = mutation_;
            }

            static BeatAMLFile parser(ASETools.HeaderizedFile<BeatAMLFile>.FieldGrabber fieldGrabber)
            {
                return new BeatAMLFile(fieldGrabber.AsLongFromSquareBracketLong("file.fileSize"), fieldGrabber.AsStringFromSquareBracketSingleQuoteString("file.platform"), fieldGrabber.AsStringFromSquareBracketSingleQuoteString("file.diagnosis"),
                    fieldGrabber.AsLong("file.modifiedOn"), fieldGrabber.AsStringFromSquareBracketSingleQuoteString("file.fileFormat"), fieldGrabber.AsStringFromSquareBracketSingleQuoteString("file.coreProjectId"),
                    fieldGrabber.AsStringFromSquareBracketSingleQuoteString("file.diagnosisClass"), fieldGrabber.AsInt("file.modifiedByPrincipalId"), fieldGrabber.AsString("file.id"), fieldGrabber.AsLong("file.createdOn"),
                    fieldGrabber.AsLongFromSquareBracketLong("file.modifiedDate"), fieldGrabber.AsStringFromSquareBracketSingleQuoteString("file.specimenType"), fieldGrabber.AsString("file.name"),
                    fieldGrabber.AsStringFromSquareBracketSingleQuoteString("file.coreSampleId"), fieldGrabber.AsStringFromSquareBracketSingleQuoteString("file.labId"), fieldGrabber.AsInt("file.versionNumber"),
                    fieldGrabber.AsString("file.parentId"), fieldGrabber.AsStringFromSquareBracketString("file.datawave"), fieldGrabber.AsStringFromSquareBracketSingleQuoteString("file.fileType"), fieldGrabber.AsString("file.eTag"),
                    fieldGrabber.AsStringFromSquareBracketSingleQuoteString("file.dataType"), fieldGrabber.AsString("file.mutation"));
            }
        }
        static void Main(string[] args)
        {
            var files = BeatAMLFile.readFile(@"f:\temp\all_files.txt");
        }
    }
}
