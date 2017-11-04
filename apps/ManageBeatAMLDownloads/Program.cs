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
            public readonly int dataWave;
            public readonly string fileType;
            public readonly string eTag;
            public readonly string dataType;
            public readonly string mutation;
        }
        static void Main(string[] args)
        {
        }
    }
}
