rem This script demonstrates which env variables need to be set under windows
rem to get SNAP to run with HDFS support.
rem
rem To run, ensure that JAVA_HOME and HADOOP_HOME are set in your shell.
rem If not, you can set them here.

rem set JAVA_HOME=c:\program files\java\jre6
rem set HADOOP_HOME=c:\hadoop\hadoop-1.2.0.1.3.0.0-0380


setlocal ENABLEEXTENSIONS
setlocal ENABLEDELAYEDEXPANSION

rem set LIBHDFS_OPTS=-verbose:jni

set LIBHDFS_CLASSPATH=

for %%i in (%HADOOP_HOME%\lib\*.jar) do (
set LIBHDFS_CLASSPATH=!LIBHDFS_CLASSPATH!;%%i
)

for %%i in (%HADOOP_HOME%\*.jar) do (
set LIBHDFS_CLASSPATH=!LIBHDFS_CLASSPATH!;%%i
)

set LIBHDFS_CLASSPATH=!LIBHDFS_CLASSPATH!;%HADOOP_HOME%\conf

..\obj\bin\Debug\x64\snap paired hdfs:///scratch/moonshot/indices/hg19-20-5 example.bam
