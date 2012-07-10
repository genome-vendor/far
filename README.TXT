INSTALLATION:
=============

If you downloaded precompiled binaries for linux or windows, just make sure that the corresponding multi-threading library "Intel threading building blocks" (TBB) is available to your system:

A)

LINUX/UNIX:
----------
To make FAR work, please copy the lib file "libtbb.so.2" from ./lib to your lib searchpath or into the folder where FAR is located.

MacOSX:
------
Copy the libtbb.dylib from ./lib into your lib searchpath

WINDOWS:
-------
To make FAR work, copy the "tbb.dll" from YOUR_FAR_DIRECTORY\lib into your c:\windows\system32 folder.

For experts: Building FAR from SVN for your own operating system 
----------------------------------------------------------------
You should be able to compile FAR for any platform. What you need is:

1) check out the SVN repository into a local directory via
svn co https://theflexibleadap.svn.sourceforge.net/svnroot/theflexibleadap YOUR_FAR_DIRECTORY

or use a graphical SVN client like tortoisesvn (http://tortoisesvn.tigris.org/) to check out the repository.

2) The intel threading library (TBB) from: (This you only have to do if you dont have a linux64,Windows or MacOSX running, precompiled lib's for these operating systems can be found in the ./lib folder!)

http://www.threadingbuildingblocks.org/file.php?fid=77
Download the latest stable release (it should work with version >= 3.0)

unpack the archive:

tar xzvf tbb30_20100406oss_lin.tgz

..and run

gmake

in the unpack-folder.

3)
Make the TBB library available in your library searchpath:
- If you have linux64/windows/MaxOSX follow the steps described in the A) section
- If you compiled TBB yourself copy the compiled lib from your ...$TBB_INSTALL_DIR/build/release folder to ./lib of the YOUR_FAR_DIRECTORY directory and to your library searchpath

4) Download cmake from:
http://www.cmake.org/

5) start cmake via gui and open YOUR_FAR_DIRECTORY or use the command line to change to mylocalfolder and type:
cmake .
(dont forget the ".")

CMake allows you to create makefiles or project files for your favorite development evironment via the -G (Generator switch, for visual studio or eclipse). 

For example if you want to build Eclipse project files type
cmake -G "Eclipse CDT4 - Unix Makefiles" .

Please refer to the cmake documentation to see how it works or type "cmake --help".

6) Depending in which way you used CMake you just have to compile the sourcecode via your choosen mechanism. Always make sure the TBB library and the seqan library (incorporated in the repository) are available to the compiler/linker (mylocalfolder/lib).

In the eclipse case this means you should start eclipse, import the project from "mylocalfolder/buidl" and compile FAR after you set libtbb.so path in your eclipse project settings.
 
FAR USAGE:
==========

Please refer to the detailed documentation on:
https://sourceforge.net/apps/mediawiki/theflexibleadap/index.php?title=Main_Page


FAR USAGE (for the impatient):
=============================

input:
------
Far needs at least the following files:
1) FASTA file containing adapter sequences (specified via --adaters)
2) FASTA/Q or CSFASTA/Q file with read (specfied via --source)

mandatory parameters are:
-------------------------

--adapters
--source
--target
--format

for example the command:

far --source reads.fa --target reads_no_adapters.fa --adapters myAdapters.fa --format fasta

will remove all the adapter sequences from myAdapters.fa found in reads.fa and store the results in reads_no_adapters.fa

Additionally FAR will create a reads_no_adapters.fa.lengthdist file, which stores the length distribution of the reads of that file

There are important parameters which you shouldnt miss to configure (like --min-overlap,--min-readlength,--cut-off) - please refer to the wiki page!


Project files/folders:
======================
the ./lib folder contains all DLL/shared libs for the most convenient platforms (windows32, linux64)the ./testset folder contains small test datasets for testing FAR (development of new algorithms, etc.)
the ./include folder contains an adapted version of the seqan framework. This is only necessary for developers. 

Developers testdataset
======================
FAR comes with a test dataset which was designed to test most of the adapter removal scenarios. If you change FAR sourcecode and want to make sure it does the same as before just run
validate.sh 
from the ./testset folder on your compiled far binary.

All test should pass always.

