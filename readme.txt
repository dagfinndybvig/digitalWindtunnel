

EPPLER - DESIGN AND ANALYSIS OF LOW SPEED AIRFOILS         /eppler/readme.txt

The files for this program are in the directory eppler 
  readme.txt      general description
  tm80210.pdf     the original report - with user instructions.
  original.src    the original source - scanned from TM 80210
  profile.f90     the complete source code in modern Fortran
  e1098.dat       a sample case (described in report)
  e1098.out       output for this special case.

The reference documents for this program may be accessed from the web page 
   https://www.pdas.com/epplerrefs.html. The description of the input to 
Profile is found in NASA TM 80210.
 
To compile this program for your computer, use the command
   gfortran  profile.f90 -o profile.exe
Linux and Macintosh users may prefer to omit the .exe on the file name.

This program asks for the name of the input file. This must be 
a file written in the input format to Profile. After reading 
the input data, the program produces a file called profile.out that may be 
printed.


This program, known in the aeronautical community as the Eppler program,
is a combination of three separate algorithms:

   1) a conformal-mapping method for design of airfoils with prescribed
      velocity-distribution characteristics.
   2) a panel method for the analysis of potential flow about given
      airfoils
   3) an integral boundary layer method.

With this combined method, airfoils with prescribed boundary layer
characteristics can be designed and airfoils with prescribed shapes
can be analyzed.

The program is documented in NASA Technical Memorandum 80210 (ref 1)
and it is the basis for the results in Eppler's book (ref 4).

In addition to this program, there is an alternate version that includes
corrections for compressibility. If you would like a copy of that program 
and some additional programs by the author, John Roncz, you may download
from    https://www.pdas.com/packages/roncz.zip

According to the wishes of John Roncz, there is absolutely no charge for 
these programs and you may pass copies on to your colleagues.

The very latest version of the program, ref 3, is available from
    Universitaet Stuttgart
    Institut A fuer Mechanik
    Pfaffenwaldring 9
    D-7000 Stuttgart 80
    Germany.

This version is listed at DM 1200 in the introduction to ref 2.

I have made no changes to the program that will affect results. I have
modified the structure of the source code somewhat to look more like a real 
Fortran 90 program. The original program used a large number of lines of 
source code for the purpose of making plots on a particular plotter in use 
at NASA Langley in the late 70s. I removed all of this with the intent of 
replacing it someday with PostScript, GnuPlot or PCL. Someday hasn't
arrived yet and the inclusion of the new program Pablo with its friendly
user interface makes it unlikely that this will happen.

This program is of considerable historical interest, but does not reflect
the current state-of-the-art in airfoil analysis and design. The programs 
of choice for airfoil studies are Xfoil and INS2D. Xfoil is a potential 
flow and boundary layer approach and INS2D is a Navier-Stokes solver for 
incompressible flow. These are not on this CD-ROM. If you want to get Xfoil 
go to the web site:
  https://raphael.mit.edu/xfoil/

INS2D is distributed through NASA Ames Research Center.

