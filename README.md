[![DOI](https://zenodo.org/badge/21126/oscaribv/fargo_tools.svg)](https://zenodo.org/badge/latestdoi/21126/oscaribv/fargo_tools)

# __fargo_tools__                         
## Written by Barragán O., Rendón F. & Álvarez R. 2015          
#### email: oscaribv@gmail.com			 
#### Updated Sep 04, 2015

* This is a set of tools which allows to work with the outputs of the
 FARGO3D code.
 
* It uses OPENMP to work with different files with several proccesors.

* A fast manual to use fargo_tools is in the fargo_tools can be found in the  
  appendix of Barragán O., 2015, Master Thesis.
  link -> https://www.researchgate.net/publication/281461299_On_the_stability_of_circumplanetary_disks
 
* It is realeased under the GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007.

##  contents				

fargo_tools root directory contains several files and directory. We will
describe this briefly in the following lines.

* README:
    * I am this file.

* GNU_General_Public_License.txt
    * The GNU GENERAL PUBLIC LICENCE downloaded from http://www.gnu.org/licenses/gpl.html

* manual.pdf
    * This is the fast manual to use the fargo_tools package. In order to use correctly the
  fargo_tools routines, you have to read first this brief manual. It includes how to use
  each tool which nices examples.

* input.dat.sample
    * Sample input.dat file used by t2csv. This is an example which can be used for the using t2csv sectionin the manual. You just need chage the name from input.dat.sample to input.dat

* input_ext.dat.sample
    * Sample inputi_ext.dat file used by extracter. This is an example which can be used for the using extracter section in the manual (yet not included!). You just need chage the name from input_ext.dat.sample to input_ext.dat

* ./sources_f90
    * It contains the source FORTRAN90 files. More details of each file and its function
  is given below in the sources_f90 section.

* ./t2csv
    * This directory contains the makefile which compiles the t2csv tool. Once you have 
  compiled this code with the make command, some .o and the t2csv excecutable are 
  created here. 

* ./extracter
    * This directory contains the makefile which compiles the extracter tool. Once you have 
  compiled this code with the make command, some .o and the extracter excecutable are 
  created here. 


##		     	sources_f90				

All f90 files are inside the sources_f90. A brief description
of each file is given in the following lines.

#### create_names.f90
It contains subroutines to work with the filenames.
It is used by t2csv.f90 and extracter.f90

#### input_file.f90
It contains subroutines to read the input files.
It is used by t2csv.f90 and extracter.f90

#### trans_coords.f90
Inside this file there are all the subroutines which transform
the binary FARGO3D output files into csv files.
A detailed description of how the coordinates are transformed will
appear in Barragán O. 2015, Master Thesis.
This file is used by t2csv.f90

#### t2csv.f90
This file transforms the default binary FARGO3D output files to csv files.
It uses create_names.f90, input_file.f90 and trans_coords.f90

#### extracter.f90
This file allows to extract a specific region from the default binary FARGO3D
output files. The output files are also binary files with a lower size.
These files can be transformed to csv format with t2csv.


##		     	PARALLELIZATION				


For high resolution simulations, there are millions of points which need to be 
transformed by the subroutines described previously. The process can take several 
minutes to work with just a file. For runs with hundred of outputs, this can take
some days of data transformation. In order to avoid that, such subroutines  
were parallelized using OPENMP to run in several CPU's at the same time. 

The parallelization has not the purpouse of transforming a file using multiple cores.
Its main aim is to manipulate multiple files at the same time,
each file with a different CPU. The parallelization is worthless if transformation is 
just to one file. But if the number of files to be transformed is greater than one, 
the paralellization will be effective. This parallelization allows a maximum speed up of 
the order of N in some cases, where N is the number of computer cores.
	

##		     	ACKNOWLEDGEMENTS			


* We thank to Dr. Erick Nagel for his support to realize this work and his
time to  read and comment these files.

* Oscar Barragán thanks to Dr. F. Masset and Dra. G. Koenigsberger for the "Escuela de
Verano y supercómputo paralelo y GPU's" in Instituto de Ciencias Físicas, UNAM, Cuernavaca,
Morelos, Méx. 2014. Financed by CONACyT, grant 129343.

## Examples

![](images/planets.png?raw=true)
![](images/csd.png?raw=true)
![](images/fluxes.png?raw=true)
![](images/nice3dview.png?raw=true)
