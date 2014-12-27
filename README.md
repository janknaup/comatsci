# comatsci #
**co**mputational **mat**erials **sci**ence utility package
Version 1.5.0

Jan M. Knaup, Bremen Center for Compuational Materials Science
Jan.Knaup@gmail.com

## Requirements ##

comatsci requires the following software/libraries to be present 
on a system:

* Python version >= 2.7 (not fully compatible with python 3)
* POSIX compatible build environment and header files to compile python c-extensions
* numpy
* h5py

Note that windows is not an officially supported platform for 
comatsci.

## Installation ##

comatsci should be installed via the supplied setup script. 

Execute python setup.py install as root, to install comatsci 
system-wide. Any other installation tree can be chosen with the 
--prefix option to install, e.g. --prefix=~ to install into 
$HOME/bin and $HOME/lib. The setup script also offers options to 
generate installers for different operating systems, depending on 
the platform and installed python version, please refer to the 
integrated documentation available by calling python setup.py 
--help.

## License ##

comatsci is provided without any warranty as open software under 
the terms of the Non-commercial Open Software License v3.0. 
Please refer to the attached LICENSE file or appendix