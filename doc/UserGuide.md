# LCS Plugin for ParaView
This package has been created to extract coherent structures of the flow. The resulting output is a manifold with consistent triangulation orientation.

![Max Strain LCS of ABC flow](https://raw.github.com/ameli/lcs/master/doc/figures/1.png "Max Strain LCS of ABC flow")
<a name="figAAAexample"></a>

# Build Source
Source code can be obtained from [github repository](https://github.com/ameli/lcs). Create a build directory and download the source code inside the build directory:

    $ git clone git://github.com/ameli/lcs
    $ cd lcs
    $ mkdir build

Now, build the code:

    $ cd build
    $ cmake ..
    $ make
    $ make install

Installation directory is in _install_ directory. For system wide installation, users may change _PROJECT_INSTALL_PATH_ in cmake. For UNIX systems a preferred system wide path is `/usr/local/bin`.

# Program Usage

The programcan be used either in terminal or as a ParaView plugin. Here are examples to use it in terminal:

    $ lcs /input-path/inputfile.vtk outputfile.vtk

For help, use `-h`, for get program information use `-i` and for see the license, use `-l` options.

    $ lcs -h
    $ lcs -i
    $ lcs -l

The input file should be a `structured point` data. Also Cauchy-Green tensors should be included in `PointData` of the file as a tensor field. A sample of input data with few number line are as following:

    # vtk DataFile Version 3.0
    t = 0.000000000
    ASCII
    
    DATASET STRUCTURED_POINTS
    DIMENSIONS 101 101 101
    ORIGIN -3.141593 -3.141593 -3.141593
    SPACING 0.062831853 0.062831853 0.062831853
    
    POINT_DATA 1030301
    TENSORS CG double
    1841.529748 -571.912532 -688.962477
    -571.912532 178.498989 212.501996
    -688.962477 212.501996 260.187902
    ...

# ParaView Plugin Usage

For using the program as a plugin, should select the `BUILD_PARAVIEW_PLUGIN` while configuring in `ccmake`. After compiling the code, the two shared libraries `libDeformationPlugin.so` and `linIdentifyStructuresPlugin.so` will be created in the `bin` folder. Depending on the operating system, the extension of shared libraries are different. In Linux, the expedition is `.so`, while in Mac and Windows are `.dll` and `.dylib` respectively.

Figures [2](#PV1),[3](#PV2) and [4](#PV3) illustrates how to load the plugin and apply it to make a pipeline in ParaView. loading the plugin in ParaView. Users may find the filter in __Filter__ menu under __Extensions__.

![Manage Plugins](https://raw.github.com/ameli/lcs/master/doc/figures/PV-1.png "Manage Plugins")
<a name="PV1"></a>

![Load Plugins](https://raw.github.com/ameli/lcs/master/doc/figures/PV-2.png "Load Plugins")
<a name="PV2"></a>

![Apply Plugins](https://raw.github.com/ameli/lcs/master/doc/figures/PV-3.png "Apply Plugins")
<a name="PV3"></a>

# Options
Users can choose which type of deformation be used for extracting LCS. The options are

* Max Strain LCS
* Min Strain LCS
* Shear LCS (TODO)

Options of plugins as shown in figure [5](#PV4).

![Setting of the Plugins](https://raw.github.com/ameli/lcs/master/doc/figures/PV-4.png "Setting of the Plugins")
<a name="PV4"></a>

# Example

Figure [6](#2) shows max strain LCS for ABC flow. The resolution of initial flowmap grid is 101x101x101.

![Max Strain LCS for ABC flow](https://raw.github.com/ameli/lcs/master/doc/figures/2.png "Max Strain LCS for ABC flow")
<a name="2"></a>

# License

This source code is provided _as-is_, without any express or implied warranty. In no event will the author be held liable for any damages arising from the use of this software.

Permission is granted to anyone to use this software for non commercial purpose, to alter it and redistribute it freely, subject to the following restrictions:

* The origin of this source code must not be misrepresented; you must not claim that you wrote the original source code. 

* You may modify or alter the source code at your own risk, provided that such modifications are extensively commented. Altered source versions must not be misrepresented as being the original source code.

* Source code may not be redistributed in any way. You may distribute the target binary for non-commercial purposes.

* If you use this source code in a non-commercial product, an acknowledgment in the product documentation would be appreciated.

* This notice may not be removed or altered from any source distribution.

# Bug Report

Any bug reports and comments are appreciated. You may report bugs at [github](https://github.com/ameli/lcs), or send email.

# Acknowledgement

This work was supported by the National Science Foundation, award number 1047963.