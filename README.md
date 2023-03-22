# LeMoS-inflowgenerator

At first you have to download an addon folder where necessary libraries for the inflow generator are stored. In the future this will be changed to come in only one item.

User name: lemos, Password: af34ed12

# Workstation installation

1.) check if you have a user-defined OpenFOAM folder: "echo $WM_PROJECT_USER_DIR" should give something like "/home/userName/OpenFOAM/userName-OFxxx"
    --> if it does not exist create it (please refer to the OpenFOAM installation documentation)
    
2.) change into your user-defined OpenFOAM directory

3.) clone into the path of the addon folder: 
    git clone https://github.com/LeMoS-uni-rostock/insightcae.git
    
4.) change directory into the newly created folder "insightcae"
    
5.) change to correct branch
    git checkout next-release

6.) pull the file of the correct branch
    git pull    

7.) change in directory src/addons
    cd src/addons

8.) download add-on: 
    git clone https://github.com/LeMoS-uni-rostock/LeMoS-inflowgenerator.git

9.) change directory to the folder "LeMoS-inflowgenerator" and check the branch
    cd LeMoS-inflowgenerator/
    git branch 
    git checkout main
    
---------------------------------------------------------------------------------

10.) cd extensions/openfoam/inflowGeneratorBC/inflowGeneratorBC/Make.OFxxxx_[...]

10.1) do not change anything in files

10.2) sudo apt update
     sudo apt upgrade

10.3) make sure "armadillo" is installed

10.4) find out which vtk-version you have (e.g. with Synaptic Package Manager)
    make sure "libvtkx-dev" and "libgsl-dev" are installed
    
10.5) in "options" change
    line 8: change according to your vtk version
    line 19: change according to your vtk version

10.6) in "options" change DOR_VERSION according to this table
    OF16ext    1.6.0    010600    extend
    fx31        1.6.1    010601    extend
    fx32        1.6.2    010602    extend
    fx30        1.6.3    010603    extend
    fx41        1.6.4    010604    extend
    of21x        2.1.0    020100    vanilla
    of22x        2.2.0    020200    vanilla
    of22eng    2.2.0    020200    engys
    of23x        2.3.0    020300    vanilla
    of24x        2.4.0    020400    vanilla
    of301        3.0.1    030001    vanilla
    ofplus        4.0.0    040000    esi
    of1806        6.0.0    060000    esi
    of1906        6.5.0    060500    esi
    of2106        6.5.0    060500    esi
    of2112        6.5.5    060505    esi
    ofdev        7.0.0    070000    vanilla


10.7) if you have python3-paraview you should remove it

11.) touch toolkit_export.h

12.) source OpenFOAM

13.) "wmake libso . Make.OFxxxx_[...]"
 
---------------------------------------------------------------------------------
 
14.) if you get this error "fatal error: vtkSmartPointer.h: No such file or directory", you did not specify the correct vtk version in the "options" file
    --> look into /usr/include/vtkx.x --> "wclean" --> repeat the wmake command
     
---------------------------------------------------------------------------------

Be sure to have "libs ("libinflowGeneratorBC.so");" in controlDict.

---------------------------------------------------------------------------------
# Test Case

1.) blockMesh

2.) Run the simulation using "solver    [...];".


---------------------------------------------------------------------------------
---------------------------------------------------------------------------------

# CLUSTER Installation (for OpenFOAM v2106)

1.) download, copy into Thirdparty-v2106 and unpack "VTK-8.2.0"
        make sure that this is the mentioned version in OpenFOAM-v2106(etc/config.sh/VTK

2.) run makeVTK in ThirdParty directory
        check if the installion was correctly (in log file)
        resource OpenFOAM (to expand your $LD_LIBRARY_PATH)
        check if the linking worked correctly e.g. ldd __PATH_TO_OPENFOAM__/ThirdParty-v2106/platforms/linux64Gcc/VTK-8.2.0/lib64/libvtkIOLegacy-8.2.so 
            if "=> not found" appears, something went wrong. 
                check with "echo $LD_LIBRARY_PATH" if the "__PATH_TO_OPENFOAM__ThirdParty-v2106/platforms/linux64Gcc/VTK-8.2.0/lib" or "__PATH_TO_OPENFOAM__ThirdParty-v2106/platforms/linux64Gcc/VTK-8.2.0/lib64" is the correct path to your libraries
                if necessary expand the $LD_LIBRARY_PATH with "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/hp222/Programme/OpenFOAM/ThirdParty-v2106/platforms/linux64Gcc/VTK-8.2.0/lib64"
        
2.) copy boost_1_75_0 in the directory mentioned in Make_XXXXX/options

3.) 




1.) copy files

2.) change Make_xxxx/options
    path to boost
            -I/home/hp222/Programme/boost_1_75_0 \
    path to VTK
            -I/home/hp222/Programme/OpenFOAM/ThirdParty-v1906/platforms/linux64Gcc/VTK-8.2.0/include/vtk-8.2 \
            -L/home/hp222/Programme/OpenFOAM/ThirdParty-v1906/platforms/linux64Gcc/VTK-8.2.0/lib64 -lvtkIOLegacy-8.2 -lvtkCommonSystem-8.2 -larmadillo
