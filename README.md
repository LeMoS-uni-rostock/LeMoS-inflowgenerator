# LeMoS-inflowgenerator


User name: lemos
Password: af34ed12


1.) checkout InsightCAE: 
    git clone https://github.com/LeMoS-HRO/insightcae.git
    
2.) change to correct branch
    git checkout next-release

3.) pull the file of the correct branch
    git pull    

4.) change in directory src/addons
    cd src/addons

5.) downlaod add-on: 
    git clone https://github.com/LeMoS-HRO/insight-inflowgenerator.git

6.) make sure that you are in branch master
    cd insight-inflowgenerator/
    git branch 
    git checkout master

7.) change directory to extensions/openfoam/inflowGeneratorBC
    cd extensions/openfoam/inflowGeneratorBC

do not change nayting in files 
change DOR_VERSION according to this table 
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
of2112        6.5.5    060505    esi
ofdev        7.0.0    070000    vanilla

Change VTK version according to your installation (recommended version )

if you dont have  gsl and vtk library

sudo apt update
sudo apt upgrade

(synaptic will show GUI of installation packages)
(if you have python3-paraview you should remove it)

apt install libgsl-dev libvtk9-dev -y

## libgsl-dev ... vtk-xxx/libvtkXXX-dev müssen installiert sein
touch toolkit_export.h (keine Ahnung warum da eine leere Datei erstellt werden muss)

(
echo $WM_PROJECT_USER_DIR
if it does not exsit creat it
normally/home/mkh/OpenFOAM/"your computer name"-v2012
)

cd OpenFOAM/"your openfoam user directory folder"/insightcae/src/addons/insight-inflowgenerator/extensions/openfoam/inflowGeneratorBC
8.) wmake libso . Make.OF......
 please remember to change the vkt include in option in Make... directory
 

libs ("libinflowGeneratorBC.so"); //in controldict    

Error :undefined symbol: wrapper2_dsyevd_
solution: 











CLUSTER Installation (for OpenFOAM v2106)

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
