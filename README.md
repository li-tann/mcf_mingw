# mcf_mingw

This repository compile GAMMA's mcf (mininum cost flow) code with compiler mingw7.3.0 (provided by QtCreator 4.8.0(Qt 5.12.0))

## mcf

project mcf does not depend on any third-party library

mcf.exe will generated under the path 'mcf.pro/../../bin/plugin' after compile

## use_mcf

project use_mcf was created to call mcf.exe more conveniently

such as load tif image & convert to binary file by GDAL library (which has not been implemented yet)
