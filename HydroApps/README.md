# README #

This folder would normally document the hydrological models, both lumped and semi-distributed.

### What includes in SNOWModel ###

* The Fortran 95/2003 code with Visual Studio
* The SNOWModel is self-contained.
* SNOW model is a 2 parameter model, used to simulate the snow melting process

* fitModel.f90 is the main code to setup and call SNOWModel
* SnowMod.f90 is the subroutine for Snow Model
* fitModel.txt is the setup for parameters, data, and warmup period
* InputData.csv is the data, includs rainfall, temperature, simulation, and date 

### What includes in AWBM ###

* The Fortran 95/2003 code with Visual Studio
* The AWBM is self-contained.
* AWBM model is a 8 parameter model, used to simulate the hydrological process

* fitModel.f90 is the main code to setup and call AWBM
* SnowMod.f90 is the subroutine for Snow Model
* fitModel.txt is the setup for parameters, data, and warmup period
* InputData.csv is the data, includs rainfall, temperature, simulation, and date
* AWBM use three paralle tank to simulate the spatial variability of soil moisutre
