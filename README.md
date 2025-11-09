# tdse
numerical algorithm for Time-Dependent Schroedinger Equation

This repository contains the following programs:

## fortran
This directory contains programs written in fortran.

(1) Numerical Integration of the time-dependent Schroedinger
equation by the explicit symetric multistep scheme

reference： T.Iitaka, Phys. Rev. E49 (1994) 4684.
https://doi.org/10.1103/PhysRevE.49.4684

(1-1) Stability Check

ST2      FOR       
ST4      FOR       
ST6      FOR       

(1-2) Time-evolution of Eigenstates

BTST2    FOR    
BTST4    FOR    
BTST6    FOR    
BTSTCH   FOR    
BTSTCN   FOR    

(1-3) Time-evolution of Scatteringstates

SC2      FOR    
SC4      FOR    
SC6      FOR    
SCCH     FOR    
SCCN     FOR    

(2) Numerical Integration of the time-dependent Schroedinger
equation by Leap Frog method.

BOX2     FOR    
BOXN     FOR    
1D       FOR    



## python_gemini
This directory contains python programs translated from the fortran programs with the help of Gemini 2.5 flash.
https://g.co/gemini/share/3b6a7390eb48

## python
This directory contains python programs further edited by human.

## doc
This directory contains pdf file of 
"Introduction to Quatum Dynamics" by Toshiaki Iitaka.
(Maruzen Publish. Co., 1994,Tokyo; Parity Physics Course, Close Up)
https://www.amazon.co.jp/dp/4621039717/

