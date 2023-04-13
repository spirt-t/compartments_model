# Compartments Model

This is the implementation of compartments drug release model, which 
takes into account the processes of drug diffusion and material surface erosion.
The model makes it possible to fit experimental data on the drug release and determine the drug diffusion coefficient and polymer erosion rate.

To run the calculations you should to:  
1. install python version > 3.0 (_for example, by official site https://www.python.org/downloads/_)
2. in the project folder install libraries with command on the command line:
   `pip3 install -r requirements.txt`   
3. in the project folder run script 
   - if you want to take into account only diffusion, use command:
   `python3 ./cmd/without_erosion.py`  
   - if you want to take into account diffusion and surface erosion, use command:   
   `python3 ./cmd/with_erosion.py`  
   - if you want to get the coefficients on the Peppas model from the previous simulation results, use command:   
   `python3 ./cmd/fit.py`  
   

Some help and instructions on how to use the model can be found on the `Manual` tab.  

Test experimental data are in the folder `./testfiles`.

