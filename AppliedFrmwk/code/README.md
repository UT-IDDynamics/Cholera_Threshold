# Source code for producing the results and figures

The code is divided between [Model](Cholera_Function.R), [Parameters](Cholera_params.R), and [the script for calling the model and making the figures](CholeraSims_call.R)

* [Model](Cholera_Function.R): contains the code for the stochastic simulation of the SEIR model for a single parameter run
* [Parameters](Cholera_params.R): contains the range for each of the parameters (min and max) that are not varied (under different scenarios) as well as parameters governing the latin hypercube sampling process
* [the script for calling the model and making the figures](CholeraSims_call.R): contains code for calling the model, including the LHS parameter draws, as well as code for processing the data and generating figures