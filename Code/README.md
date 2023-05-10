The scripts contained in the folder were used to process the data and generate the results reported in the paper.
The function of each one of them is explained in their headers, we give here an overview, following the order in which they should be run:

1. `select_model_order.m` : this scripts is used to select different model classes representing the capacitance-pressure relationship for each capacitive tactile pixel (taxel).
   A calibration step is performed, and the models retrieved are validated in terms of RMSE between the pressure estimates and the ground truth.
   At the end, a visual inspection is needed to select the best trade-off in terms of model accuracy vs model complexity.

2. `calibration_taxels.m` : after selecting the model type to calibrate (in terms of n_p, n_s and n_{ps}), the user can fix it in this code and run it, to retrieve
   the optimal coefficients for each taxel. The results of the calibration are saved such that they can be used as input in the next code.

3. `model_validation.m` : given the optimal coefficients for each taxel in the insole (which are found during the calibration phase), their performance
   is validated comparing the estimates that the model achieve against experimental data.
