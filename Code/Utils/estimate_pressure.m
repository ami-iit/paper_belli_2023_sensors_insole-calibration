function [estimated_pressure_all, broken_index] = estimate_pressure(validation, model_info, num_taxels, Phi_p, broken_index_valid)
% This function takes as inputs:
% - validation dataset (validation.C: capacitance data, validation.P:
%                       pressure data)
% - model_info: struct with information regarding the model of the taxels
% - num_taxels: amount of taxels present on the insole
% - Phi_p : optionally provide capacitance regressor to estimate the
%           pressure (set it to 0 if it is to be computed inside the
%           function)
% - broken_index_valid: optionally provide list of taxels which
%                       mulfunctioned in validation (otherwise set it to 0 
%                       to compute it here)
%
% and returns the estimation for the measured pressure for each taxel, as
% well as the indexes of the mulfunctioning taxels in the dataset.
%
% author: Italo Belli (i.belli@tudelft.nl)

% extract model info
POLYNOMIAL_ORDER = model_info.POLYNOMIAL_ORDER;
HISTORY_POLYNOMIAL_ORDER = model_info.HISTORY_POLYNOMIAL_ORDER;
HISTORY_SAMPLES = model_info.HISTORY_SAMPLES;
broken_index_calib = model_info.broken_index_calib;
coefficients_all_taxels = model_info.coefficients_all_taxels;

length_coeff_per_taxel = POLYNOMIAL_ORDER + 1 + HISTORY_SAMPLES * HISTORY_POLYNOMIAL_ORDER;

if broken_index_valid==0
    % FIND the taxels that seem broken
    broken_index_valid = [];
    [validation.C, removed_index] = filter_broken_taxels(validation.C);
    if(~isempty(removed_index))
        for j=1:size(removed_index,2)
            if ~any(broken_index_valid==removed_index(j))
                broken_index_valid = [broken_index_valid, removed_index(j)];
            end
        end
    end 
end
broken_index = sort(unique([broken_index_valid, broken_index_calib]));

%% compute the regressors containing the pressure values, to perform pressure estimation
if Phi_p ==0
    disp('- Compute regressors');
    
    Phi_p = zeros(size(validation.C,1), length_coeff_per_taxel, num_taxels-length(broken_index));
    index = 1;
    for taxel = 1 : num_taxels
      if ismember(taxel,broken_index) == 0
        Phi_p(:, :, index) = regressor(validation.C(:,taxel), POLYNOMIAL_ORDER, HISTORY_SAMPLES, HISTORY_POLYNOMIAL_ORDER);
        index = index + 1;
      end
    end
end

% Use the optimal coefficients to estimate the pressure
disp('- Estimate pressures');

estimated_pressure_all = zeros(size(Phi_p,1), num_taxels);
index=1;
for taxel = 1 : num_taxels
    if ismember(taxel,broken_index) == 0
        start_index = (taxel-1)*length_coeff_per_taxel + 1;
        end_index = taxel*length_coeff_per_taxel;
        estimated_pressure_all(:, taxel) = Phi_p(:,:,index)*coefficients_all_taxels(start_index:end_index);
        index=index+1;
    else
        estimated_pressure_all(:, taxel) = zeros(size(Phi_p,1), 1);
    end
end
