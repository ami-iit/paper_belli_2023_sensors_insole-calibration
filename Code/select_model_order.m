% This code was used to get the results described in the paper 
% "Modeling and Calibration of Pressure Sensing Insoles via a New
% Plenum-Based Chamber" (Belli et al., 2022)
%
% Please refer to the paper for the theoretical background, and to the
% README in this folder for the technical details on how to run the code.
%
%% Script for selecting at the same time: 
% - best polynomial order (n_p) 
% - best ARX order (n_s)
% - best order for the ARX polynomials, making it a non-linear ARX (n_{ps})

% Different values for those parameters can be tested, and the metric used
% to then be able to manually select the best combination is the RMSE on
% the validation dataset. 
% This script considers already multiple validation sets, to have a more
% complete evaluation of the performaces of each model.
%
% The execution time depends on the dimension of the datasets and on the
% parameters to be tried, depending on the number of parameters combination 
% and on the dimension of the datasets, the execution of the code can take
% quite some time (~1 hr)
%
% author: Italo Belli (i.belli@tudelft.nl)

clear
close all
clc

%% Set the correct paths
% set the path current folder to be the one where this script is contained
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

% getting path to other folders in this repo
% getting path to other folders in this repo
addpath(pathstr)
addpath('Utils\')
addpath('..\Data\left_insole\')
addpath('..\Data\results\')

% load configuration parameters
configurationfile;

% loading training set
training = load('calibration_dataset');
training = training.experiment;

% loading validation sets
validation1 = load('validation_dataset_1');
validation1 = validation1.experiment;

validation2 = load('validation_dataset_2');
validation2 = validation2.experiment;

validation3 = load('validation_dataset_3');
validation3 = validation3.experiment;

%% preprocessing (Training dataset)
% ALIGNING the dataset
training.size_exp = min(min(size(training.P, 1), size(training.C,1)));
training.C = training.C(1:training.size_exp, 3:end);
training.P = training.P(1:training.size_exp, :);

% FILTERING the dataset, in two steps:
% - use filter_high_variation_data, with a step size of 1, to get rid of
%   weird spikes in the pressure dataset;
% - use an exponential filter, to reduce the measuring noise
[training.P, training.C] = filter_high_variation_data(training.P, training.C, 0.05, 1);

alfa_c = 0.1;
alfa_p = 0.7;

for j=1:NUMBER_OF_TAXELS
    for k=2:size(training.C,1)
        training.C(k,j)=alfa_c*training.C(k,j)+(1-alfa_c)*training.C(k-1, j);
    end
end

for k=2:size(training.P, 1)
    training.P(k)=alfa_p*training.P(k)+(1-alfa_p)*training.P(k-1);
end

% FIND the taxels that seem broken, to exclude them from the calibration
broken_index_calib = [];

[training.C, removed_index] = filter_broken_taxels(training.C);
if(~isempty(removed_index))
    for j=1:size(removed_index,2)
        if ~any(broken_index_calib==removed_index(j))
            broken_index_calib = [broken_index_calib, removed_index(j)];
        end
    end
end

broken_index_calib = sort(broken_index_calib);

% set to 0 the value of the capacitance for those taxels. In this way the
% dimension of the datset remains untouched, but we can track back which
% sensors have been mulfunctioning.
training.C(:, broken_index_calib)= zeros(size( training.C(:, broken_index_calib)));

NUMBER_OF_TAXELS = size(training.C, 2);

% CHECK that, among the taxels that are working, we do not have excessive
% values for the capacitance (to account for the hysteresis effect) and
% remove possible negative values from the pressure dataset
ind = (CAPACITANCE_REST_CONDITION - training.C) < 0;
training.C(ind) = CAPACITANCE_REST_CONDITION;
ind = training.P < 0;
training.P(ind) = 0;

% PRUNE the dataset, given that its very big dimension can make our
% optimization problem very ill-conditioned

% The calibration dataset is halved for 5 times
for i=1:5
    training.P(1:2:end)=[];
    training.C(1:2:end, :)=[];
end

%% preprocessing (Validation dataset)
% ALIGNING the dataset 1
validation1.size_exp = min(size(validation1.P, 1), size(validation1.C,1));

validation1.C = validation1.C(1:validation1.size_exp, 3:end);
validation1.P = validation1.P(1:validation1.size_exp, :);

% FILTERING the dataset, in two steps:
% - use filter_high_variation_data, with a step size of 1, to get rid of
%   weird spikes in the pressure dataset;
% - use an exponential filter, to reduce the measuring noise
[validation1.P, validation1.C] = filter_high_variation_data(validation1.P, validation1.C, 0.05, 1);

for j=1:NUMBER_OF_TAXELS
    for k=2:size(validation1.C,1)
        validation1.C(k,j)=alfa_c*validation1.C(k,j)+(1-alfa_c)*validation1.C(k-1, j);
    end
end

for k=2:size(validation1.P, 1)
    validation1.P(k)=alfa_p*validation1.P(k)+(1-alfa_p)*validation1.P(k-1);
end

validation1.size_exp = size(validation1.P, 1);    % evaluating dimension after filtering

% ALIGNING the dataset 2
validation2.size_exp = min(size(validation2.P, 1), size(validation2.C,1));

validation2.C = validation2.C(1:validation2.size_exp, 3:end);
validation2.P = validation2.P(1:validation2.size_exp, :);

% FILTERING the dataset, in two steps:
% - use filter_high_variation_data, with a step size of 1, to get rid of
%   weird spikes in the pressure dataset;
% - use an exponential filter, to reduce the measuring noise
[validation2.P, validation2.C] = filter_high_variation_data(validation2.P, validation2.C, 0.05, 1);

for j=1:NUMBER_OF_TAXELS
    for k=2:size(validation2.C,1)
        validation2.C(k,j)=alfa_c*validation2.C(k,j)+(1-alfa_c)*validation2.C(k-1, j);
    end
end

for k=2:size(validation2.P, 1)
    validation2.P(k)=alfa_p*validation2.P(k)+(1-alfa_p)*validation2.P(k-1);
end

validation2.size_exp = size(validation2.P, 1);    % evaluating dimension after filtering

% ALIGNING the dataset 3
validation3.size_exp = min(size(validation3.P, 1), size(validation3.C,1));

validation3.C = validation3.C(1:validation3.size_exp, 3:end);
validation3.P = validation3.P(1:validation3.size_exp, :);

% FILTERING the dataset, in two steps:
% - use filter_high_variation_data, with a step size of 1, to get rid of
%   weird spikes in the pressure dataset;
% - use an exponential filter, to reduce the measuring noise
[validation3.P, validation3.C] = filter_high_variation_data(validation3.P, validation3.C, 0.05, 1);

for j=1:NUMBER_OF_TAXELS
    for k=2:size(validation3.C,1)
        validation3.C(k,j)=alfa_c*validation3.C(k,j)+(1-alfa_c)*validation3.C(k-1, j);
    end
end

for k=2:size(validation3.P, 1)
    validation3.P(k)=alfa_p*validation3.P(k)+(1-alfa_p)*validation3.P(k-1);
end

validation3.size_exp = size(validation3.P, 1);    % evaluating dimension after filtering

% concatenate the validation sets
validation.P = [validation1.P; validation2.P; validation3.P];
validation.C = [validation1.C; validation2.C; validation3.C];

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

broken_index_valid = sort(broken_index_valid);

% prune the validation set, due to memory errors
for i=1:4
    validation.P(1:2:end)=[];
    validation.C(1:2:end, :)=[];
    validation1.size_exp = floor(validation1.size_exp/2);
    validation2.size_exp = floor(validation2.size_exp/2);
    validation3.size_exp = floor(validation3.size_exp/2);
end

% set to 0 the value of the capacitance for those taxels. In this way the
% dimension of the datset remains untouched, but I can track back which
% sensors have been mulfunctioning.
validation.C(:, broken_index_valid)= zeros(size( validation.C(:, broken_index_valid)));

%% MODEL SELECTION phase
% max value of the parameters to be explored (USER CHOSEN)
maxOrder_pol = 6;
maxNumber_samples = 60;
step_samples = 10;
maxOrder_pol_history = 6;

rmse_calib=zeros(maxOrder_pol, maxNumber_samples/step_samples, maxOrder_pol_history);
rmse_min_calib=zeros(maxOrder_pol, maxNumber_samples/step_samples, maxOrder_pol_history);
rmse_valid = zeros(maxOrder_pol, maxNumber_samples/step_samples, maxOrder_pol_history);
rmse_min_valid = zeros(maxOrder_pol, maxNumber_samples/step_samples, maxOrder_pol_history);
disp("Performing optimization for:")
for POLYNOMIAL_ORDER = 1:maxOrder_pol
    aux = 1;
    for HISTORY_SAMPLES=0:step_samples:maxNumber_samples
        for HISTORY_POLYNOMIAL_ORDER = 1:maxOrder_pol_history
            fprintf("\n   - n_p = %i, n_s = %i, n_{ps} = %i \n", POLYNOMIAL_ORDER, HISTORY_SAMPLES, HISTORY_POLYNOMIAL_ORDER)
            % Compute regressors
            length_coeff_per_taxel = POLYNOMIAL_ORDER + 1 + HISTORY_SAMPLES * HISTORY_POLYNOMIAL_ORDER;
            Phi_p_train = zeros(size(training.C,1), length_coeff_per_taxel, NUMBER_OF_TAXELS-length(broken_index_calib));
            Phi_p_train4estim = Phi_p_train;
            Phi_p_valid = zeros(size(validation.C,1), length_coeff_per_taxel, NUMBER_OF_TAXELS-length(broken_index_calib));
            tic
            index = 1;
            scaling = [];
            for taxel = 1 : NUMBER_OF_TAXELS
                if ismember(taxel,broken_index_calib) == 0
                    % compute the regressors for the current taxel, and store it in Phi_p (both for training and validation)
                    Phi_p_train(:, :, index) = regressor(training.C(:,taxel), POLYNOMIAL_ORDER, HISTORY_SAMPLES, HISTORY_POLYNOMIAL_ORDER);
                    Phi_p_train4estim(:,:,index) = Phi_p_train(:,:,index);
                    Phi_p_valid(:,:,index) = regressor(validation.C(:,taxel), POLYNOMIAL_ORDER, HISTORY_SAMPLES, HISTORY_POLYNOMIAL_ORDER);
                    % apply a regularization to the regressor, by dividing each column by
                    % its maximum value (the optimal solution will also be treated similarly)
                    if ENABLE_REGULARIZATION==true
                        for k = 1 : length_coeff_per_taxel
                            scaling=[scaling; max(Phi_p_train(:, k, index))];
                            Phi_p_train(:, k, index) = Phi_p_train(:, k, index)/scaling(end);
                        end
                    end
                    index = index + 1;
                else
                    % Set a scaling of 0 for the coefficients of the meaningless taxels
                    scaling = [scaling; zeros(length_coeff_per_taxel, 1)];
                end
            end

            %Set up and solve opt problems
            index = 1;
            k = [];
            condition_number = zeros(NUMBER_OF_TAXELS, 2);

            index=1;
            for taxel = 1 : NUMBER_OF_TAXELS
                if ismember(taxel,broken_index_calib) == 0

                    % compute hessian and gradient
                    H_p_i = Phi_p_train(:, :, index)' * Phi_p_train(:, :, index);
                    condition_number(taxel,1) = cond(H_p_i);
                    if ENABLE_REGULARIZATION == true
                        H_p_i = H_p_i + LAMBDA_REG * eye(size(H_p_i,1));
                        condition_number(taxel,2) = cond(H_p_i);
                    end
                    g_p_i = Phi_p_train(:, :, index)' * training.P;

                    % solve optimization problem for taxel i
                    k_i = solveOptimizationProblem(H_p_i, g_p_i).x;
                    if ENABLE_REGULARIZATION==true
                        for j = 1 : length_coeff_per_taxel
                            k_i(j)=k_i(j)/scaling((taxel-1)*length_coeff_per_taxel+j);
                        end
                    end
                    index = index + 1;
                else
                    k_i = zeros(length_coeff_per_taxel, 1);
                    condition_number(taxel,:) = [0, 0];
                end
                k = [k; k_i];
            end

            toc

            % Use the optimal coefficients to estimate the pressure, on a validation
            % dataset. We need to build the regressors each time, since its parameters
            % change
            estimated_pressure_all= zeros(size(Phi_p_valid,1), NUMBER_OF_TAXELS);

            index=1;
            for taxel = 1 : NUMBER_OF_TAXELS
                if ismember(taxel,broken_index_calib) == 0
                    start_index = (taxel-1)*length_coeff_per_taxel + 1;
                    end_index = taxel*length_coeff_per_taxel;
                    estimated_pressure_all(:, taxel) = Phi_p_valid(:,:,index)*k(start_index:end_index);
                    index=index+1;
                else
                    estimated_pressure_all(:, taxel) = zeros(size(Phi_p_valid,1), 1);
                end
            end

            % evaluate the results
            rms_all_1 =[];
            rms_all_2 =[];
            rms_all_3 =[];
            rms_all=[];
            
            rms_relevant_1 = [];
            rms_relevant_2 = [];
            rms_relevant_3 = [];
            rms_relevant= [];
            for i=1:NUMBER_OF_TAXELS
                if(~any(broken_index_calib == i))
                    aux1 = sqrt(mean((validation.P(1:validation1.size_exp)-estimated_pressure_all(1:validation1.size_exp,i)).^2));
                    aux2 = sqrt(mean((validation.P(validation1.size_exp+1:validation1.size_exp+validation2.size_exp)-estimated_pressure_all(validation1.size_exp+1:validation1.size_exp+validation2.size_exp,i)).^2));
                    aux3 = sqrt(mean((validation.P(validation1.size_exp+validation2.size_exp+1:validation1.size_exp+validation2.size_exp+validation3.size_exp)-estimated_pressure_all(validation1.size_exp+validation2.size_exp+1:validation1.size_exp+validation2.size_exp+validation3.size_exp,i)).^2));
                    rms_all_1=[rms_all_1, aux1];
                    rms_relevant_1 = [rms_relevant_1, aux1];
                    rms_all_2=[rms_all_2, aux2];
                    rms_relevant_2 = [rms_relevant_2, aux2];
                    rms_all_3=[rms_all_3, aux3];
                    rms_relevant_3 = [rms_relevant_3, aux3];
                    
                    rms_all=[rms_all, mean([aux1, aux2, aux3])];
                    rms_relevant= [rms_relevant, mean([aux1, aux2, aux3])];
                else
                    rms_all = [rms_all, 0];
                end
            end

            minimum_rms = min(rms_relevant);
            maximum_rms = max(rms_relevant);
            avg_rms = mean(rms_relevant);
            lowest_rms_taxel = find(rms_all == minimum_rms);
            highest_rms_taxel = find(rms_all == maximum_rms);
            rmse_valid(aux, HISTORY_POLYNOMIAL_ORDER, POLYNOMIAL_ORDER) = avg_rms;
            rmse_min_valid(aux, HISTORY_POLYNOMIAL_ORDER, POLYNOMIAL_ORDER) = minimum_rms;

%             % Use again the optimal coefficients, to estimate the RMSE error in
%             % training set
%             estimated_pressure_all= zeros(size(Phi_p_train4estim,1), NUMBER_OF_TAXELS);
%         
%             index=1;
%             for taxel = 1 : NUMBER_OF_TAXELS
%                 if ismember(taxel,broken_index_calib) == 0
%                     start_index = (taxel-1)*length_coeff_per_taxel + 1;
%                     end_index = taxel*length_coeff_per_taxel;
%                     estimated_pressure_all(:, taxel) = Phi_p_train4estim(:,:,index)*k(start_index:end_index);
%                     index=index+1;
%                 else
%                     estimated_pressure_all(:, taxel) = zeros(size(Phi_p_train4estim,1), 1);
%                 end
%             end
%         
%             % evaluate the results
%         
%             rms_all =[];
%             rms_relevant = [];
%             for i=1:NUMBER_OF_TAXELS
%                 if(~any(broken_index_calib == i))
%                     rms_all=[rms_all, sqrt(mean((training.P-estimated_pressure_all(:,i)).^2))];
%                     rms_relevant = [rms_relevant, sqrt(mean((training.P-estimated_pressure_all(:,i)).^2))];
%                 else
%                     rms_all = [rms_all, 0];
%                 end
%             end
%         
%             minimum_rms = min(rms_relevant);
%             maximum_rms = max(rms_relevant);
%             avg_rms = mean(rms_relevant);
%             lowest_rms_taxel = find(rms_all == minimum_rms);
%             highest_rms_taxel = find(rms_all == maximum_rms);
%             rmse_calib(aux, HISTORY_POLYNOMIAL_ORDER, POLYNOMIAL_ORDER) = avg_rms;
%             rmse_min_calib(aux, HISTORY_POLYNOMIAL_ORDER, POLYNOMIAL_ORDER) = minimum_rms;
        end
        aux = aux+1;
    end
    disp(POLYNOMIAL_ORDER);
end

%% PLOT THE MESH REPRESENTING PERFORMANCES 
% On the vertical axis we have the RMSE, and each plane represents a
% different POLINOMIAL_ORDER
[X,Y] = meshgrid(0:step_samples:maxNumber_samples , 1:maxOrder_pol_history);

figure
hold on
for i=1:maxOrder_pol
    mesh(X, Y, rmse_valid(1:end,1:end,i)', 'FaceColor', 'blue')
end

grid on
xlabel('Past samples ARX (n_s)')
ylabel('Polynomial order ARX part (n_{ps})')
zlabel('RMSE [bar]')
hold off