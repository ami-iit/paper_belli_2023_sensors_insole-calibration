% This code was used to get the results described in the paper 
% "Modeling and Calibration of Pressure Sensing Insoles via a New
% Plenum-Based Chamber" (Belli et al., 2023)
%
% Please refer to the paper for the theoretical background, and to the
% README in this folder for the technical details on how to run the code.
%
% Here we deal with the calibration step of our hardware, by solving a 
% different optimization problem for each taxel. 
% Compared to the work in Sorrentino et al., 2020, this method is much faster 
% and the results are almost the same (when the same model class is
% considered).
%
% After having found (or supposed) a model, in terms of:
% - n_p      (polynomial order - POLYNOMIAL_ORDER)
% - n_s      (number of past samples considered - HISTORY_SAMPLES)
% - n_{ps}   (order of the history terms - HISTORY_POLINOMIAL_ORDER)
% you can use this script to find the optimal coefficients of the model, 
% to estimate the pressure.
%
% NOTE that the code takes as input an 'experiment struct', having two fields:
% - experiment.C = capacitance history for the taxels (282 colums)
% - experiment.P = recorded pressure
%
% author: Italo Belli (i.belli@tudelft.nl)


clear
close all
clc

model_name = 'model_20221108';

%% Set the correct paths
% set the path current folder to be the one where this script is contained
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

% getting path to other folders in this repo
addpath(pathstr)
addpath('Utils\')
addpath('..\Data\results\')

%% SELECT THE INSOLE CONSIDERED
% Note: the paper is based on results from the left insole

% if you consider the right insole
% shoe = 1; 
% addpath('..\Data\right_insole\');

% if you consider the left insole
shoe = 0; 
addpath('..\Data\left_insole\');

% load parameters
configurationfile;

% load calibration dataset data
load('calibration_dataset');

%% preprocessing
% ALIGNING the dataset, by cropping the longest one
size_exp = min(min(size(experiment.P, 1), size(experiment.C,1)));
experiment.C = experiment.C(1:size_exp, 3:end);       % here the first two columns are discarded (message IDs, timestamps)
experiment.P = experiment.P(1:size_exp, :);

% FILTERING the dataset, in two steps:
% - use filter_high_variation_data, with a step size of 1, to get rid of
%   weird spikes in the pressure dataset;
% - use an exponential filter, to reduce the measuring noise
[experiment.P, experiment.C] = filter_high_variation_data(experiment.P, experiment.C, 0.05, 1);

alfa_c = 0.1;
alfa_p = 0.7;   

for j=1:NUMBER_OF_TAXELS
    for k=2:size(experiment.C,1)
        experiment.C(k,j)=alfa_c*experiment.C(k,j)+(1-alfa_c)*experiment.C(k-1, j);
    end
end

for k=2:size(experiment.P, 1)
    experiment.P(k)=alfa_p*experiment.P(k)+(1-alfa_p)*experiment.P(k-1);
end

% FIND the taxels that seem broken, to exclude them from the calibration
broken_index = [];

[experiment.C, removed_index] = filter_broken_taxels(experiment.C);
if(~isempty(removed_index))
    for j=1:size(removed_index,2)   
        if ~any(broken_index==removed_index(j))
            broken_index = [broken_index, removed_index(j)];
        end
    end
end
broken_index = sort(broken_index);

% set to 0 the value of the capacitance for those taxels. In this way the
% dimension of the dataset remains untouched, but we can track back which
% sensors have been mulfunctioning for postprocessing, if needed
experiment.C(:, broken_index)= zeros(size( experiment.C(:, broken_index)));

% CHECK that, among the taxels that are working, we do not have excessive
% values for the capacitance (to account for the hysteresis effect) and
% remove possible negative values from the pressure dataset
ind = (CAPACITANCE_REST_CONDITION - experiment.C) < 0;
experiment.C(ind) = CAPACITANCE_REST_CONDITION;
ind = experiment.P < 0;
experiment.P(ind) = 0;

% PRUNE the dataset, given that its very big dimension can make our
% optimization problem very ill-conditioned

% for 5 times, the dataset is halved, in order to reduce the condition
% number of the regression matrix built at optimization stage
for i=1:5
    experiment.P(1:2:end)=[];
    experiment.C(1:2:end, :)=[];
end


%% CALIBRATION phase
% Select the model to be calibrated, to be chosen by the USER!
POLYNOMIAL_ORDER = 3;
HISTORY_POLYNOMIAL_ORDER=4;
HISTORY_SAMPLES=40;
fprintf('Calibrating model with n_s = %i, n_s = %i, n_{ps} = %i \n\n', POLYNOMIAL_ORDER, HISTORY_SAMPLES, HISTORY_POLYNOMIAL_ORDER)

rmse=[];
min_rmse=[];
disp('Compute regressors');
length_coeff_per_taxel = POLYNOMIAL_ORDER + 1 + HISTORY_SAMPLES * HISTORY_POLYNOMIAL_ORDER;


Phi_p_temp = zeros(size(experiment.C,1), length_coeff_per_taxel, NUMBER_OF_TAXELS-length(broken_index));
Phi_p_4estim = Phi_p_temp;
tic
index = 1;
scaling = [];
for taxel = 1 : NUMBER_OF_TAXELS
    if ismember(taxel,broken_index) == 0
        % compute the regressors for the current taxel, and store it in Phi_p_temp
        Phi_p_temp(:, :, index) = regressor(experiment.C(:,taxel), POLYNOMIAL_ORDER, HISTORY_SAMPLES, HISTORY_POLYNOMIAL_ORDER);
        Phi_p_4estim(:,:,index) = Phi_p_temp(:, :, index);
        % apply a regularization to the regressor, by dividing each column by
        % its maximum value (the optimal solution will also be treated similarly)
        if ENABLE_REGULARIZATION==true
            for k = 1 : length_coeff_per_taxel
                scaling=[scaling; max(Phi_p_temp(:, k, index))];
                Phi_p_temp(:, k, index) = Phi_p_temp(:, k, index)/scaling(end);
            end
        end
        index = index + 1;
    else
        % Set a scaling of 0 for the coefficients of the meaningless taxels
        scaling = [scaling; zeros(length_coeff_per_taxel, 1)];
    end
end

disp('Set up and solve opt problems');

index = 1;
k = [];
condition_number = zeros(NUMBER_OF_TAXELS, 2);

index=1;
for taxel = 1 : NUMBER_OF_TAXELS
    if ismember(taxel,broken_index) == 0

        % compute hessian and gradient
        H_p_i = Phi_p_temp(:, :, index)' * Phi_p_temp(:, :, index);
        condition_number(taxel,1) = cond(H_p_i);
        if ENABLE_REGULARIZATION == true
            H_p_i = H_p_i + LAMBDA_REG * eye(size(H_p_i,1));
            condition_number(taxel,2) = cond(H_p_i);
        end
        g_p_i = Phi_p_temp(:, :, index)' * experiment.P;

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

tOptim = toc;

% Use the optimal coefficients to estimate the pressure 
% NOTE that Phi_p_4estim should be used here, as it is non-scaled
disp('Estimate pressures');

estimated_pressure_all= zeros(size(Phi_p_4estim,1), NUMBER_OF_TAXELS);

index=1;
for taxel = 1 : NUMBER_OF_TAXELS
    if ismember(taxel,broken_index) == 0
        start_index = (taxel-1)*length_coeff_per_taxel + 1;
        end_index = taxel*length_coeff_per_taxel;
        estimated_pressure_all(:, taxel) = Phi_p_4estim(:,:,index)*k(start_index:end_index);
        index=index+1;
    else
        estimated_pressure_all(:, taxel) = zeros(size(Phi_p_4estim,1), 1);
    end
end

%% evaluate the results
rms_all =[];
rms_relevant = [];
for i=1:NUMBER_OF_TAXELS
    if(~any(broken_index == i))
        rms_all=[rms_all, sqrt(mean((experiment.P-estimated_pressure_all(:,i)).^2))];
        rms_relevant = [rms_relevant, sqrt(mean((experiment.P-estimated_pressure_all(:,i)).^2))];
    else
        rms_all = [rms_all, 0];
    end
end

minimum_rms = min(rms_relevant);
maximum_rms = max(rms_relevant);
avg_rms = mean(rms_relevant);
lowest_rms_taxel = find(rms_all == minimum_rms);
highest_rms_taxel = find(rms_all == maximum_rms);
    
%% visualize estimates
timestamps = (1:size(experiment.P,1))/ACQUISITION_FREQ*32;

figure
hold on
plot(timestamps, experiment.P)
plot(timestamps, estimated_pressure_all(:,lowest_rms_taxel))
plot(timestamps, estimated_pressure_all(:, highest_rms_taxel), 'color', '#4E8009')
low_label = "best estimate (taxel " + lowest_rms_taxel + ")";
high_label = "worst estimate (taxel " + highest_rms_taxel + ")";
legend('measured pressure', low_label, high_label);
xlabel("time [s]")
ylabel("pressure [bar]")
title('Results on training dataset')
legend
grid on

avg_condition_number_before= mean(condition_number(:,1));
avg_condition_number_after= mean(condition_number(:,2));
%% saving the results of the calibration, to validate them elsewhere
% it is important that also the index_broken is saved, so that at
% validation stage it can be taken into account
coefficients_all_taxels = k;
broken_index_calib = broken_index;

save(model_name, 'coefficients_all_taxels', 'broken_index_calib', 'alfa_p', 'alfa_c', 'POLYNOMIAL_ORDER', 'HISTORY_POLYNOMIAL_ORDER', 'HISTORY_SAMPLES', 'shoe');

%% print useful metrics
avg_rms
minimum_rms
lowest_rms_taxel
maximum_rms
highest_rms_taxel