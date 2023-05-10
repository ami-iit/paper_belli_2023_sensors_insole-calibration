% This code was used to get the results described in the paper 
% "Modeling and Calibration of Pressure Sensing Insoles via a New
% Plenum-Based Chamber" (Belli et al., 2023)
%
% Please refer to the paper for the theoretical background, and to the
% README in this folder for the technical details on how to run the code.
%
% This code implements the validation step: after running a calibration as
% in `calibration_taxels.m` and saving the results correctly (this is
% implemented in the code), we visualize here how the model found performs
% on the validation scenarios considered.

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
addpath(pathstr)
addpath('Utils\')
addpath('..\Data\')
addpath('..\Data\results\')

configurationfile;

% loading validation datasets (from Data folder)
validation1 = load('validation_dataset_1');
validation1 = validation1.experiment;

validation2 = load('validation_dataset_2');
validation2 = validation2.experiment;

validation3 = load('validation_dataset_3');
validation3 = validation3.experiment;

% load the parameters found in calibration
% This loads:
% - the various polynomial orders for the general model of the taxels
% - the optimized coefficients for each taxel
% - the taxels which were found to be broken/defectuous in calibration
% - the parameters used in the exponential filter during the processing
%   phase of the data in calibration (alfa_p, alfa_c for pressure and capacitance)
model_info = load('results\model_20221108.mat');
alfa_c = model_info.alfa_c;
alfa_p = model_info.alfa_p;
broken_index_calib = model_info.broken_index_calib;


%% PREPROCESSING
eps = 0.05;     % max variation allowed from one sample to the next (in bar)
step = 1;       % step at which the above variation is enforced (step=1= means all the data is checked for it)

% ALIGNING the dataset 1
validation1.size_exp = min(size(validation1.P, 1), size(validation1.C,1));

validation1.C = validation1.C(1:validation1.size_exp, 3:end);   % capacitance values for the experiment
validation1.P = validation1.P(1:validation1.size_exp, :);       % pressure values for the experiment

% FILTERING the dataset, in two steps:
% - use filter_high_variation_data, with a step size of 1, to get rid of
%   weird spikes in the pressure dataset;
% - use an exponential filter, to reduce the measuring noise
[validation1.P, validation1.C] = filter_high_variation_data(validation1.P, validation1.C, eps, step);

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
[validation2.P, validation2.C] = filter_high_variation_data(validation2.P, validation2.C, eps,step);

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
[validation3.P, validation3.C] = filter_high_variation_data(validation3.P, validation3.C, eps, step);

for j=1:NUMBER_OF_TAXELS
    for k=2:size(validation3.C,1)
        validation3.C(k,j)=alfa_c*validation3.C(k,j)+(1-alfa_c)*validation3.C(k-1, j);
    end
end

for k=2:size(validation3.P, 1)
    validation3.P(k)=alfa_p*validation3.P(k)+(1-alfa_p)*validation3.P(k-1);
end

validation3.size_exp = size(validation3.P, 1);    % evaluating dimension after filtering

%% ESTIMATING PRESSURE
disp('First dataset')
[estimated_pressure_1, broken_index_valid_1] = estimate_pressure(validation1, model_info, NUMBER_OF_TAXELS, 0, 0);

disp('Second dataset')
[estimated_pressure_2, broken_index_valid_2] = estimate_pressure(validation2, model_info, NUMBER_OF_TAXELS, 0, 0);

disp('Third dataset')
[estimated_pressure_3, broken_index_valid_3] = estimate_pressure(validation3, model_info, NUMBER_OF_TAXELS, 0, 0);

%% COMPARE ESTIMATION WITH GROUND TRUTH
% list all the taxels that did not work in calibration or at least one of
% the validation. They will be excluded as representing data coming from
% mulfunctioning sensors
broken_index = sort(unique([broken_index_calib, broken_index_valid_1, broken_index_valid_2, broken_index_valid_3]));

% concatenate the trhee datasets (ground-truth and estimations)
validation.P = [validation1.P; validation2.P; validation3.P];
validation.C = [validation1.C; validation2.C; validation3.C];
estimated_pressure_all = [estimated_pressure_1; estimated_pressure_2; estimated_pressure_3];

% evaluate root mean square errors (rmse)
rmse = zeros(NUMBER_OF_TAXELS,1);

for i=1:NUMBER_OF_TAXELS
    if ~any(broken_index == i)
        rmse(i) = sqrt(mean((validation.P-estimated_pressure_all(:,i)).^2));
    end
end

minimum_rms = min(rmse(rmse>0));
maximum_rms = max(rmse);
avg_rms = mean(rmse(rmse>0));

lowest_rms_taxel = find(rmse == minimum_rms);
highest_rms_taxel = find(rmse == maximum_rms);

%% Plotting
timestamps = (1:size(validation.P,1))/ACQUISITION_FREQ;

% visualize estimates
figure
hold on
plot(timestamps, validation.P, 'LineWidth', 1.3, 'color', 'black')
plot(timestamps, estimated_pressure_all(:, lowest_rms_taxel), 'LineWidth', 0.9)
plot(timestamps, estimated_pressure_all(:, highest_rms_taxel), 'LineWidth', 0.9, 'color', '#4E8009')
xlabel(" time [s]")
ylabel("pressure [bar]")
lowest_legend = "estimate best taxel";
worst_legend = "estimate worst taxel";
legend("measured pressure", lowest_legend, worst_legend);
title("Validation results");
grid on
hold off

%% 
% visualize all estimates
figure
hold on
plot(timestamps, estimated_pressure_all(:,:))
plot(timestamps, validation.P, '*')
xlabel(" time [s]")
ylabel("pressure [bar]")
title("Validation results");
grid on
hold off

% visualize some estimates
figure
hold on
plot(timestamps, validation.P, '*')
plot(timestamps, estimated_pressure_all(:,1:10))
xlabel(" time [s]")
ylabel("pressure [bar]")
title("Validation results");
grid on
hold off

%% print useful metrics
avg_rms
minimum_rms
lowest_rms_taxel
maximum_rms
highest_rms_taxel
