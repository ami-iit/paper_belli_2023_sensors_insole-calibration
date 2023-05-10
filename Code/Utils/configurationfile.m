% Configuration parameters for calibration


% Dependent on the hardware
ACQUISITION_FREQ = 50; %Hz
CAPACITANCE_REST_CONDITION = 250;
NUMBER_OF_TAXELS = 280;
AREA_TAXEL = (0.024 * 0.027 / 2) / 10;
NUM_TRIANGLES = 28;

% Regularization
ENABLE_REGULARIZATION = true;
LAMBDA_REG = 1e-3; 