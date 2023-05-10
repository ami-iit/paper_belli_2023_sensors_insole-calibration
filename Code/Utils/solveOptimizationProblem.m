function k = solveOptimizationProblem(H, g)

% Load configuration parameters
configurationfile;

lb = [];
ub = [];
A = [];
lbA = [];
ubA = [];

% Create an OSQP object
prob = osqp;

% Setup workspace
prob.setup(H, -g, A, lbA, ubA, 'max_iter', 15000, 'verbose', 0);

% solve the problem
k = prob.solve();
end