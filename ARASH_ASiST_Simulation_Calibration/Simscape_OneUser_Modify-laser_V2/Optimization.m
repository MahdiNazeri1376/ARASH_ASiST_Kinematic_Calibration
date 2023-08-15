clc 
clear 
close all

% Load the dataset containing measured values (px, py, pz) and ground truth data
 dataset = load("dataset3_for_calibration_fix_d(FK_GT).mat"); % Replace 'dataset.mat' with the actual filename

% Define the objective function (error metric)
objectiveFcn = @(params) calculateError(params, dataset);

% Define the parameter bounds
lb = [30*(pi/180), 0.2, 30*(pi/180)]; % Lower bounds for theta, L, alpha
ub = [60*(pi/180), 0.5, 60*(pi/180)]; % Upper bounds for theta, L, alpha


% %%%%%%%%%%%%%%fmincon%%%%%%%%%%%%%%%%%%
% initialGuess = [31*(pi/180), 0.29, 31*(pi/180)]; % Initial guess for theta, L, alpha
% % Run the optimization using fmincon
% options = optimoptions('fmincon', 'Display', 'iter');
% optimizedParams = fmincon(objectiveFcn, initialGuess, [], [], [], [], lb, ub, [], options);
% 

%%%%%%%%%%%% Genetic Algorithm (GA): %%%%%%%%%%%%%%%%%%%%%
% Set up GA options
gaOptions = optimoptions('ga', 'Display', 'iter');

% Run GA optimization
optimizedParams = ga(objectiveFcn, 3, [], [], [], [], lb, ub, [], gaOptions);


% %%%%%%%%%%%%% Particle Swarm Optimization (PSO) %%%%%%%%%%%%%%%%%%%%%
% % Set up PSO options
% psoOptions = optimoptions('particleswarm', 'Display', 'iter');
% 
% % Run PSO optimization
% optimizedParams = particleswarm(objectiveFcn, 3, lb, ub, psoOptions);


% %%%%%%%%%%%%% Pattern Search %%%%%%%%%%%%%%%%%%%%%
% % Set up pattern search options
% initialGuess = [31*(pi/180), 0.29, 31*(pi/180)]; % Initial guess for theta, L, alpha
% patternSearchOptions = optimoptions('patternsearch', 'Display', 'iter');
% 
% % Run pattern search optimization
% optimizedParams = patternsearch(objectiveFcn, initialGuess, [], [], [], [], lb, ub, patternSearchOptions);


% %%%%%%%%%%%% Trust-Region Reflective Algorithm: %%%%%%%%%%%%%%%%%%%%%
% % Set up trust-region reflective algorithm options
% initialGuess = [31*(pi/180), 0.29, 31*(pi/180)]; % Initial guess for theta, L, alpha
% trustRegionOptions = optimoptions('lsqnonlin', 'Algorithm', 'trust-region-reflective', 'Display', 'iter');
% 
% % Run trust-region reflective algorithm optimization
% optimizedParams = lsqnonlin(objectiveFcn, initialGuess, lb, ub, trustRegionOptions);
% 

% Extract the optimized parameters
thetaOptimized = optimizedParams(1);
LOptimized = optimizedParams(2);
alphaOptimized = optimizedParams(3);

% Function to calculate the error metric
function error = calculateError(params, dataset)
    theta = params(1);
    L = params(2);
    alpha = params(3);
    
    pxMeasured = dataset.x_GT;
    pyMeasured = dataset.y_GT;
    pzMeasured = dataset.z_GT;
    phiMeasured = dataset.phi;
    psiMeasured = dataset.psi;
    dMeasured   = dataset.d;
    
    % Calculate predicted values using the forward kinematic equations
    pxPredicted = (sin(theta).*cos(psiMeasured + alpha).*dMeasured + cos(theta).*cos(phiMeasured).*sin(psiMeasured + alpha).*dMeasured + L*sin(theta))*1000;
    pyPredicted = (dMeasured.*sin(phiMeasured).*sin(psiMeasured + alpha))*1000;
    pzPredicted = (-sin(theta).*cos(phiMeasured).*sin(psiMeasured + alpha).*dMeasured + cos(theta).*cos(psiMeasured + alpha).*dMeasured + L*cos(theta))*1000;
    

    % Calculate the error metric (e.g., mean squared error)
%     error = mse((pxMeasured - pxPredicted).^2 + (pyMeasured - pyPredicted).^2 + (pzMeasured - pzPredicted).^2);
error = 10000*rms((pxMeasured - pxPredicted).^2 + (pyMeasured - pyPredicted).^2 + (pzMeasured - pzPredicted).^2);
end