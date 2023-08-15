clc
clear

load('dataset4_for_calibration_fix_d(RCM_GT).mat');

C = [d, phi,psi,x_GT,y_GT,z_GT];
%Write CSV file
csvwrite('dataset4_for_calibration_fix_d(RCM_GT).csv',C);