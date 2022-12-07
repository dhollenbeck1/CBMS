clc
clear

% --
capacity = 25; % AH
f1 = 0.9; % cycle number correction factor
f2 = 0.9; % temperature correction factor
% 0 < V_SOC < 1
C_capacity = 3600*capacity*f1*f2;
R_selfDischarge = 1e12;


R_series = 1;


C_capacity = 3600*capacity*f1*f2;