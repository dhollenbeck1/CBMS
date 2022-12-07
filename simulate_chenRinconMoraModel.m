clc
clear

% -- Battery Lifetime
params.capacity = 25; % AH
params.f1 = 0.9; % cycle number correction factor
params.f2 = 0.9; % temperature correction factor
% 0 < V_SOC < 1
params.C_capacity = 3600*params.capacity*params.f1*params.f2;
params.R_selfDischarge = 1e12;

% -- Voltage-Current Characteristics 
fn.V_oc     = @(SOC,p)  p(1)*exp(p(2)*SOC) + p(3)*SOC.^3 + p(4)*SOC.^2 + p(5)*SOC + p(6);
fn.R_series = @(SOC,p)  p(1)*exp(p(2)*SOC) + p(3);
fn.R_S      = @(SOC,p)  p(1)*exp(p(2)*SOC) + p(3);
fn.C_S      = @(SOC,p)  p(1)*exp(p(2)*SOC) + p(3);
fn.R_L      = @(SOC,p)  p(1)*exp(p(2)*SOC) + p(3);
fn.C_L      = @(SOC,p)  p(1)*exp(p(2)*SOC) + p(3);

% -- Voltage-Current Characteristics
params.V_oc         = [-1.031,-35,0.3201,-0.1178,0.2156,3.685];
params.R_series     = [0.1562,-24.37,0.07446];
params.R_S          = [0.3208,-29.14,0.04669];
params.C_S          = [-752.9,-13.51,703.6];
params.R_L          = [6.603,-155.2,0.04984];
params.C_L          = [-6056,-27.12,4475];




