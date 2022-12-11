clc
clear

% -- NOTES
% state of charge varies between 0 < V_SOC < 1
T = 1000;
dt = 1e-1;
n = round(T/dt);
% t = linspace(0,T,n)';
% R_load = 0.4; %[t,0.4*ones(n,1)];

% -- Battery Lifetime
p.capacity = 25; % AH
p.f1 = 0.9; % cycle number correction factor
p.f2 = 0.9; % temperature correction factor
p.C_capacity = 3600*p.capacity*p.f1*p.f2;
p.R_selfDischarge = 1e9;

% -- Voltage-Current Characteristics 
fn.V_oc     = @(SOC,p)  p(1)*exp(p(2)*SOC) + p(3)*SOC.^3 + p(4)*SOC.^2 + p(5)*SOC + p(6);
fn.R_series = @(SOC,p)  p(1).*exp(p(2).*SOC) + p(3);
fn.R_1      = @(SOC,p)  p(1).*exp(p(2).*SOC) + p(3);
fn.C_1      = @(SOC,p)  p(1).*exp(p(2).*SOC) + p(3);
fn.R_2      = @(SOC,p)  p(1).*exp(p(2).*SOC) + p(3);
fn.C_2      = @(SOC,p)  p(1).*exp(p(2).*SOC) + p(3);

% -- Voltage-Current Characteristics
p.V_oc         = 12*[-1.031,-35,0.3201,-0.1178,0.2156,3.685];
p.R_series     = [0.1562,-24.37,0.07446];
p.R_1          = [0.3208,-29.14,0.04669];
p.C_1          = [-752.9,-13.51,703.6];
p.R_2          = [6.603,-155.2,0.04984];
p.C_2          = [-6056,-27.12,4475];

% -- Variables
v.SOC         = 1;
v.V_oc        = fn.V_oc(v.SOC,p.V_oc);
v.R_series    = fn.R_series(v.SOC,p.R_series);
v.R_1         = fn.R_1(v.SOC,p.R_1);
v.C_1         = fn.C_1(v.SOC,p.C_1);
v.R_2         = fn.R_2(v.SOC,p.R_2);
v.C_2         = fn.C_2(v.SOC,p.C_2);

% sim("chenRinconMoraModel.slx")
