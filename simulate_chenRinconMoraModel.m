clc
clear

% -- NOTES
% state of charge varies between 0 < V_SOC < 1
T = 1000;
dt = 1e-1;
n = round(T/dt);
t = linspace(0,T,n)';
load data_test1.mat
[C,IA,IC] =unique(data.t);
data.tc = C;
% data.tc = data.t(IC);
data.ic = data.x(IA,2);
data.vc = data.x(IA,1);
data.rc = data.vc./data.ic;
data.rc(isinf(data.rc)) = 1e9;

data.rc_interp = interp1(data.tc,data.rc,t);
simin.time = t;
simin.signals.values = data.rc_interp;
% simin.signals.dimensions = size(data.x(:,2));
% R_load = 0.4; %[t,0.4*ones(n,1)];

% -- Battery Lifetime
p.capacity = 25; % AH
p.f1 = 0.9; % cycle number correction factor
p.f2 = 0.9; % temperature correction factor
p.C_capacity = 3600*p.capacity*p.f1*p.f2;
p.R_selfDischarge = 1e9;

% ---- From Chen Rincon-Mora 
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

% % ---- From Our project 
% % -- Voltage-Current Characteristics 
% fn.V_oc     = @(SOC,p)  p(1)*exp(p(2)*SOC) + p(3)*SOC.^3 + p(4)*SOC.^2 + p(5)*SOC + p(6);
% fn.R_series = @(SOC,p)  p(1).*exp(p(2).*SOC) + p(3);
% fn.R_1      = @(SOC,p)  p(1).*exp(p(2).*SOC) + p(3);
% fn.C_1      = @(SOC,p)  p(1).*exp(p(2).*SOC) + p(3);
% fn.R_2      = @(SOC,p)  p(1).*exp(p(2).*SOC) + p(3);
% fn.C_2      = @(SOC,p)  p(1).*exp(p(2).*SOC) + p(3);
% % -- Voltage-Current Characteristics
% p.V_oc         = 12*[-1.031,-35,0.3201,-0.1178,0.2156,3.685];
% p.R_series     = [0.1562,-24.37,0.07446];
% p.R_1          = [0.3208,-29.14,0.04669];
% p.C_1          = [-752.9,-13.51,703.6];
% p.R_2          = [6.603,-155.2,0.04984];
% p.C_2          = [-6056,-27.12,4475];


% -- Variables
v.SOC         = 1;
v.V_oc        = fn.V_oc(v.SOC,p.V_oc);
v.R_series    = fn.R_series(v.SOC,p.R_series);
v.R_1         = fn.R_1(v.SOC,p.R_1);
v.C_1         = fn.C_1(v.SOC,p.C_1);
v.R_2         = fn.R_2(v.SOC,p.R_2);
v.C_2         = fn.C_2(v.SOC,p.C_2);

%% Loop
out = sim("chenRinconMoraModel_paramEstimation.slx");

% -- Extract data...
y.t = out.simout.Time;% 0: time
y.R_series = out.simout.Data(:,1);% 1: R_series
y.R_1 = out.simout.Data(:,2);% 2: R_1
y.C_1 = out.simout.Data(:,3);% 3: C_1
y.R_2 = out.simout.Data(:,4);% 4: R_2
y.C_2 = out.simout.Data(:,5);% 5: C_2
y.SOC = out.simout.Data(:,6);% 6: SOC
y.V_oc = out.simout.Data(:,7);% 7: Voc
y.i = out.simout.Data(:,8);% 8: i_bat
y.V_L = out.simout.Data(:,9);% 9: V_L
y.R_L = out.simout.Data(:,10);%10: R_L

% -- Plot data...
plot(y.t,y.V_oc,'k-')
hold on; plot(y.t,y.i,'b-'); hold off
yyaxis right
plot(y.t,y.R_L,'r-')
set(gca,'YColor','r')
