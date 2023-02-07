% clc
% clear

addpath data\ data_cut\ data_raw\ auxillary\

% -- Import Data
fn = 'data_cut\data_20220810.mat';
load(fn)
T = data.te(end);
dt = 1e-0;
n = round(T/dt);
t = linspace(0,T,n)';
data.ic_interp = interp1(data.te,data.ie,t);
data.vc_interp = interp1(data.te,data.ve,t);
data.rc_interp = interp1(data.te,data.re,t);
data.rc_interp = movmean(data.rc_interp,5);
% data.rc_interp = sgolayfilt(data.rc_interp,13,21);
rclimit = 1;
data.rc_interp(data.rc_interp<rclimit) = rclimit;

% -- Input to Simulink model
simin.time = t;
simin.signals.values = data.rc_interp;

exp_simin_v = timeseries(data.vc_interp,t);
exp_simin_i = timeseries(data.ic_interp,t);
% exp_simin_v.time = t;
% exp_simin_v.values = data.vc_interp;
% exp_simin_i.time = t;
% exp_simin_i.values = data.ic_interp;
% simin.signals.dimensions = size(data.x(:,2));
% R_load = 0.4; %[t,0.4*ones(n,1)];

% -- Battery Lifetime
p.capacity = 25; % AH
p.f1 = 0.9; % cycle number correction factor
p.f2 = 0.9; % temperature correction factor
p.C_capacity = 3600*p.capacity*p.f1*p.f2;
p.R_selfDischarge = 1e12;

% -- fitting functions
f       = @(SOC,p)  p(1)*exp(p(2)*SOC) + p(3)*SOC.^3 + p(4)*SOC.^2 + p(5)*SOC + p(6);
g       = @(SOC,p)  p(1).*exp(p(2).*SOC) + p(3);
% gradf_p   = @(SOC,p)
% gradf_SOC   = @(SOC,p)

% -- Voltage-Current hyperparams
p.R_series     = [1.63339628126812e-05          8.71313419110223     -0.000181956163055476];
p.R_1          = [-2.72130383416836e-07          12.4043553051855        0.0337869742167239];
p.C_1          = [                0                         0                        89];
p.R_2          = [0.021478839277293          0.19205191547777       -0.0115329708457113];
p.C_2          = [0                         0          6091.58194788576];
p.R_3          = 1.5*[6.84750107001574        0.0741110878565304         -6.96796134105061];
p.C_3          = [ 0*14215582.5672873+1e8        -0.558664045080813        -5e5+0*-8929606.37707713];
p.K_V_oc         = [0.153418165249012+0.05         0.273499935298638-0.02         0.835830972401159-0.05];

% -- Voltage-Current Characteristics from fitted Manufacturer data
ncells = 12;  % Number of cells used: 6S x2
p.V_oc         = ncells*[-0.6388,-48.7497,0.7153,-0.6260,0.4743,3.6115];

% -- Voltage-Current Characteristics from Chen-Rincon-Mora
% ncells = 12;
% p.V_oc         = ncells*[-1.031,-30,0.3201,-0.1178,0.2156,3.685];
% p.R_series     = [0.1562,-24.37,0.07446];
% p.R_1          = 0*[0.3208,-29.14,0.04669];
% p.C_1          = 0*[-600,-13.51,703.6];
% p.R_2          = 0*[6.603,-155.2,0.04984];
% p.C_2          = 0*[-3056,-27.12,4475];

% -- Get initial SOC
voc_0       = data.vc_interp(1)-0*0.7*ncells;  % percent*usableVolt/cell*ncells
SOC         = linspace(0,1,1000);
temp        = f(SOC,p.V_oc);
tempk       = find(temp <= voc_0,1,"last");
SOC_0       = SOC(tempk);

% -- Initialize Simulation Variables
v.SOC         = 0.86;%SOC_0;
v.V_oc        = f(v.SOC,p.V_oc);
v.R_series    = g(v.SOC,p.R_series);
v.R_1         = g(v.SOC,p.R_1);
v.C_1         = g(v.SOC,p.C_1);
v.R_2         = g(v.SOC,p.R_2);
v.C_2         = g(v.SOC,p.C_2);

%% Loop
format long g
dJ = 0;
J = 0;
ep = -10.^[-3,-2,-2,-2,-2];
nk = 1;
% R_series = linspace(0.01,0.07,nk);
tol = 1e-4;
for k=1:nk
    %     clc
    %     disp(k)
    %     p.R_series(3) = R_series(k);
    %out = sim("chenRinconMoraModel_paramEstimation.slx");
    out = sim("CRM_model.slx");
    %
    % -- Extract data...
    y.t = out.simout.Time;% 0: time
    y.R_series = out.simout.Data(:,1);% 1: R_series
    y.R_1 = out.simout.Data(:,2);% 2: R_1
    y.C_1 = out.simout.Data(:,3);% 3: C_1
    y.R_2 = out.simout.Data(:,4);% 4: R_2
    y.C_2 = out.simout.Data(:,5);% 5: C_2
    y.SOC = out.simout.Data(:,6);% 6: SOC
    y.V_oc = out.simout.Data(:,7);% 7: Voc
    y.i = out.simout.Data(:,8);  % 8: i_bat
    y.V_L = out.simout.Data(:,9);% 9: V_L
    y.R_L = simin.signals.values;%10: R_L

    ylen = length(y.R_series);
    datalen = length(data.vc_interp);
    J(k) = sum((y.V_L(2:end)-data.vc_interp(1:ylen-1)).^2 +(y.i(2:end)-data.ic_interp(1:ylen-1)).^2);
    thet = [p.R_series(3),p.R_1(3),p.C_1(3),p.R_2(3),p.C_2(3)];
%     lam = 1e-4;
%     J(k) = J(k) + lam*sum(thet);
    if k==1
        dJ = ep(1).^2*randn();
    else
        dJ = J(k)-J(k-1);
    end

    if dJ == 0
        dJ = 0*0.01.^2*randn();
    end

        %     param_dir = randi(5);
        dk = mod(k,1)+1;
        switch dk
            case 1
                %randi(3);
                temp = p.R_series(3) + ep(dk)*dJ;
                            if temp > 0
                                p.R_series(3) = temp;
                            end
            case 2
                %randi(3);
                temp = p.R_1(3) + ep(dk)*dJ;
                            if temp > 0
                                p.R_1(3) = temp;
                            end
            case 3
                 %randi(3);
                temp = p.C_1(3) + ep(dk)*dJ;
                            if temp > 0
                                p.C_1(3) = temp;
                            end
            case 4
                 %randi(3);
                temp = p.R_2(3) + ep(dk)*dJ;
                            if temp > 0
                                p.R_2(3) = temp;
                            end
            case 5
                 %randi(3);
                temp = p.C_2(3) + ep(dk)*dJ;
                            if temp > 0
                                p.C_2(3) = temp;
                            end
            otherwise
        end
        clc
        thet = [p.R_series(3),p.R_1(3),p.C_1(3),p.R_2(3),p.C_2(3)];
        disp([k,thet,dJ]')
%         disp()
%         disp(p.R_1)
%         disp(p.C_1)
%         disp(p.R_2)
%         disp(p.C_2)
%         disp(dJ)
        % disp(param_dir)

%     end
    % Plots
    % -- Plot data...
    lw = 2;
    plot(data.tc,data.vc,'k','LineWidth',lw)
    hold on;
    plot(y.t(2:end),y.V_L(2:end),'r--','LineWidth',lw)
    yyaxis left
    ylim([40,52])
    hold off;
    yyaxis right
    plot(data.tc,data.ic,'b','LineWidth',lw)
    hold on
    plot(y.t(2:end),y.i(2:end),'g--','LineWidth',lw);
    hold off

    xlim([y.t(1),y.t(end)])

    legend('V_L experiment','V_L simulated','i_L experiment','i_L simulated')
    shg
end
% yyaxis right
% plot(simin.time,simin.signals.values,'g')
% ylim([0,8])
% set(gca,'YColor','g')
