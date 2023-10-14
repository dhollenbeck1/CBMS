%
clear all
clc
% Copyright 2016-2018 The MathWorks, Inc.
%
% Using MathWorks tools, estimation techniques, and measured lithium-ion 
% or lead acid battery data, you can generate parameters for the 
% Equivalent Circuit Battery block. The Equivalent Circuit Battery 
% block implements a resistor-capacitor (RC) circuit battery with open 
% circuit voltage, series resistance, and 1 through N RC pairs. The 
% number of RC pairs reflects the number of time constants that 
% characterize the battery transients. Typically, the number of RC pairs 
% ranges from 1 through 5.

% Step 1: Load and Preprocess Data
close all
%
% Create a Pulse Sequence object, which represents a pulse sequence
% experiment
psObj = Battery.PulseSequence;
% disp(psObj)

% Load in data
% Specify the file name
%FileName = 'Synthetic_LiPo_PulseDischarge.mat';
%fn = 'capture_20230815_1.txt';
fn = '20230919_batteryTest.txt';
temp = importdata(fn,',');
t = temp(:,1); voltage = temp(:,2); current = temp(:,3); 
t = datevec(t); time = t(:,4) + t(:,5)/60 + t(:,6)/3600;
time = (time - time(1))*3600;

k = [1;diff(time)];
kk = (k==0);
idx = 1:length(time);

plot(idx,time,'k',idx(kk),time(kk),'r.')
% time(kk)=[];
% voltage(kk)=[];
% current(kk)=[];
% clear k kk

tc_idx = 95741:95748;
M = (time(tc_idx(end)+1)-time(tc_idx(1)-1))/(length(tc_idx)+2);
temp2 = M*ones(length(tc_idx),1);
temp2 = cumsum(temp2);
temp2 = temp2 + time(tc_idx(1)-1); 
time(tc_idx) = temp2;

% hold on; plot(tc_idx,time(tc_idx),'b--'); hold off
clear t temp

%len = length(voltage);
%vn = 0; cn = 0;
% voltage_new = voltage + vn^2*randn(len,1);
% current_new = current + cn^2*randn(len,1);
%voltage_new = voltage + 0.05*sign(randn(len,1)) + vn^2*randn(len,1);
%current_new = current + 0.05*sign(randn(len,1))+ cn^2*randn(len,1);

% -- Smooth data
if 1
switch 4
    case 1
        ord = 1; wlen = 11;
        voltage_new = sgolayfilt(voltage,ord,wlen);
        current_new = sgolayfilt(current,ord,wlen);
        ord = 1; wlen = 9;
        voltage_new = sgolayfilt(voltage_new,ord,wlen);
        current_new = sgolayfilt(current_new,ord,wlen);
        ord = 1; wlen = 7;
        voltage_new = sgolayfilt(voltage_new,ord,wlen);
        current_new = sgolayfilt(current_new,ord,wlen);
    case 2
        wlen = 20;
        voltage_new = movmean(voltage,wlen);
        current_new = movmean(current,wlen);
%         wlen = 8;
%         voltage_new = movmean(voltage_new,wlen);
%         current_new = movmean(current_new,wlen);
%         wlen = 12;
%         voltage_new = movmean(voltage_new,wlen);
%         current_new = movmean(current_new,wlen);
    case 3
        ord = 13; wlen = 101;
        voltage_new = sgolayfilt(voltage,ord,wlen);
        current_new = sgolayfilt(current,ord,wlen);
    case 4
        wlen = 11; smoothtype = 'movmean';
        voltage_new = smoothdata(voltage,smoothtype,wlen);
        current_new = smoothdata(current,smoothtype,wlen);
    otherwise
        voltage_new = voltage;
        current_new = current;
end
figure(1)
plot(time,voltage,'k',time,voltage_new,'r--',time,current,'k',time,current_new,'b--')
voltage = voltage_new;
current = current_new;
end

% r = 3;
% time = decimate(time,r,"fir");
% voltage = decimate(voltage,r,'fir');
% current = decimate(current,r,'fir');
% voltage(time<0) = [];
% current(time<0) = [];
% time(time<0) = [];

% Read the raw data
%[time,voltage,current] = Battery.loadDataFromMatFile(FileName);
% load test2_pulseDischarge_data.mat
% kin = 73865; time = time(1:kin); voltage = voltage(1:kin); current = current(1:kin);

% Add the data to the PulseSequence
psObj.addData(time,voltage,current);

% Review the data
psObj.plot();

% Break up the data into Battery.Pulse objects
%
% Create the Pulse objects within the PulseSequence. This creates a
% parameter object that contains look-up tables of the correct size, given
% the number of pulses
psObj.createPulses(...
    'CurrentOnThreshold',12,... %minimum current magnitude to identify pulse events
    'NumRCBranches',3,... %how many RC pairs in the model
    'RCBranchesUse2TimeConstants',false,... %do RC pairs have different time constant for discharge and rest?
    'PreBufferSamples',10,... %how many samples to include before the current pulse starts
    'PostBufferSamples',15); %how many samples to include after the next pulse starts

% Correct for pulses that are not normal 
% psObj.NumPulses = 10;
% psObj.idxEdge = psObj.idxEdge(1:40,:);
% psObj.idxLoad = psObj.idxLoad(1:10,:);
% psObj.idxRelax = psObj.idxRelax(1:10,:);
% psObj.Pulse = psObj.Pulse(1:10);
if 0
    idxRemove = [3,4,5,6,7,9,10,11,13,15,16,17,21];
    psObj.removePulses(idxRemove)
end

% Specify the Simulink model that matches the number of RC branches and
% time constants:
psObj.ModelName = 'BatteryEstim3RC_PTBS';

% Plot the identified pulses in the data
psObj.plotIdentifiedPulses();

%Note: if for some reason we want to exclude some of the pulses, you can do
%something like this:
% psObj.removePulses(indexToRemove);
% psObj.plotIdentifiedPulses();


%% Step 2: Determine the Number of RC Pairs
%
% This step helps us decide how many RC pairs should be used in the model.
% More RC pairs add complexity and might over-fit the data. Too few 
% increases the fit error. Note: If you decide to change the number of
% pairs, you need to rerun createPulses above and change NumRCBranches to
% the new value.

% Pick a pulse near the beginning, middle, and end. Note: you could run all
% of them if you want.
PulsesToTest = [1 floor(psObj.NumPulses/2), psObj.NumPulses-1];
% PulsesToTest = [1 3 6 9];

% Perform the comparison
psObj.Pulse(PulsesToTest).compareRelaxationTau();

%% Step 3: Estimate Parameters
% Set settings
% Pull out the parameters. Only update parameters once because it changes
% the history each time they update.
Params = psObj.Parameters;

% Set Em constraints and initial guesses (or don't and try the defaults)
Params.Em(:) = 4.4*12;
Params.EmMin(:) = 3.4*12;
Params.EmMax(:) = 4.4*12;

% Set R0 constraints and initial guesses (or don't and try the defaults)
Params.R0(:) = 0.001;
Params.R0Min(:) = 1e-5;
Params.R0Max(:) = 1;

% Set Tx constraints and initial guesses. It is important to default each
% time constant Tx to different values, otherwise the optimizer may not
% pull them apart.
Params.Tx(1,:,:) = 5;
Params.Tx(2,:,:) = 50;
Params.Tx(3,:,:) = 200;
% Params.Tx(4,:,:) = 500;

Params.TxMin(1,:,:) = 1;
Params.TxMax(1,:,:) = 50;

Params.TxMin(2,:,:) = 20;
Params.TxMax(2,:,:) = 1000;

Params.TxMin(3,:,:) = 100;
Params.TxMax(3,:,:) = 3500; %don't set this bigger than the relaxation time available


% Params.TxMin(4,:,:) = 200;
% Params.TxMax(4,:,:) = 10000;
% Set Rx constraints and initial guesses (or don't and try the defaults)
Params.Rx(:) = 0.01;
Params.RxMin(:) = 1e-5;
Params.RxMax(:) = 0.5;

% Update parameters
psObj.Parameters = Params;

% Estimate initial Em and R0 values
%
% This step inspects the voltage immediately before and after the current
% is applied and removed at the start and end of each pulse. It uses that
% for a raw calculation estimating what the open-circuit voltage (Em) and
% the series resistance R0 should be.

psObj.estimateInitialEmR0(...
    'SetEmConstraints',false,... %Update EmMin or EmMax values based on what we learn here
    'EstimateEm',true,... %Keep this on to perform Em estimates
    'EstimateR0',true); %Keep this on to perform R0 estimates

% Plot results
psObj.plotLatestParameters();

% Get initial Tx (Tau) values
%
% This step performs curve fitting on the pulse relaxation to estimate the
% RC time constant at each SOC.

psObj.estimateInitialTau(...
    'UpdateEndingEm',false,... %Keep this on to update Em estimates at the end of relaxations, based on the curve fit
    'ShowPlots',true,... %Set this true if you want to see plots while this runs
    'ReusePlotFigure',true,... %Set this true to overwrite the plots in the same figure
    'UseLoadData',false,... %Set this true if you want to estimate Time constants from the load part of the pulse, instead of relaxation
    'PlotDelay',0.5); %Set this to add delay so you can see the plots 

% Plot results
psObj.plotLatestParameters(); %See what the parameters look like so far
psObj.plotSimulationResults(); %See what the result looks like so far

%
% Get initial Em and Rx values using a linear system approach - pulse by 
% pulse
%
% This step takes the data for each pulse and treats it as a linear system
% It attempts to fit the Rx values for each RC branch. Optionally, you can
% allow it to adjust the Em and R0 values, and if these are adjusted, you
% also have the option whether to retain the optimized values of these or
% to discard them.

psObj.estimateInitialEmRx(...
    'IgnoreRelaxation',false,... %Set this true if you want to ignore the relaxation periods during this step
    'ShowPlots',true,...  %Set this true if you want to see plots while this runs
    'ShowBeforePlots',true,... %Set this true if you want to see the 'before' value on the plots
    'PlotDelay',0.5,... %Set this to add delay so you can see the plots 
    'EstimateEm',true,... %Set this true to allow the optimizer to change Em further in this step
    'RetainEm',true,... %Set this true keep any changes made to Em in this step
    'EstimateR0',true,... %Set this true to allow the optimizer to change R0 further in this step
    'RetainR0',true); %Set this true keep any changes made to R0 in this step

% Plot results
psObj.plotLatestParameters(); %See what the parameters look like so far
psObj.plotSimulationResults(); %See what the result looks like so far


% Perform SDO Estimation
SDOOptimizeOptions = sdo.OptimizeOptions(...
    'OptimizedModel',psObj.ModelName,...
    'Method','lsqnonlin',...
    'UseParallel','always');

psObj.estimateParameters(...
    'CarryParamToNextPulse',true,... %Set this true to use the final parameter values from the prior pulse and SOC as initial values for the next pulse and SOC
    'SDOOptimizeOptions',SDOOptimizeOptions,... %Specify the SDO options object
    'ShowPlots',true,... %Set this true if you want to see plots while this runs
    'EstimateEm',true,... %Set this true to allow the optimizer to change Em further in this step
    'RetainEm',true,... %Set this true keep any changes made to Em in this step
    'EstimateR0',true,... %Set this true to allow the optimizer to change R0 further in this step
    'RetainR0',true); %Set this true keep any changes made to R0 in this step

% Plot results
psObj.plotLatestParameters(); %See what the parameters look like so far
psObj.plotSimulationResults(); %See what the result looks like so far

%% Step 4: Set Equivalent Circuit Battery Block Parameters
%
% The experiment was run at ambient temperature (303°K) only.  Repeat 
% the tables across the operating temperature range.  If the discharge 
% experiment was run at 2 different constant temperatures, then include 
% these in the tables below. 

load 20230919_data.mat

EmPrime = repmat(Em,2,1)';
R0Prime = repmat(R0,2,1)';
SOC_LUTPrime = SOC_LUT;
TempPrime = [303 315.15];
CapacityAhPrime = [CapacityAh CapacityAh];

R1Prime = repmat(Rx(1,:),2,1)';
C1Prime = repmat(Tx(1,:)./Rx(1,:),2,1)';
R2Prime = repmat(Rx(2,:),2,1)';
C2Prime = repmat(Tx(2,:)./Rx(2,:),2,1)';
R3Prime = repmat(Rx(3,:),2,1)';
C3Prime = repmat(Tx(3,:)./Rx(3,:),2,1)';
% R4Prime = repmat(Rx(4,:),2,1)';
% C4Prime = repmat(Tx(4,:)./Rx(4,:),2,1)';

open_system('BatteryEstim3RC_PTBS_EQ');

