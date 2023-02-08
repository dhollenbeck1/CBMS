clc
clear
clf

addpath data\ data_cut\ data_raw\ auxillary\

% -- Import Data
fn = 'data\2022-08-10 09-40-57.tlog.mat';
fn2 = 'data_cut\data_20220810.mat';
T = 5000;
dt = 1e-0;
n = round(T/dt);
t = linspace(0,T,n)';
load(fn)

ve = voltage_battery_mavlink_sys_status_t;
ie = current_battery_mavlink_sys_status_t;
wp = seq_mavlink_mission_current_t(:,2);
wp_t = seq_mavlink_mission_current_t(:,1);
te = ve(:,1);
ve = ve(:,2)/1e3;
ie = ie(:,2)/1e2;
k_wp = find(wp==212,1,'first');
wp(k_wp:end) = 50;
wp(1:k_wp-1) = 0;
wp = interp1(wp_t,wp,te);
wp(isnan(wp)) = 0;
ind = 3733;
te = te(ind:end);
ve = ve(ind:end);
ie = ie(ind:end);
ie(ie<1e-2) = 1e-2;
wp = wp(ind:end);
re = ve./ie;

temp = datevec(te);
te = temp(:,4)*3600 + temp(:,5)*60 + temp(:,6);
te = te - te(1);

data.te = te;
data.ve = ve;
data.ie = ie;
data.re = re;
data.wp = wp;

save(fn2,'data')

plot(te,ve,'k',te,ie,'r',te,re,'g')
