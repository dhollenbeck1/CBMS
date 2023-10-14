clc
clear

fn = '2022-07-29 10-38-55.tlog.mat';    % -- test flight of demo
% fn = '2022-08-04 10-58-30.tlog.mat';   aidx = [2198,4978;6773,10185]; % -- aircraft tuning flight
% fn = '2022-08-08 10-44-35.tlog.mat'; aidx = [6776,19062];   % -- test flight endurance
% fn = '2022-08-10 09-40-57.tlog.mat';    % -- test flight of demo (accident)
% fn = '2022-09-26 11-49-08.tlog.mat'; aidx = [2740,6020];   % -- survey flight
% fn = '2022-10-21 13-31-20.tlog.mat'; aidx = [2060,11897];   % -- demo

load(fn) 
% -- Settings
lw = 1.5;
fntsze = 14;
fnum = 1;


% -- Extract the position, voltage and current information
t = x_mavlink_local_position_ned_t(2:end,1);
x = x_mavlink_local_position_ned_t(2:end,2);
y = y_mavlink_local_position_ned_t(2:end,2);
alt = -z_mavlink_local_position_ned_t(2:end,2);
i = current_battery_mavlink_sys_status_t(2:end,2);
i_t = current_battery_mavlink_sys_status_t(2:end,1);
i = interp1(i_t,i,t)/100;
v = voltage_battery_mavlink_sys_status_t(2:end,2);
v = interp1(i_t,v,t)/1000;

% -- convert to seconds of the day
t = datevec(t);
t = t(:,6) + t(:,5)*60 + t(:,4)*3600;

% -- shift so time starts at t0 = 0;
t = t-t(1); 

% -- convert to hours
dt = [0;diff(t)]/3600; 

% -- integrate the current to estimate the energy usage
idt = i.*dt;
idt(isnan(idt)) = 0;    % -- clean up any NaN's
I = cumsum(idt);        % -- compute the cumlative sum

% -- get mode change heartbeat
mode_t = custom_mode_mavlink_heartbeat_t(:,1);
mode_t = datevec(mode_t);
mode_t = mode_t(:,6) + mode_t(:,5)*60 + mode_t(:,4)*3600;
mode_t = mode_t - mode_t(1);
mode = custom_mode_mavlink_heartbeat_t(:,2);

% -- Get the flight controller mode
basemode_t = base_mode_mavlink_heartbeat_t(:,1);
basemode_t = datevec(basemode_t);
basemode_t = basemode_t(:,6) + basemode_t(:,5)*60 + basemode_t(:,4)*3600;
basemode_t = basemode_t - basemode_t(1);
basemode = base_mode_mavlink_heartbeat_t(:,2);

% setmode_t = base_mode_mavlink_set_mode_t(:,1);
% setmode_t = datevec(setmode_t);
% setmode_t = setmode_t(:,6) + setmode_t(:,5)*60 + setmode_t(:,4)*3600;
% setmode_t = setmode_t - setmode_t(1);
% setmode = base_mode_mavlink_set_mode_t(:,2);

% -- if adjusted idx is avail. trim dataset
if exist('aidx')
t = t(aidx(fnum,1):aidx(fnum,2));
x = x(aidx(fnum,1):aidx(fnum,2));
y = y(aidx(fnum,1):aidx(fnum,2));
alt = alt(aidx(fnum,1):aidx(fnum,2));
i = i(aidx(fnum,1):aidx(fnum,2));
v = v(aidx(fnum,1):aidx(fnum,2));
I = I(aidx(fnum,1):aidx(fnum,2));

ai = find(mode_t>=t(1));
af = find(mode_t<=t(end));
mode_t = mode_t(ai(1):af(end));
mode = mode(ai(1):af(end));
end

% -- plot results to verify the regions for 
figure(1)
clf
ax1 = nexttile;
plot(y,x,'k','LineWidth',lw)
title('Flight Path')
xlabel('m'); ylabel('m')

ax2 = nexttile;
plot(t,alt,'k','LineWidth',lw)
title('Flight Altitude')
xlabel('time, s'); ylabel('m')
ylim([0,150])
yyaxis right
plot(mode_t,mode,'r')
% hold on; plot(basemode_t,basemode/10,'g--'); hold off
ax2.YAxis(2).Color = 'r';

ax3 = nexttile;
plot(t, i,'k','LineWidth',lw) 
title('Battery System Current Draw and Voltage')
xlabel('time, s'); ylabel('A')
ylim([0,100])
yyaxis right
plot(t, v, 'LineWidth',lw,'Color','b')
ylabel('V')
ax3.YAxis(2).Color = 'b';
ylim([38,54])

ax4 = nexttile;
plot(t,I,'k','LineWidth',lw)
title('Battery Energy Consumption')
xlabel('time, s'); ylabel('Ah')
ylim([0,25])

linkaxes([ax2,ax3,ax4],'x')

set(ax1,'FontSize',fntsze)
set(ax2,'FontSize',fntsze)
set(ax3,'FontSize',fntsze)
set(ax4,'FontSize',fntsze)

saveas(gcf,[fn(1:end-9),'_f',num2str(fnum)],'png')
saveas(gcf,[fn(1:end-9),'_f',num2str(fnum)],'svg')
%% figure 2

k0 = find(mode == 21);      % -- find the landing points.
k1 = find(t<mode_t(k0(end))); 
k2 = find(t>=mode_t(k0(1)));
k3 = find(alt >= 90);       % -- find points at cruise alt

switch 1                    % -- number of flights
    case 1
        % -- take first point during mode 21
        idx(1,1) = k2(1); 
        % -- take last point during mode 21 and add a few 
        % samples to ensure the aircraft landed
        idx(1,2) = k1(end)+4; 
        % -- take last point near 90m before descend
        idx(1,3) = k3(end);
    case 2
        idx(1,1) = k2(end);
        idx(1,2) = k1(end)+4;
        idx(2,1) = k2(end);
        idx(2,2) = k1(end)+4; 
    otherwise
end

figure(2)
% clf
% mode_t = custom_mode_mavlink_heartbeat_t(:,1);
% mode = custom_mode_mavlink_heartbeat_t(:,2);
% plot(mode_t,mode)
n = 1;
r1 = idx(n,3):idx(n,2);
r2 = idx(n,1):idx(n,2);

clf
ax1 = nexttile;
plot(y,x,'k','LineWidth',lw); hold on
plot(y(r1),x(r1),'r','LineWidth',lw)
plot(y(r2),x(r2),'g','LineWidth',lw)
hold off
title('Flight Path')
xlabel('m'); ylabel('m')

ax2 = nexttile;
plot(t,alt,'k','LineWidth',lw); hold on
plot(t(r1),alt(r1),'r','LineWidth',lw);
plot(t(r2),alt(r2),'g','LineWidth',lw);
hold off
title('Flight Altitude')
xlabel('time, s'); ylabel('m')
ylim([0,150])
% yyaxis right
% plot(mode_t,mode,'r')
% hold on; plot(basemode_t,basemode/10,'g--'); hold off
% ax2.YAxis(2).Color = 'r';

ax3 = nexttile;
plot(t, i,'k','LineWidth',lw); hold on
plot(t(r1), i(r1),'r','LineWidth',lw);
plot(t(r2), i(r2),'g','LineWidth',lw);
hold off
title('Battery System Current Draw and Voltage')
xlabel('time, s'); ylabel('A')
ylim([0,100])
yyaxis right
plot(t, v, 'LineWidth',lw,'Color','b')
ylabel('V')
ax3.YAxis(2).Color = 'b';
ylim([38,54])

ax4 = nexttile;
plot(t,I,'k','LineWidth',lw); hold on
plot(t(r1), I(r1),'r','LineWidth',lw);
plot(t(r2), I(r2),'g','LineWidth',lw);
hold off
title('Battery Energy Consumption')
xlabel('time, s'); ylabel('Ah')
ylim([0,25])

linkaxes([ax2,ax3,ax4],'x')

set(ax1,'FontSize',fntsze)
set(ax2,'FontSize',fntsze)
set(ax3,'FontSize',fntsze)
set(ax4,'FontSize',fntsze)

saveas(gcf,[fn(1:end-9),'_f',num2str(fnum),'_2'],'png')
saveas(gcf,[fn(1:end-9),'_f',num2str(fnum),'_2'],'svg')
%% Figure 3
figure(3)
buf = 150;
tbuf = 20;
clf
ax1 = nexttile;
plot(y,x,'k','LineWidth',lw); hold on
plot(y(r1),x(r1),'r','LineWidth',lw)
plot(y(r2),x(r2),'g','LineWidth',lw)
hold off
title('Flight Path')
xlabel('m'); ylabel('m')
axis([min(y(r1))-buf,max(y(r1))+buf,min(x(r1))-buf,max(x(r1))+buf])

ax2 = nexttile;
plot(t,alt,'k','LineWidth',lw); hold on
plot(t(r1),alt(r1),'r','LineWidth',lw);
plot(t(r2),alt(r2),'g','LineWidth',lw);
hold off
title('Flight Altitude')
xlabel('time, s'); ylabel('m')
ylim([0,150])
xlim([min(t(r1))-tbuf,max(t(r1))+tbuf])
% yyaxis right
% plot(mode_t,mode,'r')
% hold on; plot(basemode_t,basemode/10,'g--'); hold off
% ax2.YAxis(2).Color = 'r';

ax3 = nexttile;
plot(t, i,'k','LineWidth',lw); hold on
plot(t(r1), i(r1),'r','LineWidth',lw);
plot(t(r2), i(r2),'g','LineWidth',lw);
hold off
title('Battery System Current Draw and Voltage')
xlabel('time, s'); ylabel('A')
ylim([0,100])
xlim([min(t(r1))-tbuf,max(t(r1))+tbuf])
yyaxis right
plot(t, v, 'LineWidth',lw,'Color','b')
ylabel('V')
ax3.YAxis(2).Color = 'b';
ylim([38,54])

ax4 = nexttile;
plot(t,I,'k','LineWidth',lw); hold on
plot(t(r1), I(r1),'r','LineWidth',lw);
plot(t(r2), I(r2),'g','LineWidth',lw);
hold off
title('Battery Energy Consumption')
xlabel('time, s'); ylabel('Ah')
% ylim([0,25])
xlim([min(t(r1))-tbuf,max(t(r1))+tbuf])

linkaxes([ax2,ax3,ax4],'x')

set(ax1,'FontSize',fntsze)
set(ax2,'FontSize',fntsze)
set(ax3,'FontSize',fntsze)
set(ax4,'FontSize',fntsze)

saveas(gcf,[fn(1:end-9),'_f',num2str(fnum),'_3'],'png')
saveas(gcf,[fn(1:end-9),'_f',num2str(fnum),'_3'],'svg')
%% Save variables
save([fn(1:end-9),'_f',num2str(fnum)],'t','x','y','alt',...
    'i','v','I','mode_t','mode','r1','r2');