clc
clear

if 1

%     load('2022-08-04 10-58-30.tlog.mat')
%     ind = [2276,2331; 2331,2463; 2463,4628; 4628,4767];
%     ind2 = [6999,7053; 7053,7196; 7196,9976; 9976,10113];
%     ind = [ind;ind2];

    load('2022-08-10 09-40-57.tlog.mat');
    ind = [3734,3794; 3794,3874; 3874,17878; 17878,17948];
    
    % get data from log file
    I = current_battery_mavlink_battery_status_t;
    ti = I(:,1); I = I(:,2)/100;
    % C = current_consumed_mavlink_battery_status_t(:,2);
    tv = voltage_battery_mavlink_sys_status_t(:,1);
    V = voltage_battery_mavlink_sys_status_t(:,2)/1000;
    Iv = interp1(ti,I,tv);
    P = V.*Iv;
    tv = datevec(tv);
    tv = tv(:,4)*3600 + tv(:,5)*60 + tv(:,6);
    dt = [0;diff(tv)];
    Pint = [cumsum(P(1:end,1).*dt,'omitnan')/3600/45.6]; %(capacity over time)

    % save what you need and clear rest
    save simpleScriptVariables ti tv I V Iv P dt Pint ind
    clear
end
load simpleScriptVariables

%% Plot basic results
figure(1)
subplot 131
yyaxis left
plot(tv,V,'k')
title('Voltage/Current vs time')
xlabel('time, s')
ylabel('V','Rotation',0)
yyaxis right
plot(tv,Iv,'r')
ylabel('A','Rotation',0)
legend('Voltage','Current')

ax1 = gca;

subplot 132
plot(tv,P)
title('Power vs time')
xlabel('time, s')
ylabel('W','Rotation',0)

ax2 = gca;

subplot 133
plot(tv,Pint)
title('Capacity vs time')
xlabel('time, s')
ylabel('Ah','Rotation',0)

ax3 = gca;

linkaxes([ax1,ax2,ax3],'x')

%% Change of capacity (rate of use)
figure(2)
dCapdt = 3600*[0;diff(Pint)]./dt;
% i = 2;
% j = 2;
%
% k = 1;
% subplot(i,j,k)
plot(tv,dCapdt)

%% Change of capacity (rate of use) - 2022 08 04
if 0
    figure(3)
    dCapdt = 3600*[0;diff(Pint)]./dt;
    i = 2;
    j = 4;

    k = 1;
    subplot(i,j,k)
    histogram(dCapdt(ind(k,1):ind(k,2)),40,'Normalization','pdf')
    title('Takeoff')
    xlabel('Ah/h')
    ylabel('Probability')

    k = 2;
    subplot(i,j,k)
    histogram(dCapdt(ind(k,1):ind(k,2)),40,'Normalization','pdf')
    title('Transistion to Fixed-wing')
    xlabel('Ah/h')
    ylabel('Probability')

    k = 3;
    subplot(i,j,k)
    histogram(dCapdt(ind(k,1):ind(k,2)),40,'Normalization','pdf')
    title('Cruising')
    xlabel('Ah/h')
    ylabel('Probability')
    % xlim([2,20])

    k = 4;
    subplot(i,j,k)
    histogram(dCapdt(ind(k,1):ind(k,2)),40,'Normalization','pdf')
    title('Landing')
    xlabel('Ah/h')
    ylabel('Probability')
    % xlim([30,70])
    % ylim([0,1e-3])

    k = 5;
    subplot(i,j,k)
    histogram(dCapdt(ind(k,1):ind(k,2)),40,'Normalization','pdf')
    title('Takeoff')
    xlabel('Ah/h')
    ylabel('Probability')

    k = 6;
    subplot(i,j,k)
    histogram(dCapdt(ind(k,1):ind(k,2)),40,'Normalization','pdf')
    title('Transistion to Fixed-wing')
    xlabel('Ah/h')
    ylabel('Probability')

    k = 7;
    subplot(i,j,k)
    histogram(dCapdt(ind(k,1):ind(k,2)),40,'Normalization','pdf')
    title('Cruising')
    xlabel('Ah/h')
    ylabel('Probability')
    % xlim([2,20])

    k = 8;
    subplot(i,j,k)
    histogram(dCapdt(ind(k,1):ind(k,2)),40,'Normalization','pdf')
    title('Landing')
    xlabel('Ah/h')
    ylabel('Probability')

end

tphase = tv(ind);
tphase = tphase(:,2)-tphase(:,1);
%% Change of capacity (rate of use) - 2022 08 10
if 1
    figure(4)
    dCapdt = 3600*[0;diff(Pint)]./dt;
    i = 2;
    j = 2;

    k = 1;
    subplot(i,j,k)
    histogram(dCapdt(ind(k,1):ind(k,2)),40,'Normalization','pdf')
    title('Takeoff')
    xlabel('Ah/h')
    ylabel('Probability')

    k = 2;
    subplot(i,j,k)
    histogram(dCapdt(ind(k,1):ind(k,2)),40,'Normalization','pdf')
    title('Transistion to Fixed-wing')
    xlabel('Ah/h')
    ylabel('Probability')

    k = 3;
    subplot(i,j,k)
    histogram(dCapdt(ind(k,1):ind(k,2)),40,'Normalization','pdf')
    title('Cruising')
    xlabel('Ah/h')
    ylabel('Probability')
    % xlim([2,20])

    k = 4;
    subplot(i,j,k)
    histogram(dCapdt(ind(k,1):ind(k,2)),40,'Normalization','pdf')
    title('Landing')
    xlabel('Ah/h')
    ylabel('Probability')
    % xlim([30,70])
    % ylim([0,1e-3])
end

tphase = tv(ind);
tphase = tphase(:,2)-tphase(:,1);
%% testing
% winLen = 200;
% for i=winLen+1:length(V)
%     % take best fit line
%     p = polyfit(tv(i-winLen:i),V(i-winLen:i),1);
%     Vm = polyval(p,tv(i-winLen:winLen));
%
%
%     % record slope
%     m(i,1) = p(1);
%     shg
%
% end
% figure(2)
% plot(tv,m,tv,movstd(V,winLen),tv,movmad(V,winLen))
%
% ax4 = gca;
%
% linkaxes([ax1,ax2,ax3,ax4],'x')
%% Position variance versus local position

% x1 = pos_horiz_variance_mavlink_ekf_status_report_t(:,2);
% x2 = pos_vert_variance_mavlink_ekf_status_report_t(:,2);
% t1 = pos_horiz_variance_mavlink_ekf_status_report_t(:,1);
% t2 = pos_vert_variance_mavlink_ekf_status_report_t(:,1);
%
% ax1 = subplot(1,2,1);
% plot(t1,x1,t2,x2)
% legend('horz var','vert var')
%
% subplot 122
% plot(y_mavlink_local_position_ned_t(:,1),y_mavlink_local_position_ned_t(:,2))
% hold on
% plot(x_mavlink_local_position_ned_t(:,1),x_mavlink_local_position_ned_t(:,2))
% hold off
% ax2 = gca;
% legend('y pos','x pos')
%
% linkaxes([ax1,ax2],'x' )
