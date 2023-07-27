clear
clc

% temp = readtable('20230718_pdc_1_corrected.txt');
% temp = readtable('20230718_sim_1_corrected.txt');
temp = readtable('capturebatjulypulse.txt');
data = table2array(temp(:,1:67));
data = data(1:end-1,:);

T = data(:,4:67);
T = mean(T')';

% temp2 = table2array(temp(:,1:2));
% temp3_ = table2cell(temp(:,3));
% temp5 = table2array(temp(:,4:66));
% temp6 = table2array(temp(:,67));
% temp7 = table2array(temp(:,68:end));

% for i=1:size(temp3_,1)
%     s = char(temp3_(i));
%     k = find(s=='.');
%     temp3(i) = str2double(s(1:(k(1)+2)));
%     temp4(i) = str2double(s(k(1)+3:end));
% end

% data = [temp2,temp3',temp4'];%,temp5];%,temp7];


%%
figure(1)
clf

ax1 = nexttile;
plot(data(:,1),data(:,2))
xlabel('matlab time')
ylabel('voltage (V)')
title('Voltage vs time')
%legend('Voltage')
ylim([42,55])

ax2 = nexttile;
plot(data(:,1),data(:,3))
xlabel('matlab time')
ylabel('current (A)')
title('Current vs time')
%legend('Voltage')

ax3 = nexttile;
z0 = reshape(data(1,4:67),[8,8]);
surf(z0)
shading interp
hold on;
zf = reshape(data(end,4:67),[8,8]);
surf(zf)
shading interp
hold off
zlim([0,50])
caxis([20,50])
xlabel('pixels')
ylabel('pixels')
title('Temperature vs time')
legend('Initial temperature','Final temperature')

ax4 = nexttile;
plot(data(:,1),T)
xlabel('matlab time')
ylabel('Temperature (C)')
title('Temperature vs time')
%legend('Voltage')

linkaxes([ax1,ax2,ax4],'x')