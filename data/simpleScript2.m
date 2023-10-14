clc
clear
% load '2022-08-10 09-40-57.tlog.mat'
load '2022-07-29 10-38-55.tlog.mat'

t = lng_mavlink_camera_feedback_t(:,1);
x = lat_mavlink_camera_feedback_t(:,2);
y = lng_mavlink_camera_feedback_t(:,2);
alt = alt_msl_mavlink_camera_feedback_t(:,2);

% figure(1)
% i = 11; j = 1;
% yyaxis left
% plot(t(i:end-j),x(i:end-j),'o-')
% yyaxis right
% plot(t(i:end-j),y(i:end-j),'o-')


xnew = x(i:end-j);
ynew = y(i:end-j);
tnew = t(i:end-j);
altnew = alt(i:end-j);

dx = diff(xnew);
dy = diff(ynew);
dt = diff(tnew);

normdxdy = sqrt(dx.^2 + dy.^2);

% thres = 2900;
% xnew(normdxdy<thres) = [];
% ynew(normdxdy<thres) = [];

thres = 1.52545*1e-5;
xnew(dt<thres) = [];
ynew(dt<thres) = [];
altnew(dt<thres) = [];


figure(2)
%plot(tnew(2:end),normdxdy)
scatter(xnew,ynew);

writematrix([xnew,ynew,altnew],'20220810_cam_updated.csv')