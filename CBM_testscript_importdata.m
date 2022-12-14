clc
clear

fn = 'data_test1';
x = importdata('test1.txt',',');
t = datevec(x(:,1));
x = x(:,3:4);
t = 3600*t(:,4) + 60*t(:,5) + t(:,6);

x(t==0,:) = [];
t(t==0) = [];

t = t-t(1);
T = t(end)/3600;

plot(t,x(:,1),'k')
yyaxis right
plot(t,x(:,2),'r')

data.x = x; 
data.t = t;
data.T = T;

save([fn,'.mat'],'data')

writematrix([t,x],[fn,'.csv'])