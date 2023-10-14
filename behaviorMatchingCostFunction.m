function J = behaviorMatchingCostFunction(x)
x
% -- load in global params
global mpg
mp = mpg;
% load(char(mp.fn(mp.j)));

% -- find important params
mpg.Em0 = mpg.v(5);%;max(v);
mpg.SOC_0 = interp1(mpg.Em,mpg.SOC_LUT,mpg.Em0);
mpg.CapacityAh =  x(1);%24.1;%p(1,j);
mpg.CapacityAhPrime = x(1)*ones(1,2);
mpg.InitialChargeDeficitAh = x(2)*(1-mpg.SOC_0)*mpg.CapacityAh;%5.75;

% -- run simulation
sim(mpg.simulinkModel);

% -- extract the data
mpg.t1 = yout.getElement(1).Values.Time;
mpg.v1 = yout.getElement(1).Values.Data;
mpg.t2 = yout.getElement(4).Values.Time;
mpg.v2 = yout.getElement(4).Values.Data;
mpg.X1 = yout.getElement(2).Values.Data;
mpg.X2 = yout.getElement(3).Values.Data;

% -- compute cost function
mpg.v1c = interp1(mpg.t1,mpg.v1,mpg.t);
switch 1
    case 1
        J = mpg.J(mpg.v,mpg.v1c);
    case 2
        wlen = 10;
        J = mpg.J(movmean(mpg.v,wlen),movmean(mpg.v1c,wlen));
end

if mpg.PLOT
    figure(5)
    clf
    axn = nexttile;
    yyaxis left
    %plot(mp.t,movmean(v,wlen),'k'); hold on
    %plot(t1,movmean(v1,wlen),'g--')
    plot(mpg.t,mpg.v,'k'); hold on
    plot(mpg.t1,mpg.v1,'g--')
    plot(mpg.t2,mpg.v2,'b--')
    hold off
    ylim([39,51]);% xlim([1800,t(end)]);
    legend('Measured','Parameterized','Lookup Table')
    xlabel('time, s'); ylabel('Voltage, V')
    title('Real flight voltage vs model performance')
    yyaxis right
    plot(mpg.t,mpg.i,'r')
    ylabel('Current, A')
    set(gca,'FontSize',14)
    ax = gca;
    ax.YAxis(2).Color = 'r';

    % axm = nexttile;
    % plot(1:loopcount,Jvec)

    shg
end

end