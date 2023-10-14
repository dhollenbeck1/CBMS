clear
clc

% -- Algorithm
% 1. load in experimental data, lookup table data, and fitted parameter data
% 2. Clean up database of unsued parameters
% 3. Create the simin data using the voltage and current.
% 4. Create output data and expand params for plotting


addpath data\

load 20230919_data_fit.mat

% clear ax0 ax1 ax2 ax3 ax4 ax5 ax6 current_new current idx logsout
%%
% load data\'2022-07-29 10-38-55_f1.mat'
% load data\'2022-08-04 10-58-30_f1.mat'
% load data\'2022-08-04 10-58-30_f2.mat'
% load data\'2022-08-08 10-44-35_f1.mat'
% load data\'2022-08-10 09-40-57_f1.mat'
% load data\'2022-09-26 11-49-08_f1.mat'
% load data\'2022-10-21 13-31-20_f1.mat'

fn = [{'2022-07-29 10-38-55_f1.mat'},
    {'2022-08-04 10-58-30_f1.mat'}
    {'2022-08-04 10-58-30_f2.mat'}
    {'2022-08-08 10-44-35_f1.mat'}
    {'2022-08-10 09-40-57_f1.mat'}
    {'2022-09-26 11-49-08_f1.mat'}
    {'2022-10-21 13-31-20_f1.mat'}];

CAh = 24*ones(1,length(fn));
DefCo = ones(1,length(fn));
p = [CAh;DefCo];
J = @(x,y) sum((x-y).^2,'omitnan');
Jold = ones(1,length(fn));
Jnew = ones(1,length(fn));
Jvec = [];
dJ = ones(1,length(fn));
thres = 1e-4;
looplim = 100;
lr = [0.006; -0.0003];

for j=3%length(fn)
    % -- load data for jth set
    Fn = char(fn(j));
    load(Fn)

    % -- clean up current nans
    i(isnan(i))=0;

    % -- find important params
    Em0 = v(2);%;max(v);
    %InitialCapVoltage = 0*ones(1,3);
    SOC_0 = interp1(Em,SOC_LUT,Em0);
    CapacityAh =  24.1;%p(1,j);
    CapacityAhPrime = CapacityAh*ones(1,2);
%     InitialChargeDeficitAh = p(2,j)*(1-SOC_0)*CapacityAh;%5.75;
    InitialChargeDeficitAh = 0.91*(1-SOC_0)*CapacityAh;%5.75;
    t = t - t(1);

    % EmPrime = Em;
    % SOC_LUTPrime = SOC_LUT;
    % TempPrime = [303,315];
    % CapacityAhPrime = CapacityAh;
    R0prime = f(SOC_LUT,thet0);
    R1prime = f(SOC_LUT,thet1);
    R2prime = f(SOC_LUT,thet2);
    R3prime = f(SOC_LUT,thet3);
    C1prime = f(SOC_LUT,thet4);
    C2prime = f(SOC_LUT,thet5);
    C3prime = f(SOC_LUT,thet6);

    % -- optimization loop
    loopcount = 0;
    Jvec = [];
    while true
        loopcount = loopcount + 1;
        clc
        disp(num2str(loopcount))
        sim("BatteryEstim3RC_PTBS_EQ_v1.slx")

        t1 = yout.getElement(1).Values.Time;
        v1 = yout.getElement(1).Values.Data;
        t2 = yout.getElement(4).Values.Time;
        v2 = yout.getElement(4).Values.Data;
        X1 = yout.getElement(2).Values.Data;
        X2 = yout.getElement(3).Values.Data;

        if loopcount >= looplim
            disp('loop limit reached')
            break;
        else
            v1c = interp1(t1,v1,t);
            
            wlen = 1;
            Jnew(j) = J(movmean(v,wlen),movmean(v1c,wlen));
            dJ(j) = Jnew(j) - Jold(j);
            Jold(j) = Jnew(j);
            Jvec(loopcount) = Jnew(j);
            pdir = mod(loopcount,2)+1;%ceil(rand()*2);
            p(pdir,j) = p(pdir,j) + lr(pdir)*dJ(j);
            p(:,j) = p(:,j) + lr(:).*dJ(j);
            p(:,j) = p(:,j).*(1-sign(dJ(j))*0.01^2*randn(2,1));
            CapacityAh = p(1,j); %21.3;
            CapacityAhPrime = CapacityAh*ones(1,2);
            InitialChargeDeficitAh = p(2,j)*(1-SOC_0)*CapacityAh;
            disp(['p dir is ',num2str(pdir)])
            disp(['p',num2str(pdir),' is ',num2str(p(pdir,j))])
        end

        if abs(dJ(j)) <= thres
            disp(['Threshold reached ',num2str(dJ(j)), ', count is ',num2str(loopcount)])
            break;
        end

        figure(5)
        clf
        axn = nexttile;
        yyaxis left
        plot(t,movmean(v,wlen),'k'); hold on
        plot(t1,movmean(v1,wlen),'g--')
        plot(t2,v2,'b--')
        hold off
        ylim([39,51]);
        % xlim([1800,t(end)]);
        legend('Measured','Parameterized','Lookup Table')
        xlabel('time, s'); ylabel('Voltage, V')
        title('Real flight voltage vs model performance')

        yyaxis right
        plot(t,i,'r')
        ylabel('Current, A')
        set(gca,'FontSize',14)
        ax = gca;
        ax.YAxis(2).Color = 'r';
        % open_system('BatteryEstim3RC_PTBS_EQ_v1.slx');
        
        axm = nexttile;
        plot(1:loopcount,Jvec)
        
        shg
%         pause(0.5)
    end


    % end



end