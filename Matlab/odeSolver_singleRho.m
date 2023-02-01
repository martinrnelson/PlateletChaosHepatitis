%% Run ODE solver to simulate the system (1) from the manuscript
% (This code reproduces Figure 3. Pass the output to plotAttractor.m to
% reproduce Figures 5 and 6)

function sol=odeSolver_singleRho()

    % Plots on/off
    plotsOn=1;

    % Parameters
    p=setBaselineParams();
    p.rho=4;

    % Initial conditions;
    n0 = 0;
    a0 = 0;
    m10 = 0;
    m20 = 0;
    c0 = 0.04;
    g0 = 0;
    h0 = 1;
    ha0 = 0;
    e0 = 0;
    s0 = 1;
    sa0 = 0;

    ICs=[n0,a0,m10,m20,c0,g0,h0,ha0,e0,s0,sa0];

    % Integration range
    tspan=[0 200];

    % Call ODE solver
    sol = ode45(@(t,y)odeSystem_singleRho(t,y,p),tspan,ICs);
    t=sol.x;
    y=sol.y;

    % Plot timecourses (Figure 3)
    if(plotsOn)
        close all
        figure(1);

        subplot(2,3,1);
        plot(t,y(1,:),'r');
        hold on
        plot(t,y(2,:),'r--');
        xlabel('t')
        ylabel('Neutrophils');

        subplot(2,3,2);
        plot(t,y(3,:),'g');
        hold on
        plot(t,y(4,:),'g--');
        xlabel('t')
        ylabel('Macrophages');

        subplot(2,3,3);
        plot(t,y(5,:),'r');
        xlabel('t')
        ylabel('Pro mediators');

        subplot(2,3,4);
        plot(t,y(6,:),'b');
        xlabel('t')
        ylabel('Anti mediators');

        subplot(2,3,5);
        plot(t,y(7,:),'Color',[0.8500 0.3250 0.0980],'LineStyle','-');
        hold on
        plot(t,y(8,:),'Color',[0.8500 0.3250 0.0980],'LineStyle','--');
        plot(t,y(9,:),'k');
        xlabel('t')
        ylabel('Hepatocytes/ECM');
        ylim([-0.1 1.1]);

        subplot(2,3,6);
        plot(t,y(10,:),'b');
        hold on
        plot(t,y(11,:),'b--');
        xlabel('t')
        ylabel('Stellate cells');
        ylim([-0.1 1.1]);
    end

end
