%% Plot 1D trajectories and 2D (c,g) slices of attractor
% Input: solution structure produced by odeSolver_singleRho.m
% (Produces Figures 5 and 6)

function plotAttractor(sol)

    t=sol.x;
    y=sol.y;

    close all;

    figure(1);
    plot(t,y(2,:),'r');
    hold on
    plot(t,y(1,:),'b');
    xlabel('t')
    ylabel('Neutrophils');
    ylim([0 0.1]);

    figure(2);
    plot(t,y(3,:),'b');
    hold on
    plot(t,y(4,:),'r');
    xlabel('t')
    ylabel('Macrophages');

    figure(3)
    plot(t,y(5,:),'b');
    xlabel('t')
    ylabel('Pro mediators');
    ylim([0 3]);

    figure(4);
    hold on
    plot(t,y(6,:),'r');
    xlabel('t')
    ylabel('Anti mediators');
    ylim([0 40]);

    figure(5);
    plot(t,y(7,:),'b');
    hold on
    plot(t,y(8,:),'r','Color',[0.8500 0.3250 0.0980]);
    xlabel('t')
    ylabel('Hepatocytes');

    figure(6);
    plot(t,y(10,:),'b');
    hold on
    plot(t,y(11,:),'r');
    xlabel('t')
    ylabel('Stellate cells');

    figure(7);
    plot(t,y(9,:),'b');
    xlabel('t')
    ylabel('ECM');

    figure(8);
    cline(y(5,1000:end),y(6,1000:end),y(7,1000:end));

    for i=1:7, figure(i);xlim([0 t(end)]);end
    for i=1:8, figure(i);box on;end

end
