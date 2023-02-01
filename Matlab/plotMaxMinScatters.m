%% Plot the maxima and minima of the attractor as a function of rho
% (This code produces Figures 4(c,d).)

clear

tspan=[0 500];

% Starting ICs
load('rho_ICs','ICs_baseline');
ICs=ICs_baseline;

% Parameters
p=setBaselineParams();

% Range of rho values of interest
rhoVec=1:0.025:50;

% Initialise variables
peakArray=cell(size(rhoVec));
minima=zeros(size(rhoVec));
num_pks=0;

% Run each rho value through ODE solver (parallelised for speed)
parfor i=1:numel(rhoVec)
    
    params_i=p;
    params_i.rho=rhoVec(i);
    
    % Remove transient dynamics to converge to attractor
    [transientT,transientY] = ode45(@(t,y)odeSystem_singleRho(t,y,params_i),[0 1000],ICs);
    
    % Integrate around attractor
    sol = ode45(@(t,y)odeSystem_singleRho(t,y,params_i),tspan,transientY(end,:));
    t=sol.x;
    y=sol.y;
    
    % Compute and store maxima and minima
    peakArray{i}=findpeaks(sol.y(5,:));
    minima(i)=min(sol.y(5,:));
    num_peaks=num_pks+numel(peakArray{i});
    
end

% Reshape data for plotting
peakMatrix=zeros(num_pks,2);
counter=0;
for i=1:numel(rhoVec)
   peakMatrix(counter+1:counter+numel(peakArray{i}),1)=rhoVec(i); 
   peakMatrix(counter+1:counter+numel(peakArray{i}),2)=peakArray{i}';
   counter=counter+numel(peakArray{i});
end

% Plot
hold off
plot(peakMatrix(:,1),peakMatrix(:,2),'k.');
hold on
plot(rhoVec,minima,'k.');
drawnow;

xlim([1 max(rhoVec)]);
xlabel('rho');
ylabel('c');


