%% Metropolis-Hastings algorithm to explore 7D rho-space, seeking maximal LLE.
% (Produces upper-right panels of Figure 9)
%
% Inputs:
%   nSteps_in = number of steps to perform
%   filename = datafile storing a batch of previous steps, in order
%              to augment with further steps (see line 38 for details)


function [rhoMatrix,LLE,nAccepts,nSteps,p]=mcmc_manyRhos(nSteps_in,varargin)

    % Lyapunov parameters
    timeStep=0.5;
    maxIntegrationTime=100;

    % MCMC parameters
    stepSize=10;

    % Starting ICs
    load('rho_ICs','ICs_baseline');
    ICs=ICs_baseline([1:8,10]);

    % Parameters
    p=setBaselineParams();

    % To start from scratch
    if(nargin==1)
        nSteps=nSteps_in;
        rhoMatrix=zeros(nSteps,7);
        LLE=zeros(nSteps,1);
        nAccepts=0;
        
        iStart=1;
        p.rho=2*ones(1,7);
        rhoMatrix(1,:)=p.rho;
        [~,Res]=lyapunov(9,@(t,y)odeSystemLyapunov_manyRhos(t,y,p),@ode45,0,timeStep,maxIntegrationTime,ICs,0);
        LLE(1)=max(Res(end,:));
    else
        % To extend existing run
        filename=varargin{1};
        varsToLoad={'rhoMatrix','LLE','nAccepts','nSteps','p'};
        dataLoaded=load(filename,varsToLoad{:});
        
        nSteps=dataLoaded.nSteps+nSteps_in;
        rhoMatrix=zeros(nSteps,7);
        LLE=zeros(nSteps,1);
        
        p=dataLoaded.p;
        rhoMatrix(1:dataLoaded.nSteps,:)=dataLoaded.rhoMatrix;
        LLE(1:dataLoaded.nSteps)=dataLoaded.LLE;
        nAccepts=dataLoaded.nAccepts;
        iStart=dataLoaded.nSteps;
    end
    
    for i=iStart+1:nSteps

        proposed=p;
        proposed.rho=p.rho+stepSize*randn(1,7);

        % Constrain the search to [0,40]^7 
        if(min(proposed.rho)<0 || max(proposed.rho)>40)
            rhoMatrix(i,:)=p.rho;
            LLE(i)=LLE(i-1);         
            continue
        end

        [~,Res]=lyapunov(9,@(t,y)odeSystemLyapunov_manyRhos(t,y,proposed),@ode45,0,timeStep,maxIntegrationTime,ICs,0);
        LLE_proposed=max(Res(end,:));

        % Accept or reject
        if(LLE_proposed>LLE(i-1) || log(LLE_proposed)-log(LLE(i-1))<log(rand(1)))
            rhoMatrix(i,:)=proposed.rho;
            LLE(i)=LLE_proposed;
            p=proposed;
            nAccepts=nAccepts+1;
        else
            rhoMatrix(i,:)=p.rho;
            LLE(i)=LLE(i-1);        
        end

    end
    
end