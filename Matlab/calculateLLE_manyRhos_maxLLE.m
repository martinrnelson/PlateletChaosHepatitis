%% Load the results from mcmc_manyRhos.m, find the rho-values that attained the maximal LLE of those samples, and run local variations around these values to produce the lower-left panels of Figure 9

%% Load datafile

load('MCMC_results.mat');

nParams=size(rhoMatrix,2);

%% Define default parameter set
[LLE_sorted,ind]=sort(LLE,'descend');
rhoMatrix_sorted=rhoMatrix(ind,:);

% Parameters
p=setBaselineParams();
p.rho=rhoMatrix_sorted(1,:);

%% LLE preamble

timeStep=0.5;
maxIntegrationTime=1000;

% Starting ICs
load('rho_ICs','ICs_baseline');
ICs=ICs_baseline([1:8,10]);

%% Range for plots
rhoVec=linspace(0,40,3);%100);
nRhos=numel(rhoVec);

%% Create list of parameter combinations
nCombos=nParams*(nParams-1)/2;
combos=zeros(nCombos,2);
counter=1;
for i=1:nParams-1
    for j=i+1:nParams
        combos(counter,:)=[i j];
        counter=counter+1;
    end
end

%% Initialise arrays
LLE=zeros(nCombos,nRhos,nRhos);
KY=zeros(nCombos,nRhos,nRhos);

parfor parCounter=1:size(combos,1)
    
    param1=combos(parCounter,1);
    param2=combos(parCounter,2); %#ok<*PFBNS>
    LLE_this=zeros(nRhos,nRhos);
    KY_this=zeros(nRhos,nRhos);
    
    tic
    disp(strcat('Running: ',num2str(param1),' ',num2str(param2)));
        
    params_i=p;

    for i=1:nRhos           
        for j=1:nRhos
            
            params_i.rho(param1)=rhoVec(i);
            params_i.rho(param2)=rhoVec(j);
            
            [~,Res]=lyapunov(9,@(t,y)odeSystemLyapunov_manyRhos(t,y,params_i),@ode45,0,timeStep,maxIntegrationTime,ICs,0);
            
            LLE_this(i,j)=max(Res(end,:))
            
            % Compute KY-Dim
            LLElistj=sort(Res(end,:),'descend');
            sumLLEj=zeros(9,1);
            for k=1:size(Res,2)
                sumLLEj(k)=sum(LLElistj(1:k));
            end
            index1=find(sumLLEj>=0);
            index2=find(sumLLEj(2:end)<0);
            KYindex=max(intersect(index1,index2));
            if(KYindex<numel(LLElistj))
                KY_this(i,j)=KYindex+sumLLEj(KYindex)/abs(LLElistj(KYindex+1));
            else
                KY_this(i,j)=0;
            end
            
        end
                    
    end       
    disp(strcat('Complete: ',num2str(param1),' ',num2str(param2)));
    
    LLE(parCounter,:,:)=LLE_this;
    KY(parCounter,:,:)=KY_this;
    toc
        
end

top
save('lyapunovCalc_maxMCMC.mat');


