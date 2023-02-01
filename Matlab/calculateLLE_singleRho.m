clear

timeStep=0.5;
maxIntegrationTime=1000;

% Starting ICs
load('rho_ICs','ICs_baseline');
ICs=ICs_baseline([1:8,10]);

% ICs=rand(9,1);
% ICs=ICs.*(ICs>0);

% Parameters
p=setBaselineParams();
p.phi=0.01;
p.nu2=0.01;

% filename='LLE_baseline_rho_nu2.mat';
% rhoVec=1:0.1:30;
% nu2Vec=0.01:0.01:1;

% filename='LLE_baseline_rho_phi.mat';
% rhoVec=1:0.1:30;
% phiVec=0.001:0.001:0.04;

filename='LLE_baseline_rho_km.mat';
rhoVec=1:0.1:30;
kmVec=0.00001:0.00001:0.001;


nParam1s=numel(rhoVec);
nParam2s=numel(kmVec);

LLE=zeros(nParam1s,nParam2s);
KY=zeros(nParam1s,nParam2s);

parfor i=1:nParam1s
   
   params_i=p;
   params_i.rho=rhoVec(i);
   
   LLEj=zeros(1,nParam2s);
   KYj=zeros(1,nParam2s);
    
   for j=1:nParam2s
            
        params_i.km=kmVec(j); %#ok<PFBNS>

        [T,Res]=lyapunov(9,@(t,y)odeSystemLyapunov_singleRho(t,y,params_i),@ode45,0,timeStep,maxIntegrationTime,ICs,0);
%         plot(T,Res);
%         title('Dynamics of Lyapunov exponents');
%         xlabel('Time'); ylabel('Lyapunov exponents');

        LLEj(j)=max(Res(end,:));
        
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
            KYj(j)=KYindex+sumLLEj(KYindex)/abs(LLElistj(KYindex+1));
        else
            KYj(j)=0;
        end
        
   end
   
   LLE(i,:)=LLEj(:);
   KY(i,:)=KYj(:);
end

save(filename);

%% Example code to plot results
% 
% load(filename)
%
% caxis_LLE=[-0.01 0.03];
% 
% figure(5);
% imagesc(rhoVec,kmVec,LLE'); % or phiVec or nu2Vec in place of kmVec
% set(gca,'YDir','normal');
% caxis(caxis_LLE);
% colormap('jet');
% % colorbar;
% xlabel('r');
% ylabel('k');
% ylim([0 max(kmVec)]);