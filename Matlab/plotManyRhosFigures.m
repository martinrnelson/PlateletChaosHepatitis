%% Load datafile
clear
load('MCMC_results.mat');

nIter=size(rhoMatrix,1);
nParams=size(rhoMatrix,2);

%% 2D scatter plots

figure(2);

labels={'A','B','C','D','E','F','G'};

[LLE_sorted,ind]=sort(LLE,'descend');
rhoMatrix_sorted=rhoMatrix(ind,:);

% Which steps to plot
percentToPlot=0.533951;
numToPlot=floor((percentToPlot/100)*nIter);

for i=1:nParams
    for j=i+1:nParams      
        subplot(nParams,nParams,nParams*(i-1)+j);
        hold on;
        plot(rhoMatrix_sorted(1:numToPlot,i),rhoMatrix_sorted(1:numToPlot,j),'.','Color',[.7 .7 .7]);
        
        plot(rhoMatrix_sorted(1,i),rhoMatrix_sorted(1,j),'ks','MarkerFaceColor','k');
        xlabel(labels{i});
        ylabel(labels{j});
        box on;
    end
end

% Which steps to plot
percentToPlot=0.05;
numToPlot=floor((percentToPlot/100)*nIter);

for i=1:nParams
    for j=i+1:nParams      
        subplot(nParams,nParams,nParams*(i-1)+j);
        hold on;
        plot(rhoMatrix_sorted(1:numToPlot,i),rhoMatrix_sorted(1:numToPlot,j),'r.');
        
        plot(rhoMatrix_sorted(1,i),rhoMatrix_sorted(1,j),'ks','MarkerFaceColor','k');
%         xlabel(labels{i});
%         ylabel(labels{j});
        box on;
        axis equal
        xlim([0 40]);
        ylim([0 40]);
        set(gca,'XTick',[0 40]);
        set(gca,'YTick',[0 40]);
    end
end

%% LLE colour plots for max LLE baseline params

clear;
load('lyapunovCalc_maxMCMC.mat');
load('colormaps.mat')
labels={'A','B','C','D','E','F','G'};

caxis_LLE=[-0.01 0.03];

figure(2);

for i=1:nParams
    for j=i+1:nParams      
        subplot(nParams,nParams,nParams*(j-1)+i);
        hold on;
        
        thisSlice_idx=find(combos(:,1)==i & combos(:,2)==j);
        thisSlice=reshape(LLE(thisSlice_idx,:,:),nRhos,nRhos);
        
        imagesc(rhoVec,rhoVec,thisSlice');
        set(gca,'YDir','normal');
        caxis(caxis_LLE)
        colormap('hot');
%         colormap(cmap_LLE);
        xlabel(labels{i});
        ylabel(labels{j});
        box on
        set(gca,'XTick',[0 40]);
        set(gca,'YTick',[0 40]);
        axis equal
        xlim([0 40]);
        ylim([0 40]);
        
        %plot(rhoMatrix_sorted(1,i),rhoMatrix_sorted(1,j),'ws','MarkerFaceColor','w');
    end
end
