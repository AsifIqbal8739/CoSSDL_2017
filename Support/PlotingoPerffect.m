% function to generate the plots for paper
% load('GoodData.mat'); close all;
% Ground truth
Disp_Act2(TC,SM,Mask);
suptitle('Ground Truth');
nSub = size(DD,2) - 1;
for i = 1:nSub+1
    XX{i}(abs(XX{i})<0.5) = 0;
end
Comm = size(DD,2);
% % Individual Maps
nTC = size(TC,2);   nV = sqrt(size(Mask,2));
Mask2 = 0.6*double(edge(reshape(Mask,nV,nV)));

figure();
for i = 1:nTC
    jj = (i-1)*3;
    if i <= 3
        CCT = corr(DD{Comm},TC(:,i));      CCS = corr(full(XX{Comm})',SM(i,:)');
        [valT,indT] = max(abs(CCT));    [valS,indS] = max(abs(CCS));
        im = Mask2 + reshape(sign(CCS(indS)).*(XX{7}(indS,:)),nV,nV);
        atom = sign(CCT(indT)).*DD{7}(:,indT);
    else
        CCT = corr(DD{i-3},TC(:,i));      CCS = corr(full(XX{i-3})',SM(i,:)');
        [valT,indT] = max(abs(CCT));    [valS,indS] = max(abs(CCS));
        im = Mask2 + reshape(sign(CCT(indT)).*(XX{i-3}(indT,:)),nV,nV);
        atom = sign(CCT(indT)).*DD{i-3}(:,indT);
    end        
        subplot(nTC,3,jj+1); 
        imagesc(flipud(im));
        title(sprintf('SM {%d}:(%0.3f)',i,valS));  % or use imshow
        set(gca,'YTick',[]); set(gca,'XTick',[]);   colormap jet; axis image
        subplot(nTC,3,jj+[2,3]);
        plot(atom,'LineWidth',1.5);  title(sprintf('TC_{%d}:(%0.3f)',i,valT)); axis tight;
        set(gca,'XTick',[]); 
end

suptitle('Recovered TCs and SMs')



%% Plots for COrrelations
load('Trials100.mat');      % load the .mat file with Mean results here
figure();
imagesc(MeanTC);
set(gca, 'XTick', 1:param.nSub+1,'TickDir','both'); % center x-axis ticks on bins
set(gca, 'YTick', 1:nComp); % center y-axis ticks on bins
Xlabel = {'D_1','D_2','D_3','D_4','D_5','D_6','D_0'};
Ylabel = {'TC_1','TC_2','TC_3','TC_4','TC_5','TC_6','TC_7','TC_8','TC_9'};
set(gca, 'XTickLabel', Xlabel,'FontSize',16); % set x-axis labels
set(gca, 'YTickLabel', Ylabel); % set y-axis labels
colorbar;   colormap hot;
title('Mean TC correlations')

figure();
imagesc(MeanSM);
set(gca, 'XTick', 1:param.nSub+1,'TickDir','both'); % center x-axis ticks on bins
set(gca, 'YTick', 1:nComp); % center y-axis ticks on bins
Xlabel = {'X_1','X_2','X_3','X_4','X_5','X_6','X_0'};
Ylabel = {'S_1','S_2','S_3','S_4','S_5','S_6','S_7','S_8','S_9'};
set(gca, 'XTickLabel', Xlabel,'FontSize',16); % set x-axis labels
set(gca, 'YTickLabel', Ylabel); % set y-axis labels
colorbar;   colormap hot;
title('Mean SM correlations')
% 
