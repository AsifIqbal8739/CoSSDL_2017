% To display all the activations and time series from a dictionary and
% Sparse Coefficient matrix
% D = Dictionary;   X = Coeff Matrix
function [h] = Disp_Act2(D,X,Mask)

%% making activations positive
ss = sign(mean(X,2));
X = diag(ss) * X;
D = D * diag(ss);

Mask = 0.6*double(edge(reshape(Mask,sqrt(size(Mask,2)),sqrt(size(Mask,2)))));

n = size(D,2);      
p = sqrt(size(X,2));
h = figure('Position',[800,400,600,500]);
% h = figure();
m = 1:3*n;
for i = 1:n
    s = m((3*(i-1))+1:(3*(i-1))+3);
    subplot(n,3,s(1));
    imagesc(Mask + reshape((X(i,:)),p,p)); title(strcat('Source ',num2str(i)));  % or use imshow
    set(gca,'YTick',[]); set(gca,'XTick',[]);   colormap jet; axis image
    subplot(n,3,s(2:3));
    plot(D(:,i),'LineWidth',1.5);  title(strcat('TC_',num2str(i))); axis tight;
    %set(gca,'YTick',[]); 
    set(gca,'XTick',[]);  
end
end