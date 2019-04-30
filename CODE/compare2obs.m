clc
clear all
close all
addpath('/DFS-L/DATA/primeau/weilewang/DATA/');
addpath('/DFS-L/DATA/primeau/weilewang/my_func');
addpath('/DFS-L/SCRATCH/moore/weilewang/get_dot_mat')
load transport_v4.mat
load raw_po4obs_90x180x24.mat
load po4obs_90x180x24.mat
load Ref_Th234t_90x180x24.mat
load M3d90x180x24v2.mat 
load Th_DIP_DOP_POP_v2.mat
grd = grid;
iwet = find(M3d(:));
nwet = length(iwet);
MSK = MSKS;
ARC = MSK.ARC;
iARC = find(ARC(:)==1);

Th234T = Th234t; 
Th234T(:,:,10:end) = 0;
Th234T(iARC) = 0;

% ikeep = find(Th234T(iwet)>0);
% pd = prctile(Th234T(iwet(ikeep)),[0.25 97.5]);
% ismall = find(Th234T(iwet(ikeep))<=pd(1)); Th234T(iwet(ikeep(ismall))) = nan;
% ibig   = find(Th234T(iwet(ikeep))>=pd(2)); Th234T(iwet(ikeep(ibig)))   = nan;

% % convert unit from mBq/kg(Th230) to dpm/m^3; 
ThM = Th234d+Th234p;
% ikeep = find(ThM(iwet)>0);
% pd = prctile(ThM(iwet(ikeep)),[0.25 97.5]);
% ismall = find(ThM(iwet(ikeep))<=pd(1)); ThM(iwet(ikeep(ismall))) = nan;
% ibig   = find(ThM(iwet(ikeep))>=pd(2)); ThM(iwet(ikeep(ibig)))   = nan;

Th234mod = ThM/61.7;
Th234obs = Th234T/61.7;

iThT = find(Th234obs(iwet)>0 & Th234mod(iwet)>0);
Th_var = (Th234obs(iwet(iThT))*0.05).^2;
W = (dVt(iwet(iThT))./sum(dVt(iwet(iThT)))).*(1./Th_var);

cr = 5:5:95;
data = [Th234obs(iwet(iThT)),Th234mod(iwet(iThT))];
O = Th234obs(iwet(iThT));
M = Th234mod(iwet(iThT));
rsquare(O,M)
[bandwidth,density,X,Y] = mykde2d(data,200,[15 15],[50 50],W);

figure(1)
dx = X(3,5)-X(3,4); 
dy = Y(4,2)-Y(3,2);
[q,ii] = sort(density(:)*dx*dy,'descend');
D = density;
D(ii) = cumsum(q);
subplot('position',[0.2 0.2 0.6 0.6])
contourf(X,Y,100*(1-D),cr); hold on
contour(X,Y,100*(1-D),cr);

caxis([5 95])
%set(gca,'FontSize',16);
grid on
axis square
xlabel('Observed ^2^3^4Th (mBq/L)');
ylabel('Model ^2^3^4Th (mBq/L)');
% title('model V.S. observation')
plot([15 50],[15 50],'r--','linewidth',2);

subplot('position',[0.82 0.2 0.05 0.6]);
contourf([1 2],cr,[cr(:),cr(:)],cr); hold on
contour([1 2],cr,[cr(:),cr(:)],cr);
hold off
%set(gca,'FontSize',14);
set(gca,'XTickLabel',[]);
set(gca,'YAxisLocation','right');
set(gca,'TickLength',[0 0])
ylabel('(percentile)')

%%%%%%%%%%%%%%%%% compare DIP  %%%%%%%%%%%%%%
figure(2)
DIPobs = po4raw; 
% DIPobs(:,:,10:end) = 0;
ipo4 = find(DIPobs(iwet)>0);

% DIP_var = SE2(iwet(ipo4));
W = (dVt(iwet(ipo4))./sum(dVt(iwet(ipo4))));

cr = 5:5:95;
data = [DIPobs(iwet(ipo4)),DIP(iwet(ipo4))];
O = DIPobs(iwet(ipo4));
M = DIP(iwet(ipo4));
rsquare(O,M)
[bandwidth,density,X,Y] = mykde2d(data,200,[0 0],[4 4],W);

dx = X(3,5)-X(3,4); 
dy = Y(4,2)-Y(3,2);
[q,ii] = sort(density(:)*dx*dy,'descend');
D = density;
D(ii) = cumsum(q);
subplot('position',[0.2 0.2 0.6 0.6])
contourf(X,Y,100*(1-D),cr); hold on
contour(X,Y,100*(1-D),cr);

caxis([5 95])
%set(gca,'FontSize',16);
grid on
axis square
xlabel('Observed DIP (mmol/m^3)');
ylabel('Model DIP (mmol/m^3)');
% title('model V.S. observation')
plot([0 4],[0 4],'r--','linewidth',2);

subplot('position',[0.82 0.2 0.05 0.6]);
contourf([1 2],cr,[cr(:),cr(:)],cr); hold on
contour([1 2],cr,[cr(:),cr(:)],cr);
hold off
%set(gca,'FontSize',14);
set(gca,'XTickLabel',[]);
set(gca,'YAxisLocation','right');
set(gca,'TickLength',[0 0])
ylabel('(percentile)')


