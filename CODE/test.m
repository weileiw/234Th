%load ~/Dropbox/MOCM/DATA/omega2_ad1e-05_ai1000.mat
%M3d = output.M3d;
%grd = output.grid;
addpath('~/Dropbox/MOCM/WEILEI/myfunc')
load ~/Dropbox/MOCM/DATA/omega2_simple.mat
%load ~/Dropbox/MOCM/DATA/transport_v4.mat
% load ~/Dropbox/Manuscripts/Thorium_export/DATA/totC_exp_PNPP_91x180.mat
load ~/Desktop/Thorium_export/Supp_methods/TOC_POC_v1.mat
grd = grid;
M3d = M3d;

eTOC = TOC_exp;
ePOC = POC_exp;
eDOC = eTOC-ePOC;
ineg = find(eDOC(:)<0);
eDOC(ineg) = 0;

DATA = eDOC./eTOC;

% DATA = eTOC;
fname = sprintf('eDOC_per_91x180');
layer = 'layer_to_be_smoothed: ';
I_layer = input(layer);

tmp = M3d(:,:,I_layer);
iwet = find(tmp(:)==1);

data = DATA;
y = data(iwet);
y = inpaint_nans(y);
data(iwet) = y;

n = 5;
r = 1e5;
Cexp = smooth2d(data,r,n,grd,M3d,I_layer); % increase r and/or n for soother field

contourf(Cexp);colorbar

save(fname,'Cexp')

