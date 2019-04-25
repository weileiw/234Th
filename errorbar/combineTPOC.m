clc
clear all
close all
addpath('/DFS-L/DATA/primeau/weilewang/DATA/');
addpath('/DFS-L/DATA/primeau/weilewang/my_func/');
load transport_v4.mat 
load teng_regions.mat R  % regions based on Teng et al.
load ../xhat_v2.mat
grd = grid;
% constants define.
p.R1  = M3d*0; tmp = find(R(:) == 1);  p.R1(tmp)  = 1;
p.R2  = M3d*0; tmp = find(R(:) == 2);  p.R2(tmp)  = 1;
p.R3  = M3d*0; tmp = find(R(:) == 3);  p.R3(tmp)  = 1;
p.R4  = M3d*0; tmp = find(R(:) == 4);  p.R4(tmp)  = 1;
p.R5  = M3d*0; tmp = find(R(:) == 5);  p.R5(tmp)  = 1;
p.R6  = M3d*0; tmp = find(R(:) == 6);  p.R6(tmp)  = 1;
p.R7  = M3d*0; tmp = find(R(:) == 7);  p.R7(tmp)  = 1;
p.R8  = M3d*0; tmp = find(R(:) == 8);  p.R8(tmp)  = 1;
p.R9  = M3d*0; tmp = find(R(:) == 9);  p.R9(tmp)  = 1;
p.R10 = M3d*0; tmp = find(R(:) == 10); p.R10(tmp) = 1;
p.R11 = M3d*0; tmp = find(R(:) == 11); p.R11(tmp) = 1;
p.R12 = M3d*0; tmp = find(R(:) == 12); p.R12(tmp) = 1;

b1      = RT.xhat(1);  % Martin curve exponential.
b2      = RT.xhat(2);  % Martin curve exponential.
b3      = RT.xhat(3);  % Martin curve exponential.
b4      = RT.xhat(4);  % Martin curve exponential.
b5      = RT.xhat(5);  % Martin curve exponential.
b6      = RT.xhat(6);  % Martin curve exponential.
b7      = RT.xhat(7);  % Martin curve exponential.
b8      = RT.xhat(8);  % Martin curve exponential.
b9      = RT.xhat(9);  % Martin curve exponential.
b10     = RT.xhat(10); % Martin curve exponential.
b11     = RT.xhat(11); % Martin curve exponential.
b12     = RT.xhat(12); % Martin curve exponential.

load POC.mat
load TOC.mat
POC1 = POC_exp;
TOC1 = TOC_exp;
load POC2.mat
POC2 = POC_exp;
TOC2 = TOC_exp;
grd = grid;

for ji = 1:500
    tmpPOCexp = zeros(90,180);
    tmpPOCexp = smoothit(grd,M3d,POC1(:,:,ji),3,1e5);
    TOCexp(:,:,ji) = smoothit(grd,M3d,TOC1(:,:,ji),3,1e5);
    POCexp(:,:,ji) = ...
        (tmpPOCexp.*p.R1(:,:,3)).*(100/grd.zw(4)).^(-b1)+...
        (tmpPOCexp.*p.R2(:,:,3)).*(100/grd.zw(4)).^(-b2)+...
        (tmpPOCexp.*p.R3(:,:,3)).*(100/grd.zw(4)).^(-b3)+...
        (tmpPOCexp.*p.R4(:,:,3)).*(100/grd.zw(4)).^(-b4)+...
        (tmpPOCexp.*p.R5(:,:,3)).*(100/grd.zw(4)).^(-b5)+...
        (tmpPOCexp.*p.R6(:,:,3)).*(100/grd.zw(4)).^(-b6)+...
        (tmpPOCexp.*p.R7(:,:,3)).*(100/grd.zw(4)).^(-b7)+...
        (tmpPOCexp.*p.R8(:,:,3)).*(100/grd.zw(4)).^(-b8)+...
        (tmpPOCexp.*p.R9(:,:,3)).*(100/grd.zw(4)).^(-b9)+...
        (tmpPOCexp.*p.R10(:,:,3)).*(100/grd.zw(4)).^(-b10)+...
        (tmpPOCexp.*p.R11(:,:,3)).*(100/grd.zw(4)).^(-b11)+...
        (tmpPOCexp.*p.R12(:,:,3)).*(100/grd.zw(4)).^(-b12);

    C_exp_anu = POCexp(:,:,ji).*grd.DXT3d(:,:,3).*grd.DYT3d(:,:,3);
    sum_Cexp(ji) = nansum(C_exp_anu(:))*365*1e-18;
    clear C_exp_anu;
end

for jj = 1:500
    tmpPOCexp = zeros(90,180);
    tmpPOCexp = smoothit(grd,M3d,POC2(:,:,jj),3,1e5);
    TOCexp(:,:,jj+500) = smoothit(grd,M3d,TOC2(:,:,jj),3,1e5);
    
    POCexp(:,:,jj+500) = ...
        (tmpPOCexp.*p.R1(:,:,3)).*(100/grd.zw(4)).^(-b1)+...
        (tmpPOCexp.*p.R2(:,:,3)).*(100/grd.zw(4)).^(-b2)+...
        (tmpPOCexp.*p.R3(:,:,3)).*(100/grd.zw(4)).^(-b3)+...
        (tmpPOCexp.*p.R4(:,:,3)).*(100/grd.zw(4)).^(-b4)+...
        (tmpPOCexp.*p.R5(:,:,3)).*(100/grd.zw(4)).^(-b5)+...
        (tmpPOCexp.*p.R6(:,:,3)).*(100/grd.zw(4)).^(-b6)+...
        (tmpPOCexp.*p.R7(:,:,3)).*(100/grd.zw(4)).^(-b7)+...
        (tmpPOCexp.*p.R8(:,:,3)).*(100/grd.zw(4)).^(-b8)+...
        (tmpPOCexp.*p.R9(:,:,3)).*(100/grd.zw(4)).^(-b9)+...
        (tmpPOCexp.*p.R10(:,:,3)).*(100/grd.zw(4)).^(-b10)+...
        (tmpPOCexp.*p.R11(:,:,3)).*(100/grd.zw(4)).^(-b11)+...
        (tmpPOCexp.*p.R12(:,:,3)).*(100/grd.zw(4)).^(-b12);
    
    C_exp_anu = POCexp(:,:,jj+500).*grd.DXT3d(:,:,3).*grd.DYT3d(:,:,3);
    sum_Cexp(jj+500) = nansum(C_exp_anu(:))*365*1e-18;
    clear C_exp_anu;
end

mPOCexp = nanmean(POCexp,3)*12;
mTOCexp = nanmean(TOCexp,3);
mDOCexp = mTOCexp-mPOCexp;
TPOCexp = mean(sum_Cexp);


figure
subplot(3,1,1)
pcolor(mTOCexp);colorbar;shading flat;
subplot(3,1,2)
pcolor(mPOCexp);colorbar;shading flat;
subplot(3,1,3)
pcolor(mDOCexp);colorbar;shading flat;
save meanCexp mPOCexp mTOCexp mDOCexp

