function [sum_Cexp,C_exp_3] = POC_export_errorbar(p,M3d,grd)

u238 = p.u238;
iwet = find(M3d(:));
lambda4 = log(2)/24.1;  % Th234 decay time scale [day^-1];

dzt   = grd.dzt;

% relationship based on Owen et al 2011;
Def = u238-p.ThM;
ineg = find(Def(:)<0);
Def(ineg) = 0;

INT(:,:,1) = lambda4*(Def(:,:,1))*grd.dzt(1);
INT(:,:,2) = lambda4*(Def(:,:,2))*grd.dzt(2);
INT(:,:,3) = lambda4*(Def(:,:,3))*grd.dzt(3);

sigma = 0.25;
C2Th3 = 135.3*grd.zw(4).^(-0.795);
C2Th3 = normrnd(C2Th3,sigma);

C_exp_3 = nansum(INT,3)*C2Th3*1e-3*12; %[mg/m2/day]

zonal = nansum(C_exp_3,2);
C_exp_anu = C_exp_3.*grd.DXT3d(:,:,3).*grd.DYT3d(:,:,3);
sum_Cexp = nansum(C_exp_anu(:))*365*1e-18;


