clc
clear all
close all
addpath('/DFS-L/DATA/primeau/weilewang/DATA')
addpath('/DFS-L/DATA/primeau/weilewang/my_func')
load ../../../add_superlabile_DOP/N2P_Nmodel.mat n2p_exp
load teng_regions.mat R  % regions based on Teng et al.
load Sobs_90x180x24.mat      % woa2013 salinity data.
load tempobs_90x180x24.mat
load Mouw_POC_90x180x24.mat  % sediment trap data MOUW
load npp_90x180.mat
load po4obs_90x180x24.mat
load tempobs_90x180x24.mat
load transport_v4.mat
load Th_DIP_DOP_POP_v2.mat
load xhat_v2.mat
% load tmp_xhat_v2.mat
% RT.xhat = xii;

grd  = grid;
iwet = find(M3d(:));
nwet = length(iwet);
I    = speye(nwet);  
TRdiv = -TR;

dAt = grd.DXT3d.*grd.DYT3d;
dVt = dAt.*grd.DZT3d;

sigma = 0.10;
gamma = 0;

dpa  = 365;  spd  = 24*60^2;   spa  = dpa*spd;

R1  = M3d*0; tmp = find(R(:) == 1);  R1(tmp)  = 1;
R2  = M3d*0; tmp = find(R(:) == 2);  R2(tmp)  = 1;
R3  = M3d*0; tmp = find(R(:) == 3);  R3(tmp)  = 1;
R4  = M3d*0; tmp = find(R(:) == 4);  R4(tmp)  = 1;
R5  = M3d*0; tmp = find(R(:) == 5);  R5(tmp)  = 1;
R6  = M3d*0; tmp = find(R(:) == 6);  R6(tmp)  = 1;
R7  = M3d*0; tmp = find(R(:) == 7);  R7(tmp)  = 1;
R8  = M3d*0; tmp = find(R(:) == 8);  R8(tmp)  = 1;
R9  = M3d*0; tmp = find(R(:) == 9);  R9(tmp)  = 1;
R10 = M3d*0; tmp = find(R(:) == 10); R10(tmp) = 1;
R11 = M3d*0; tmp = find(R(:) == 11); R11(tmp) = 1;
R12 = M3d*0; tmp = find(R(:) == 12); R12(tmp) = 1;

lambda4 = log(2)/24.1;  % Th234 decay time scale [day^-1];
dzt   = grd.dzt;

b1 = RT.xhat(1);   b2 = RT.xhat(2);   b3 = RT.xhat(3);
b4 = RT.xhat(4);   b5 = RT.xhat(5);   b6 = RT.xhat(6);
b7 = RT.xhat(7);   b8 = RT.xhat(8);   b9 = RT.xhat(9);
b10 = RT.xhat(10); b11 = RT.xhat(11); b12 = RT.xhat(12);
kappa4p = RT.xhat(13); % DOP remineralization rate constant. 
alpha   = RT.xhat(14); % linear parameter of npp to DIP assimilation function. 
beta    = RT.xhat(15); % exponential parameter of npp to DIN assimilation function.
kappa_p = 1/(720*60^2);   % POP disolution constant [s^-1];

u238 = 78.6*Sobs-315; % dissolved U-238 (dmp/m^3)
                      % relationship based on Owen et al 2011;
Th234 = Th234d+Th234p;

Def = u238-Th234;
ineg = find(Def(:)<0);
Def(ineg) = 0;

INT(:,:,1) = lambda4*(Def(:,:,1))*grd.dzt(1);
INT(:,:,2) = lambda4*(Def(:,:,2))*grd.dzt(2);
INT(:,:,3) = lambda4*(Def(:,:,3))*grd.dzt(3);

C2Th3 = 135.3*grd.zw(4).^(-0.795);

POCexp = nansum(INT,3)*C2Th3*1e-3*12; %[mg/m2/day]

zonal = nansum(POCexp,2);
C_exp_anu = POCexp.*grd.DXT3d(:,:,3).*grd.DYT3d(:,:,3);
sum_Cexp = nansum(C_exp_anu(:))*365*1e-18;
fprintf('annual POC export is %2.3f \n',sum_Cexp)

%%%%%%%%%% propogate flux to deeper layer %%%%%%%5
fPOC = M3d+nan;
for ji = 4:24
    fPOC(:,:,ji) = ...
        POCexp.*R1(:,:,3).*(nansum(grd.DZT3d(:,:,4:ji),3)/grd.zw(4)).^-b1+...
        POCexp.*R2(:,:,3).*(nansum(grd.DZT3d(:,:,4:ji),3)/grd.zw(4)).^-b2+...
        POCexp.*R3(:,:,3).*(nansum(grd.DZT3d(:,:,4:ji),3)/grd.zw(4)).^-b3+...
        POCexp.*R4(:,:,3).*(nansum(grd.DZT3d(:,:,4:ji),3)/grd.zw(4)).^-b4+...
        POCexp.*R5(:,:,3).*(nansum(grd.DZT3d(:,:,4:ji),3)/grd.zw(4)).^-b5+...
        POCexp.*R6(:,:,3).*(nansum(grd.DZT3d(:,:,4:ji),3)/grd.zw(4)).^-b6+...
        POCexp.*R7(:,:,3).*(nansum(grd.DZT3d(:,:,4:ji),3)/grd.zw(4)).^-b7+...
        POCexp.*R8(:,:,3).*(nansum(grd.DZT3d(:,:,4:ji),3)/grd.zw(4)).^-b8+...
        POCexp.*R9(:,:,3).*(nansum(grd.DZT3d(:,:,4:ji),3)/grd.zw(4)).^-b9+...
        POCexp.*R10(:,:,3).*(nansum(grd.DZT3d(:,:,4:ji),3)/grd.zw(4)).^-b10+...
        POCexp.*R11(:,:,3).*(nansum(grd.DZT3d(:,:,4:ji),3)/grd.zw(4)).^-b11+...
        POCexp.*R12(:,:,3).*(nansum(grd.DZT3d(:,:,4:ji),3)/grd.zw(4)).^-b12;
end
fPOC(:,:,3) = POCexp;
%%%%%%%% compare to sediment trap %%%%%%%%%%%%
ikeep = find(POC_flux(:)>0 & fPOC(:)>0);
loglog(POC_flux(ikeep),fPOC(ikeep),'o','Markerfacecolor','b')
rsquare(POC_flux(ikeep),fPOC(ikeep))
rsquare(log10(POC_flux(ikeep)),log10(fPOC(ikeep)))

hold on
plot([0.1 1000],[0.1 1000],'r','linewidth',2)
hold off
xlim([0.1 1000])
ylim([0.1 1000])
xlabel('POC flux from Mouw et al database (mg/m^2/day)')
ylabel('POC flux from the inverse model (mg/m^2/day)')

%%%%%%%% calculate TOCexp  %%%%%%%%%%%%
c2p1 = R1(:,:,1)*355+R2(:,:,1)*81+...
       R3(:,:,1)*163+R4(:,:,1)*91+...
       R5(:,:,1)*115+R6(:,:,1)*103+...
       R7(:,:,1)*138+R8(:,:,1)*83+...
       R9(:,:,1)*176+R10(:,:,1)*86+...
       R11(:,:,1)*0 +R12(:,:,1)*63;

p2c = 0.0069*po4obs+0.006;

DIP  = DIP(iwet);

inan = find(isnan(npp(:))); 
npp(inan) = min(npp(:));
npp = npp/(12*spd);

npp1 = (1/2)*npp*(1/grd.dzt(1)).*p2c(:,:,1);
npp2 = (1/2)*npp*(1/grd.dzt(2)).*p2c(:,:,2);

Lambda = M3d*0;
Lambda(:,:,1) = 1./(1e-9+po4obs(:,:,1));
Lambda(:,:,2) = 1./(1e-9+po4obs(:,:,2));

junk = M3d;
junk(:,:,4:end) = 0;
isf = find(junk(:));

% DIP assimilation
Lambda(:,:,1) = (npp1.^beta).*Lambda(:,:,1);
Lambda(:,:,2) = (npp2.^beta).*Lambda(:,:,2);
L = alpha*d0(Lambda(iwet)); % per second

% preparation for adjoint method.
W = d0(dVt(iwet));
p.kappa_p = kappa_p;

% n2p = 12.41+6.584*(1-tanh(po4obs)); % based on Wang et al Nature
N2P = n2p_exp.*(M3d(:,:,3)-R11(:,:,3));
c2p2 = N2P.*117/16;
c2p = 0.5*(c2p1+c2p2);
% c2p = c2p2;
nn = 3;
%%%%%%%%%----------------%%%%%%%%%%%%%%
% calculate model primary production.
NPP = M3d*0;
NPP(iwet) = L*DIP; % primary production in P unit.
Int_CNPP = 0*M3d(:,:,1);
Int_PNPP = 0*M3d(:,:,1);

for ij = 1:nn
    Int_CNPP = Int_CNPP+NPP(:,:,ij).*grd.dzt(ij).*c2p*12;
    Int_PNPP = Int_PNPP+NPP(:,:,ij).*grd.dzt(ij); 
end
PNPP = Int_PNPP*spa*1e-3;
CNPP = Int_CNPP*spa*1e-3; % convert production from mg C/m^3/s to g
                         % C/m^2/year;
tem_CNPP = CNPP.*dAt(:,:,1)*1e-15;
Sum_CNPP = nansum(tem_CNPP(:));
fprintf('Model NPP is %3.3e \n',Sum_CNPP);
%%%%%%%%% -------------- %%%%%%%%%%%%%%%

PFD_r1    = buildPFD(M3d,p,grd,b1); 
PFD_r2    = buildPFD(M3d,p,grd,b2);
PFD_r3    = buildPFD(M3d,p,grd,b3);
PFD_r4    = buildPFD(M3d,p,grd,b4);
PFD_r5    = buildPFD(M3d,p,grd,b5);
PFD_r6    = buildPFD(M3d,p,grd,b6);
PFD_r7    = buildPFD(M3d,p,grd,b7);
PFD_r8    = buildPFD(M3d,p,grd,b8);
PFD_r9    = buildPFD(M3d,p,grd,b9);
PFD_r10 = buildPFD(M3d,p,grd,b10);
PFD_r11 = buildPFD(M3d,p,grd,b11);
PFD_r12 = buildPFD(M3d,p,grd,b12);

% calculate total export.
PFdiv_p  = ...
    PFD_r1*d0(R1(iwet))+PFD_r2*d0(R2(iwet))+PFD_r3*d0(R3(iwet))+...
    PFD_r4*d0(R4(iwet))+PFD_r5*d0(R5(iwet))+PFD_r6*d0(R6(iwet))+...
    PFD_r7*d0(R7(iwet))+PFD_r8*d0(R8(iwet))+PFD_r9*d0(R9(iwet))+...
    PFD_r10*d0(R10(iwet))+PFD_r11*d0(R11(iwet))+PFD_r12*d0(R12(iwet));

F_diag_p = inv(W)*PFdiv_p'*W;
T_diag   = inv(W)*TRdiv'*W;

junk = M3d;
junk(:,:,1:nn) = 0;
Omega = junk(iwet);

% adjoint method.
Jex_P = kappa4p*d0(DIP)*L*(sigma*I+kappa_p*(1-sigma-gamma)* ...
        inv(F_diag_p+kappa_p*I))*((T_diag+kappa4p*I)\Omega); 

P3d = M3d+nan;
P3d(iwet) = Jex_P; 

Int_p = 0*M3d(:,:,1);

for ij = 1:nn
    Int_p = Int_p+P3d(:,:,ij).*grd.dzt(ij).*c2p*12;
end
% convert P export from mmol P/m^3/s to mg C/m^2/day;
TOCexp = Int_p*spd; 
tem_Cexp = TOCexp.*dAt(:,:,3);
Sum_Cexp = nansum(tem_Cexp(:))*365*1e-18;
fprintf('Model C export is %3.3e \n\n',Sum_Cexp);

%%%%%%%%%%%%%%%%%%% compare to ANCP %%%%%%%%%%%%%%%%%%%%%%%%%%%%
TOCexp = smoothit(grd,M3d,TOCexp,3,1e5);
POCexp = smoothit(grd,M3d,POCexp,3,1e5);
DOCexp = TOCexp-POCexp;

Lat_HOTS = 22+45/60; Lon_HOTS = mod(-158,360);
Lat_BATS = 31+40/60; Lon_BATS = mod((-64-10/60),360);
Lat_OSP  = 50+1/60;  Lon_OSP  = mod((-144-9/60),360);


indx_hots = length(find(grd.xt<Lon_HOTS));
indy_hots = length(find(grd.yt<Lat_HOTS));

indx_bats = length(find(grd.xt<Lon_BATS));
indy_bats = length(find(grd.yt<Lat_BATS));

indx_osp = length(find(grd.xt<Lon_OSP));
indy_osp = length(find(grd.yt<Lat_OSP));

% find ANCP at specific location and convert unit from
% mg/m2/day to mol/m2/year;
TOCexp_HOTS = TOCexp(indy_hots,indx_hots)/12/1000*365;
TOCexp_BATS = TOCexp(indy_bats,indx_bats)/12/1000*365;
TOCexp_OSP  = TOCexp(indy_osp,indx_osp)/12/1000*365;
fprintf('TOC export at HOT is %2.2f mol/m2/year\n', TOCexp_HOTS)
fprintf('TOC export at BATS is %2.2f mol/m2/year \n',TOCexp_BATS)
fprintf('TOC export at OSP is %2.2f mol/m2/year \n\n', TOCexp_OSP)

msk_tropical = M3d(:,:,1)*0;
msk_tropical(length(find(grd.yt<-15)):length(find(grd.yt<15)),:) = 1;

junk1 = M3d(:,:,1)*0;
junk1(length(find(grd.yt<-30)):length(find(grd.yt<30)),:) = 1;
msk_subtro = junk1-msk_tropical;

junk2  = M3d(:,:,1)*0;
junk2(length(find(grd.yt<-45)):length(find(grd.yt<45)),:) = 1;
msk_subtro_subpo = junk2-junk1;

junk3 =  M3d(:,:,1)*0;
junk3(length(find(grd.yt<-60)):length(find(grd.yt<60)),:) = 1;
msk_subpolar = junk3-junk2;

% units mg/m2/day;
TOCexp_tropical = msk_tropical.*TOCexp;
TOCexp_subtro = msk_subtro.*TOCexp;
TOCexp_subtro_subpo = msk_subtro_subpo.*TOCexp;
TOCexp_subpolar = msk_subpolar.*TOCexp;

TOCexp_tropical(TOCexp_tropical(:)==0) = nan;
TOCexp_subtro(TOCexp_subtro(:)==0) = nan;
TOCexp_subtro_subpo(TOCexp_subtro_subpo(:)==0) = nan;
TOCexp_subpolar(TOCexp_subpolar(:)==0) = nan;

mean_TOC_tropical = nanmean(TOCexp_tropical(:))/12/1000*365;
mean_TOC_subtro = nanmean(TOCexp_subtro(:))/12/1000*365;
mean_TOC_subtro_subpo = nanmean(TOCexp_subtro_subpo(:))/12/1000*365;
mean_TOC_subpolar = nanmean(TOCexp_subpolar(:))/12/1000*365;
fprintf('mean TOC export at tropical is %2.2f mol/m2/yr\n', mean_TOC_tropical)
fprintf('TOC export at subtropical is %2.2f mol/m2/yr \n',mean_TOC_subtro)
fprintf('TOC export at subtropical-subpolar is %2.2f mol/m2/yr \n', mean_TOC_subtro_subpo)
fprintf('TOC export at subpolar is %2.2f mol/m2/year \n\n', mean_TOC_subpolar)

% calculate DOC to TOC export ratio for the four biomes.
D2T = DOCexp./TOCexp;
D2T_tropical = (msk_tropical.*D2T);
D2T_subtro = (msk_subtro.*D2T);
D2T_subtro_subpo = (msk_subtro_subpo.*D2T);
D2T_subpolar = (msk_subpolar.*D2T);

D2T_tropical(D2T_tropical(:)==0) = nan;
D2T_subtro(D2T_subtro(:)==0) = nan;
D2T_subtro_subpo(D2T_subtro_subpo(:)==0) = nan;
D2T_subpolar(D2T_subpolar(:)==0) = nan;

mean_D2T_tropical = nanmean(D2T_tropical(:));
mean_D2T_subtro = nanmean(D2T_subtro(:));
mean_D2T_subtro_subpo = nanmean(D2T_subtro_subpo(:));
mean_D2T_subpolar = nanmean(D2T_subpolar(:));

fprintf('tropial zonal mean DOC to TOC export ratio is %2.2f percent\n', ...
        mean_D2T_tropical*100)
fprintf('subtropical zonal mean DOC to TOC export ratio is %2.2f percent \n', ...
        mean_D2T_subtro*100)
fprintf('subtropical subpolar zonal mean DOC to TOC export ratio is %2.2f percent \n', ...
        mean_D2T_subtro_subpo*100)
fprintf('subpolar zonal mean DOC to TOC export ratio is %2.2f percent\n\n', ...
        mean_D2T_subpolar*100)

