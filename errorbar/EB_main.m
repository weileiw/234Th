clc
clear all
close all
addpath('/DFS-L/DATA/primeau/weilewang/DATA/');
addpath('/DFS-L/DATA/primeau/weilewang/my_func/');
load transport_v4.mat
load teng_regions.mat R  % regions based on Teng et al.
load npp_90x180.mat       % npp data.
load Sobs_90x180x24.mat   % woa2013 salinity data.
load po4obs_90x180x24.mat % woa2013 phosphate data.
load N2P_Nmodel.mat n2p_exp
load ../xhat_v2.mat
N2P = n2p_exp;
% constants define.
spa  = 365*24*60^2;     % second per year;
spd  = 24*60^2;         % second per day;
iwet = find(M3d(:));
nwet = length(iwet);
TRdiv = -TR;     % transport operator (s^-1).
grd   = grid;    % OCIM grid.

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

p.po4obs = po4obs;
% carbon to phosphate ratio;
p.p2c = 0.0069*po4obs+0.006;
inan = find(isnan(npp(:))); 
npp(inan) = min(npp(:));   npp = npp/(12*spd);

p.npp1 = (1/2)*npp*(1/grd.dzt(1)).*p.p2c(:,:,1);
p.npp2 = (1/2)*npp*(1/grd.dzt(2)).*p.p2c(:,:,2);

p.Lambda = M3d*0;
p.Lambda(:,:,1) = 1./(1e-9+po4obs(:,:,1));
p.Lambda(:,:,2) = 1./(1e-9+po4obs(:,:,2));
% dissolved U-238 (dmp/m^3)
p.u238 = 78.6*Sobs-315;
% Th234 decay time scale [s^-1];
p.lambda4 = log(2)/(24.1*spd); 

% global averaged DIP conc. [mmol m^-3];
p.DIPbar = nansum(dVt(iwet).*po4obs(iwet))./nansum(dVt(iwet));

p.PO4 = po4obs(iwet);       % po4 observation [mmol m^-3];

p.sigma = 0.10;             % portion of production to DOP[unitless]
p.kappa_p = 1/(720*60^2);   % POP disolution constant [s^-1];
p.kappa_g = 1/(1e6*spa);    % DIP geological restoring constant [s^-1];

iter = 1000;
POC_exp = zeros(90,180,iter);
TOC_exp = zeros(90,180,iter);
w_all = sparse(nwet,iter);
POP_all = sparse(nwet,iter);

R = R(:,:,3);
r1 = find(R(:)==1);
r2 = find(R(:)==2);
r3 = find(R(:)==3);
r4 = find(R(:)==4);
r5 = find(R(:)==5);
r6 = find(R(:)==6);
r7 = find(R(:)==7);
r8 = find(R(:)==8);
r9 = find(R(:)==9);
r10 = find(R(:)==10);
r11 = find(R(:)==11);
r12 = find(R(:)==12);

for ji = 1:10
    if mod(ji,100) == 0
        fprintf('current iteration is %d ....:-) \n',ji);
    end
    x0 = log(RT.xhat);
    Hessian = inv(RT.Hessian);
    Hessian = (Hessian+Hessian.')/2;
    x0 = mvnrnd(x0,Hessian);
    xii = exp(x0);
    p.b1      = xii(1);  % Martin curve exponential.
    p.b2      = xii(2);  % Martin curve exponential.
    p.b3      = xii(3);  % Martin curve exponential.
    p.b4      = xii(4);  % Martin curve exponential.
    p.b5      = xii(5);  % Martin curve exponential.
    p.b6      = xii(6);  % Martin curve exponential.
    p.b7      = xii(7);  % Martin curve exponential.
    p.b8      = xii(8);  % Martin curve exponential.
    p.b9      = xii(9);  % Martin curve exponential.
    p.b10     = xii(10); % Martin curve exponential.
    p.b11     = xii(11); % Martin curve exponential.
    p.b12     = xii(12); % Martin curve exponential.
    p.kappa_d = xii(13); % DOP remineralization rate constant. 
    p.alpha   = xii(14); % linear parameter of npp to DIP assimilation function. 
    p.beta    = xii(15); % exponential parameter of npp to DIN assimilation function.
    p.kappa1  = xii(16); % adsoprtion rate constant.
    p.kappa2  = xii(17); % desorption rate constant.
    p.aa = 1;   p.bb = 0;

    % c2p1 is based on Teng et al;
    c2p  = zeros(90,180,2);
    c2p(:,:,1) = ...
        p.R1(:,:,1)*(355+(rand*2-1)*59)+...
        p.R2(:,:,1)*(081+(rand*2-1)*18)+...
        p.R3(:,:,1)*(163+(rand*2-1)*40)+...
        p.R4(:,:,1)*(091+(rand*2-1)*09)+...
        p.R5(:,:,1)*(115+(rand*2-1)*35)+...
        p.R6(:,:,1)*(103+(rand*2-1)*26)+...
        p.R7(:,:,1)*(138+(rand*2-1)*33)+...
        p.R8(:,:,1)*(083+(rand*2-1)*13)+...
        p.R9(:,:,1)*(176+(rand*2-1)*30)+...
        p.R10(:,:,1)*(086+(rand*2-1)*20)+...
        p.R11(:,:,1)*(000+(rand*2-1)*00)+...
        p.R12(:,:,1)*(063+(rand*2-1)*20);

    % c2p1 is based on Wang et al;
    N2P = N2P.*(M3d(:,:,1)-p.R11(:,:,1));
    c2p(:,:,2) = N2P*(117/16);
    p.c2p = nanmean(c2p,3);

    [p.PFD_r1,w1]  = buildPFD(M3d,p,grd,p.b1); 
    [p.PFD_r2,w2]  = buildPFD(M3d,p,grd,p.b2);
    [p.PFD_r3,w3]  = buildPFD(M3d,p,grd,p.b3);
    [p.PFD_r4,w4]  = buildPFD(M3d,p,grd,p.b4);
    [p.PFD_r5,w5]  = buildPFD(M3d,p,grd,p.b5);
    [p.PFD_r6,w6]  = buildPFD(M3d,p,grd,p.b6);
    [p.PFD_r7,w7]  = buildPFD(M3d,p,grd,p.b7);
    [p.PFD_r8,w8]  = buildPFD(M3d,p,grd,p.b8);
    [p.PFD_r9,w9]  = buildPFD(M3d,p,grd,p.b9);
    [p.PFD_r10,w10] = buildPFD(M3d,p,grd,p.b10);
    [p.PFD_r11,w11] = buildPFD(M3d,p,grd,p.b11);
    [p.PFD_r12,w12] = buildPFD(M3d,p,grd,p.b12);

    mean_w = nanmean(w1(:,:,2:3),3).*p.R1(:,:,3)+...
             nanmean(w2(:,:,2:3),3).*p.R2(:,:,3)+...
             nanmean(w3(:,:,2:3),3).*p.R3(:,:,3)+...
             nanmean(w4(:,:,2:3),3).*p.R4(:,:,3)+...
             nanmean(w5(:,:,2:3),3).*p.R5(:,:,3)+...
             nanmean(w6(:,:,2:3),3).*p.R6(:,:,3)+...
             nanmean(w7(:,:,2:3),3).*p.R7(:,:,3)+...
             nanmean(w8(:,:,2:3),3).*p.R8(:,:,3)+...
             nanmean(w9(:,:,2:3),3).*p.R9(:,:,3)+...
             nanmean(w10(:,:,2:3),3).*p.R10(:,:,3)+...
             nanmean(w11(:,:,2:3),3).*p.R11(:,:,3)+...
             nanmean(w12(:,:,2:3),3).*p.R12(:,:,3);;
    mean_w = -mean_w*60*60*24; % m/day

    [p.DIP,POP,p.ThM] = neglogpost(p,grd,M3d,TRdiv,dVt);
    
    mPOP = nanmean(POP(:,:,1:3),3);
    ePOP = POP.*mean_w;
    
    [sPOC,ePOC] = POC_export_errorbar(p,M3d,grd);
    ERR.sPOC(ji) = sPOC;
    POC_exp(:,:,ji) = ePOC/12;
    % particulate POC (from 234Th) to POP (P model) ratio
    c2p = (ePOC/12)./ePOP;

    ibad = find(c2p>500 | c2p<0);
    c2p(ibad) = nan;

    ERR.c2p_R1(ji) = nanmean(c2p(r1));
    ERR.c2p_R2(ji) = nanmean(c2p(r2));
    ERR.c2p_R3(ji) = nanmean(c2p(r3));
    ERR.c2p_R4(ji) = nanmean(c2p(r4));
    ERR.c2p_R5(ji) = nanmean(c2p(r5));
    ERR.c2p_R6(ji) = nanmean(c2p(r6));
    ERR.c2p_R7(ji) = nanmean(c2p(r7));
    ERR.c2p_R8(ji) = nanmean(c2p(r8));
    ERR.c2p_R9(ji) = nanmean(c2p(r9));
    ERR.c2p_R10(ji) = nanmean(c2p(r10));
    ERR.c2p_R11(ji) = nanmean(c2p(r11));
    ERR.c2p_R12(ji) = nanmean(c2p(r12));
    
    [sTOC,eTOC,NPP,TC_HOTS,TC_BATS,TC_OSP,RTOC,fD2T] = ...
        totC_export_errorbar(p,grd,M3d,TRdiv,dVt,ePOC);

    ERR.sTOC(ji) = sTOC;
    ERR.NPP(ji) = NPP;
    TOC_exp(:,:,ji) = eTOC;

    ERR.TC_HOTS(ji) = TC_HOTS;
    ERR.TC_BATS(ji) = TC_BATS;
    ERR.TC_OSP(ji) = TC_OSP;

    ERR.RTOC_tropical(ji) = RTOC.mean_TOC_tropical;
    ERR.RTOC_subtro(ji) = RTOC.mean_TOC_subtro;
    ERR.RTOC_subtro_subpo(ji) = RTOC.mean_TOC_subtro_subpo;
    ERR.RTOC_subpolar(ji) = RTOC.mean_TOC_subpolar;
    
    ERR.fD2T_tropical(ji) = fD2T.mean_D2T_tropical;
    ERR.fD2T_subtro(ji) = fD2T.mean_D2T_subtro;
    ERR.fD2T_subtro_subpo(ji) = fD2T.mean_D2T_subtro_subpo;
    ERR.fD2T_subpolar(ji) = fD2T.mean_D2T_subpolar;
    
    save tmp ERR
    % save TOC2 TOC_exp
    % save POC2 POC_exp
end
fprintf('-----END------\n')