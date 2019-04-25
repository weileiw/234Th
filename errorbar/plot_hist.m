clc
clear all
close all
load EB.mat
sPOC = sort(ERR.sPOC);
idx_LB = floor(length(sPOC)*0.025)+1;
idx_HB = floor(length(sPOC)*0.975);

LB = sPOC(idx_LB);
HB = sPOC(idx_HB);

figure(1)
subplot(2,2,1)
hist(sPOC,50);
hold on
line([LB LB],[0 70],'color','red','linewidth',1.5);
hold on
line([HB HB],[0,70],'color','red','linewidth',1.5);
hold off
xlabel('POC export flux (Pg C yr^-^1)')
fprintf('Median POC is %2.3f Pg C per year \n', median(sPOC))
fprintf('POC export lower bound is %2.3f \n',LB)
fprintf('POC export higher bound is %2.3f \n\n',HB)

subplot(2,2,2)
sTOC = sort(ERR.sTOC);
LB = sTOC(idx_LB); HB = sTOC(idx_HB);
hist(sTOC,50);
hold on
line([LB LB],[0 70],'color','red','linewidth',1.5);
hold on
line([HB HB],[0,70],'color','red','linewidth',1.5);
hold off
xlabel('TOC export flux (Pg C yr^-^1)')
fprintf('Median TOC is %2.3f Pg C per year\n', median(sTOC))
fprintf('TOC export lower bound is %2.3f \n',LB)
fprintf('TOC export higher bound is %2.3f \n\n',HB)

subplot(2,2,3)
sDOC = sTOC-sPOC;  sDOC = sort(sDOC);
LB = sDOC(idx_LB); HB = sDOC(idx_HB);
hist(sDOC,50);
hold on
line([LB LB],[0 70],'color','red','linewidth',1.5);
hold on
line([HB HB],[0,70],'color','red','linewidth',1.5);
hold off
xlabel('DOC export flux (Pg C yr^-^1)')
fprintf('Median DOC is %2.3f Pg C per year\n', median(sDOC))
fprintf('DOC export lower bound is %2.3f \n',LB)
fprintf('DOC export higher bound is %2.3f \n \n',HB)

subplot(2,2,4)
R = sDOC./sTOC; R = sort(R);
LB = R(idx_LB); HB = R(idx_HB);
hist(R,50);
hold on
line([LB LB],[0 70],'color','red','linewidth',1.5);
hold on
line([HB HB],[0,70],'color','red','linewidth',1.5);
hold off
xlabel('DOC export ratio')
fprintf('Median DOC to TOC ratio is %2.3f \n', median(R))
fprintf('DOC to TOC export lower bound is %2.3f \n',LB)
fprintf('DOC to TOC export higher bound is %2.3f \n\n',HB)

TC_HOTS = sort(ERR.TC_HOTS);
LB = TC_HOTS(idx_LB); HB = TC_HOTS(idx_HB);
fprintf('Median TOC export at HOT is %2.3f \n', median(TC_HOTS))
fprintf('TOC export at HOT lower bound is %2.3f \n',LB)
fprintf('TOC export at HOT higher bound is %2.3f \n\n',HB)

TC_BATS = sort(ERR.TC_BATS);
LB = TC_BATS(idx_LB); HB = TC_BATS(idx_HB);
fprintf('Median TOC export at BATS is %2.3f \n', median(TC_BATS))
fprintf('TOC export at BATS lower bound is %2.3f \n',LB)
fprintf('TOC export at BATS higher bound is %2.3f \n\n',HB)

TC_OSP = sort(ERR.TC_OSP);
LB = TC_OSP(idx_LB); HB = TC_OSP(idx_HB);
fprintf('Median TOC export at OSP is %2.3f \n', median(TC_OSP))
fprintf('TOC export at OSP lower bound is %2.3f \n',LB)
fprintf('TOC export at OSP higher bound is %2.3f \n\n',HB)

rTOC = sort(ERR.RTOC_tropical);
LB = rTOC(idx_LB); HB = rTOC(idx_HB);
fprintf('Median DOC:TOC at tropical ocean is %2.3f \n', median(rTOC))
fprintf('TOC lower bound is %2.3f \n',LB)
fprintf('TOC higher bound is %2.3f \n\n',HB)

rTOC = sort(ERR.RTOC_subtro);
LB = rTOC(idx_LB); HB = rTOC(idx_HB);
fprintf('Median TOC export at subtro is %2.3f \n', median(rTOC))
fprintf('TOC lower bound is %2.3f \n',LB)
fprintf('TOC higher bound is %2.3f \n\n',HB)

rTOC = sort(ERR.RTOC_subtro_subpo);
LB = rTOC(idx_LB); HB = rTOC(idx_HB);
fprintf('Median TOC export at subtro-subpolar is %2.3f \n', median(rTOC))
fprintf('TOC lower bound is %2.3f \n',LB)
fprintf('TOC higher bound is %2.3f \n\n',HB)

rTOC = sort(ERR.RTOC_subpolar);
LB = rTOC(idx_LB); HB = rTOC(idx_HB);
fprintf('Median TOC at subpolar is %2.3f \n', median(rTOC))
fprintf('TOC lower bound is %2.3f \n',LB)
fprintf('TOC higher bound is %2.3f \n\n',HB)

fD2T = sort(ERR.fD2T_tropical);
LB = fD2T(idx_LB); HB = fD2T(idx_HB);
fprintf('Median DOC:TOC at tropical ocean is %2.3f \n', median(fD2T))
fprintf('DOC:TOC lower bound is %2.3f \n',LB)
fprintf('DOC:TOC higher bound is %2.3f \n\n',HB)

fD2T = sort(ERR.fD2T_subtro);
LB = fD2T(idx_LB); HB = fD2T(idx_HB);
fprintf('Median TOC export at subtro is %2.3f \n', median(fD2T))
fprintf('DOC:TOC lower bound is %2.3f \n',LB)
fprintf('DOC:TOC higher bound is %2.3f \n\n',HB)

fD2T = sort(ERR.fD2T_subtro_subpo);
LB = fD2T(idx_LB); HB = fD2T(idx_HB);
fprintf('Median TOC export at subtro-subpolar is %2.3f \n', median(fD2T))
fprintf('DOC:TOC lower bound is %2.3f \n',LB)
fprintf('DOC:TOC higher bound is %2.3f \n\n',HB)

fD2T = sort(ERR.fD2T_subpolar);
LB = fD2T(idx_LB); HB = fD2T(idx_HB);
fprintf('Median DOC:TOC at subpolar is %2.3f \n', median(fD2T))
fprintf('DOC:TOC lower bound is %2.3f \n',LB)
fprintf('DOC:TOC higher bound is %2.3f \n',HB)
