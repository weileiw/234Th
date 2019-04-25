clc
clear all
close all
load errorbar.mat
ERR1 = ERR;
load errorbar2.mat
ERR2 = ERR;
clear  ERR 
ERR.NPP  = [ERR1.NPP(1:500), ERR2.NPP];
ERR.sPOC = [ERR1.sPOC(1:500), ERR2.sPOC];
ERR.sTOC = [ERR1.sTOC(1:500), ERR2.sTOC];

ERR.c2p_R1 = [ERR1.c2p_R1(1:500), ERR2.c2p_R1];
ERR.c2p_R2 = [ERR1.c2p_R2(1:500), ERR2.c2p_R2];
ERR.c2p_R3 = [ERR1.c2p_R3(1:500), ERR2.c2p_R3];
ERR.c2p_R4 = [ERR1.c2p_R4(1:500), ERR2.c2p_R4];
ERR.c2p_R5 = [ERR1.c2p_R5(1:500), ERR2.c2p_R5];
ERR.c2p_R6 = [ERR1.c2p_R6(1:500), ERR2.c2p_R6];
ERR.c2p_R7 = [ERR1.c2p_R7(1:500), ERR2.c2p_R7];
ERR.c2p_R8 = [ERR1.c2p_R8(1:500), ERR2.c2p_R8];
ERR.c2p_R9 = [ERR1.c2p_R9(1:500), ERR2.c2p_R9];
ERR.c2p_R10 = [ERR1.c2p_R10(1:500), ERR2.c2p_R10];
ERR.c2p_R11 = [ERR1.c2p_R11(1:500), ERR2.c2p_R11];
ERR.c2p_R12 = [ERR1.c2p_R12(1:500), ERR2.c2p_R12];

ERR.TC_HOTS = [ERR1.TC_HOTS(1:500), ERR2.TC_HOTS];
ERR.TC_BATS = [ERR1.TC_BATS(1:500), ERR2.TC_BATS];
ERR.TC_OSP  = [ERR1.TC_OSP(1:500), ERR2.TC_OSP];
ERR.RTOC_tropical = [ERR1.RTOC_tropical(1:500), ...
                    ERR2.RTOC_tropical];

ERR.RTOC_subtro = [ERR1.RTOC_subtro(1:500), ERR2.RTOC_subtro];
ERR.RTOC_subtro_subpo = [ERR1.RTOC_subtro_subpo(1:500), ...
                    ERR2.RTOC_subtro_subpo];

ERR.RTOC_subpolar = [ERR1.RTOC_subpolar(1:500),ERR2.RTOC_subpolar];
ERR.fD2T_tropical = [ERR1.fD2T_tropical(1:500),ERR2.fD2T_tropical];
ERR.fD2T_subtro   = [ERR1.fD2T_subtro(1:500),ERR2.fD2T_subtro];

ERR.fD2T_subtro_subpo = ...
    [ERR1.fD2T_subtro_subpo(1:500), ...
     ERR2.fD2T_subtro_subpo];

ERR.fD2T_subpolar = [ERR1.fD2T_subpolar(1:500),ERR2.fD2T_subpolar];
