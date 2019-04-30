function [Sum_Cexp,Cexp,Sum_CNPP, TC_HOTS,TC_BATS,TC_OSP,RTOC,fD2T] = ...
    totC_export_errorbar(p,grd,M3d,TRdiv,dVt,ePOC)
    dpa  = 365;  spd  = 24*60^2;   spa  = dpa*spd;
    iwet = find(M3d(:));  nwet = length(iwet);
    I    = speye(nwet);
    sigma = 0.10;   gamma = 0;
    dAt = grd.DXT3d.*grd.DYT3d;  
    DIP = p.DIP(iwet);
    c2p = p.c2p;
    % grid location of HOTS BATS, and BATS
    Lat_HOTS = 22+45/60; Lon_HOTS = mod(-158,360);
    Lat_BATS = 31+40/60; Lon_BATS = mod((-64-10/60),360);
    Lat_OSP  = 50+1/60;  Lon_OSP  = mod((-144-9/60),360);
    
    indx_hots = length(find(grd.xt<Lon_HOTS));
    indy_hots = length(find(grd.yt<Lat_HOTS));

    indx_bats = length(find(grd.xt<Lon_BATS));
    indy_bats = length(find(grd.yt<Lat_BATS));

    indx_osp = length(find(grd.xt<Lon_OSP));
    indy_osp = length(find(grd.yt<Lat_OSP));

    %%%%%%%%%%%%%
    % DIP assimilation
    Lambda = M3d*0;
    Lambda(:,:,1) = (p.npp1.^p.beta).*p.Lambda(:,:,1);
    Lambda(:,:,2) = (p.npp2.^p.beta).*p.Lambda(:,:,2);
    L = p.alpha*d0(Lambda(iwet)); % per second

    % preparation for adjoint method.
    W = d0(dVt(iwet));
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
    CNPP = Int_CNPP*spa*1e-3;
    % convert production from mg C/m^3/s to gC/m^2/year;
    tem_CNPP = CNPP.*dAt(:,:,1)*1e-15;
    Sum_CNPP = nansum(tem_CNPP(:));
    
    %%%%%%%%% -------------- %%%%%%%%%%%%%%
    % calculate total export.
    PFdiv_p  = ...
        p.PFD_r1*d0(p.R1(iwet))+p.PFD_r2*d0(p.R2(iwet))+...
        p.PFD_r3*d0(p.R3(iwet))+p.PFD_r4*d0(p.R4(iwet))+...
        p.PFD_r5*d0(p.R5(iwet))+p.PFD_r6*d0(p.R6(iwet))+...
        p.PFD_r7*d0(p.R7(iwet))+p.PFD_r8*d0(p.R8(iwet))+...
        p.PFD_r9*d0(p.R9(iwet))+p.PFD_r10*d0(p.R10(iwet))+...
        p.PFD_r11*d0(p.R11(iwet))+p.PFD_r12*d0(p.R12(iwet));
    
    F_diag_p = inv(W)*PFdiv_p'*W;
    T_diag   = inv(W)*TRdiv'*W;
    
    junk = M3d;
    junk(:,:,1:nn) = 0;
    Omega = junk(iwet);
    kappa4p = p.kappa_d;
    % adjoint method.
    Jex_P = kappa4p*d0(DIP)*L*...
            (sigma*I+p.kappa_p*(1-sigma-gamma)* ...
             inv(F_diag_p+p.kappa_p*I))*((T_diag+kappa4p*I)\Omega); 
    
    P3d = M3d+nan;
    P3d(iwet) = Jex_P; 
    
    Int_p = 0*M3d(:,:,1);
    
    for ij = 1:nn
        Int_p = Int_p+P3d(:,:,ij).*grd.dzt(ij).*c2p*12;
    end
    
    Cexp = Int_p*spd; % convert P export from mmol P/m^3/s to mg
                      % C/m^2/day;
    tem_Cexp = Cexp.*dAt(:,:,3);
    Sum_Cexp = nansum(tem_Cexp(:))*365*1e-18;

    TOCexp = smoothit(grd,M3d,Cexp,3,1e5);
    
    TC_HOTS = TOCexp(indy_hots,indx_hots)/12/1000*365;
    TC_BATS = TOCexp(indy_bats,indx_bats)/12/1000*365;
    TC_OSP  = TOCexp(indy_osp,indx_osp)/12/1000*365;

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
    % units mol/m2/year;
    RTOC.mean_TOC_tropical = nanmean(TOCexp_tropical(:))/12*365/1000;
    RTOC.mean_TOC_subtro = nanmean(TOCexp_subtro(:))/12*365/1000;
    RTOC.mean_TOC_subtro_subpo = nanmean(TOCexp_subtro_subpo(:))/12*365/1000;
    RTOC.mean_TOC_subpolar = nanmean(TOCexp_subpolar(:))/12*365/1000;
    
    DOCexp = TOCexp-ePOC;
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

    fD2T.mean_D2T_tropical = nanmean(D2T_tropical(:));
    fD2T.mean_D2T_subtro = nanmean(D2T_subtro(:));
    fD2T.mean_D2T_subtro_subpo = nanmean(D2T_subtro_subpo(:));
    fD2T.mean_D2T_subpolar = nanmean(D2T_subpolar(:));

    % fprintf('tropial zonal mean DOC to TOC export ratio is %2.2f percent\n', ...
            % mean_D2T_tropical*100)
    % fprintf('subtropical zonal mean DOC to TOC export ratio is %2.2f percent \n', ...
            % mean_D2T_subtro*100)
    % fprintf('subtropical subpolar zonal mean DOC to TOC export ratio is %2.2f percent \n', ...
            % mean_D2T_subtro_subpo*100)
    % fprintf('subpolar zonal mean DOC to TOC export ratio is %2.2f percent\n', ...
            % mean_D2T_subpolar*100)

