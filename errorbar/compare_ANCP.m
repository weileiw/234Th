addpath('/DFS-L/DATA/primeau/weilewang/DATA/');
addpath('/DFS-L/DATA/primeau/weilewang/my_func/');
load M3d90x180x24v2.mat MSKS
load Reginal_Cexp_90x180.mat
load transport_v4.mat
grd = grid;

TOCexp = smoothit(grd,M3d,TOCexp,3,1e5);
POCexp = smoothit(grd,M3d,POCexp,3,1e5);

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

mean_TOC_tropical = nanmean(TOCexp_tropical(:));
mean_TOC_subtro = nanmean(TOCexp_subtro(:));
mean_TOC_subtro_subpo = nanmean(TOCexp_subtro_subpo(:));
mean_TOC_subpolar = nanmean(TOCexp_subpolar(:));

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
fprintf('subpolar zonal mean DOC to TOC export ratio is %2.2f percent\n', ...
        mean_D2T_subpolar*100)

