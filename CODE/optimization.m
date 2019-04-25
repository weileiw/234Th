clc
clear all
close all
addpath('/DFS-L/DATA/primeau/weilewang/DATA/');
addpath('/DFS-L/DATA/primeau/weilewang/my_func/');
load transport_v4.mat
load M3d90x180x24v2.mat 
load teng_regions.mat R  % regions based on Teng et al.
load npp_90x180.mat       % npp data.
load Sobs_90x180x24.mat   % woa2013 salinity data.
load po4obs_90x180x24.mat % woa2013 phosphate data.
load raw_po4obs_90x180x24.mat
load Ref_Th234t_90x180x24.mat VAR Th234t
MSK = MSKS;
ARC = MSK.ARC;
iARC = find(ARC(:)==1);
Th234t(iARC) = 0;

% constants define.
spa  = 365*24*60^2;     % second per year;
spd  = 24*60^2;         % second per day;

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

p.Th234T = Th234t;

ibad = find(po4obs(iocn)<0.01);
po4obs(iocn(ibad)) = 0.01;
p.po4obs = po4obs;
ibad = find(po4raw(iocn)<0.01);
po4raw(iocn(ibad)) = nan;
p.po4raw = po4raw;
p.po4var = SE2;
% carbon to phosphate ratio;
p.p2c = 0.0069*po4obs+0.006;
inan = find(isnan(npp(:))); 
npp(inan) = min(npp(:));   npp = npp/(12*spd);

p.npp1 = 0.5*npp*(1/grd.dzt(1)).*p.p2c(:,:,1);
p.npp2 = 0.5*npp*(1/grd.dzt(2)).*p.p2c(:,:,2);

p.Lambda = M3d*0;
p.Lambda(:,:,1) = 1./(1e-9+po4obs(:,:,1));
p.Lambda(:,:,2) = 1./(1e-9+po4obs(:,:,2));
% dissolved U-238 (dmp/m^3)
p.u238 = 78.6*Sobs-315; 
% Th234 decay time scale [s^-1];
p.lambda4 = log(2)/(24.1*spd); 

% global averaged DIP conc. [mmol m^-3];
p.DIPbar = nansum(dVt(iocn).*po4obs(iocn))./nansum(dVt(iocn));

p.sigma = 0.10;              % portion of production to DOP[unitless]
p.kappa_p = 1/(720*60^2);   % POP disolution constant [s^-1];
p.kappa_g = 1/(1e6*spa);    % DIP geological restoring constant [s^-1];

b       = 0.76;         % Martin curve exponential [unitless];
alpha   = 2.075e-02;    % npp scaling factor for DIP uptake rate
beta    = 7.884e-01;    % npp scaling exponent for DIP uptake rate
kappa_d = 4.77e-8;      % DOP remineralization constant [s^-1].
kappa1  = 2.881e-05;    % adsorption rate constant [s^-1].
kappa2  = 1.422e-06;    % desorption rate constant [s^-1].
aa      = 1;
bb      = 0;

load xhat_v1.mat
xp = [RT.xhat(1:15)];
x  = [xp; RT.xhat(end-1:end)];
% xp = [b*ones(12,1);kappa_d;alpha;beta]; % P only parameters;
% x = [xp;kappa1;kappa2];  % total parameters;
x0 = log(x);

p.ip  = 1:length(xp);
p.nip = length(xp);

p.ix  = 1:length(x);
p.nix = length(p.ix);
n = (p.nix+1)*p.nix/2;
pos = tril(ones(p.nix),0);
pos(pos==1) = 1:n;
p.pos = pos';

L = @(x) neglogpost(x,p,M3d,TRdiv,grd,dVt);
options = optimoptions(@fminunc, ...
                       'Algorithm','trust-region', ...
                       'GradObj','on', ...
                       'Hessian','on', ...
                       'MaxFunEvals',1000, ...
                       'MaxIter',1000, ...
                       'TolX',1e-10, ...
                       'TolFun',1e-10, ...
                       'DerivativeCheck','off', ...
                       'FinDiffType','central', ...
                       'PrecondBandWidth',Inf);

test_the_code = 0;
if(test_the_code)
    dx = sqrt(-1)*eps.^3*eye(p.nix);
    for ii = 12:17
        x  = real(x0)+dx(:,ii);
        [f,fx,fxx] = neglogpost(x,p,M3d,TRdiv,grd,dVt);
        fprintf('%i %e \n',ii,real(fx(ii))-imag(f)/eps^3);
        fprintf('%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n', ...
                real(fxx(ii,:))-imag(fx(:))'/eps^3);
        fprintf('\n')
    end
    keyboard
end

[xhat,fval,exitflag] = fminunc(L,x0,options);

[f,dfdx,d2fdx2,HH] = neglogpost(xhat,p,M3d,TRdiv,grd,dVt);

fprintf('\n\n')
error = sqrt(diag(inv(HH)));

RT.upbar  = [exp(xhat+error)-exp(xhat)]
RT.lowbar = [exp(xhat)-exp(xhat-error)];

RT.xhat = exp(xhat);
RT.f = f;
RT.d2fdx2  = d2fdx2;
RT.Hessian = HH;

save xhat_v3 RT

name = {'b1';'b2';'b3';'b4';'b5';'b6';'b7';'b8';'b9';'b10';'b11';....
        'b12';'kd';'alpha';'beta';'kappa1';'kappa2'};

xhat = RT.xhat;
upbar = RT.upbar;
lowbar = RT.lowbar;
T = table(xhat,upbar,lowbar,'RowNames',name)

fprintf('----------end------------\n')