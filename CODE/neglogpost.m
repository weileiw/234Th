function [f,fx,fxx,HH] = neglogpost(x,p,M3d,TRdiv,grd,dVt)
nip = p.nip;   nix = p.nix;
xii = exp(x);
save tmp_xhat_v3 xii
for i1 = 1:nix
    fprintf('%3.3e, ',xii(i1));
end
fprintf('\n')

iwet = find(M3d(:));     nwet = length(iwet);

Th234T = p.Th234T;       
Th234T(:,:,10:end) = 0;  

iThT = find(Th234T(iwet)>0);
nThT = length(iThT); 
Th234OBS  = Th234T(iwet(iThT));
% assign 5% relative error based on measurement data.
Th_var = (0.05*Th234T(iwet(iThT))).^2;

DIPOBS = p.po4raw; DIPSE2 = p.po4var; 
iPO4   = find(DIPOBS(iwet)>0 & DIPSE2(iwet)>0);
DIPOBS = DIPOBS(iwet(iPO4));
DIPSE2 = DIPSE2(iwet(iPO4));
nDIP   = length(iPO4);

W1 = d0(dVt(iwet(iPO4))./sum(dVt(iwet(iPO4))))*d0(1./DIPSE2)*(1/nDIP);
W2 = d0(dVt(iwet(iThT))./sum(dVt(iwet(iThT))))*d0(1./Th_var)*(1/nThT);

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
p.aa      = 1;
p.bb      = 0;

% build particle flux divergence.
[p.PFD_r1,p.dPFDdb_r1,p.d2PFDdb2_r1]    = buildPFD(M3d,p,grd,p.b1); 
[p.PFD_r2,p.dPFDdb_r2,p.d2PFDdb2_r2]    = buildPFD(M3d,p,grd,p.b2);
[p.PFD_r3,p.dPFDdb_r3,p.d2PFDdb2_r3]    = buildPFD(M3d,p,grd,p.b3);
[p.PFD_r4,p.dPFDdb_r4,p.d2PFDdb2_r4]    = buildPFD(M3d,p,grd,p.b4);
[p.PFD_r5,p.dPFDdb_r5,p.d2PFDdb2_r5]    = buildPFD(M3d,p,grd,p.b5);
[p.PFD_r6,p.dPFDdb_r6,p.d2PFDdb2_r6]    = buildPFD(M3d,p,grd,p.b6);
[p.PFD_r7,p.dPFDdb_r7,p.d2PFDdb2_r7]    = buildPFD(M3d,p,grd,p.b7);
[p.PFD_r8,p.dPFDdb_r8,p.d2PFDdb2_r8]    = buildPFD(M3d,p,grd,p.b8);
[p.PFD_r9,p.dPFDdb_r9,p.d2PFDdb2_r9]    = buildPFD(M3d,p,grd,p.b9);
[p.PFD_r10,p.dPFDdb_r10,p.d2PFDdb2_r10] = buildPFD(M3d,p,grd,p.b10);
[p.PFD_r11,p.dPFDdb_r11,p.d2PFDdb2_r11] = buildPFD(M3d,p,grd,p.b11);
[p.PFD_r12,p.dPFDdb_r12,p.d2PFDdb2_r12] = buildPFD(M3d,p,grd,p.b12);  
%
%
%%%%%%%%%%%%%%%%%%%%========P model=========%%%%%%%%%%%%%%%%;
%
tic
pos = p.pos;  ip  = p.ip;   nip = p.nip;
[P,Px,Pxx] = eqPcycle(p,grd,M3d,TRdiv,ip,x);
fprintf('Solving P field takes %3.3f min \n',toc/60)

p.Px = Px; p.Pxx = Pxx;
DIP = M3d+nan;   DIP(iwet) = P(1:nwet);
POP = M3d+nan;   POP(iwet) = P(nwet+1:2*nwet);
DOP = M3d+nan;   DOP(iwet) = P(2*nwet+1:end);

p.DIP = DIP;   p.POP = POP;   p.DOP = DOP;
DIPMOD = DIP(iwet(iPO4));

px = zeros(nwet,nix);  px(:,1:nip) = Px(1:nwet,:);  px = px(iPO4,:);

n = (nip+1)*nip/2;
posP = tril(ones(nip),0);  posP(posP==1) = 1:n;  posP = posP';

pxx = zeros(nwet,pos(end,end));
for i1 = 1:nip
    for i2 = i1:nip
        pxx(:,pos(i1,i2)) = Pxx(1:nwet,posP(i1,i2));
    end
end
pxx = pxx(iPO4,:);
%
%
%%%%%%%%%%%%%%%%%%%%========Th model=========%%%%%%%%%%%%%%%%
%
% get All model parameters
tic
ith = 1:length(x);
[Th,Thx,Thxx] = eqThcycle(p,grd,M3d,TRdiv,x,ith);
fprintf('Solving Th field takes %3.3f min \n',toc/60)

Th234d = M3d+nan;            Th234p = M3d+nan;
Th234d(iwet) = Th(1:nwet);   Th234p(iwet) = Th(nwet+1:end);
Th234MOD = Th234d(iwet(iThT)) + Th234p(iwet(iThT));

ThDx  = Thx(1:nwet,:);  ThPx  = Thx(1+nwet:end,:);
ThDxx = Thxx(1:nwet,:); ThPxx = Thxx(nwet+1:end,:);

thdx = ThDx(iThT,:);    thpx = ThPx(iThT,:);    thx  = thdx+thpx;
thdxx = ThDxx(iThT,:);  thpxx = ThPxx(iThT,:);  thxx  = thdxx+thpxx;

%%%%%%%%%%%%====== calculate f fx fxx =======%%%%%%%%%%%%%%%%
% objective function.
ep = DIPMOD-DIPOBS;     et = Th234MOD-Th234OBS;

f = 0.5*((ep.' * W1 * ep) + (et.' * W2 * et));

% calculate gradient;
if (nargout>1)
    fx = ep.' * W1 * px + et.' * W2 * thx;
end

% calculate hessian;
if (nargout>2)
    fxx = zeros(nix,nix);
    for i1 = 1:nix
        for i2 = i1:nix
            fxx(i1,i2) = px(:,i1).'*W1*px(:,i2) + ...
                ep.'*W1*pxx(:,pos(i1,i2))+ ...
                (thx(:,i1).'*W2*thx(:,i2))+ ...
                (et.'*W2*thxx(:,pos(i1,i2)));
            fxx(i2,i1) = fxx(i1,i2);
        end
    end
end
fprintf('objective function value f = %3.3e \n\n',f);

fname4 = sprintf('Th_DIP_DOP_POP_v3');
save(fname4,'Th234d','Th234p','DIP','DOP','POP');

if nargout > 3
    sig = 0.5*f/(nDIP+nThT);
    HH  = fxx/sig;
end
% keyboard
%%%%%%%%%%%%%%%%%%%%======== end =========%%%%%%%%%%%%%%%%