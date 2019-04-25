function [DIP,POP,ThM] = neglogpost(p,grd,M3d,TRdiv,dVt)

iwet = find(M3d(:));     nwet = length(iwet);
%%%%%%%%%%%%%%%%%%%%========P model=========%%%%%%%%%%%%%%%%;
%
%
[P] = eqPcycle(p,grd,M3d,TRdiv);

DIP = M3d+nan;   DIP(iwet) = P(1:nwet);
POP = M3d+nan;   POP(iwet) = P(nwet+1:2*nwet);
DOP = M3d+nan;   DOP(iwet) = P(2*nwet+1:end);

p.DIP = DIP;   p.POP = POP;   p.DOP = DOP;
%
%%%%%%%%%%%%%%%%%%%%========Th model=========%%%%%%%%%%%%%%%%
%
% get All model parameters
[Th] = eqThcycle(p,grd,M3d,TRdiv);

Th234d = M3d+nan;            Th234p = M3d+nan;
Th234d(iwet) = Th(1:nwet);   Th234p(iwet) = Th(nwet+1:end);
ThM = Th234d + Th234p;

function [P] = eqPcycle(parm,grd,M3d,TRdiv);

iwet = find(M3d(:));        % wet point index;
nwet = length(iwet);        % number of wet points;
I = speye(nwet);            % make an identity matrix;
                            % fixed parameters
DIPbar = M3d(iwet)*parm.DIPbar;% gobal arerage PO4 conc.[mmol m^-3]; 
kappa_g = parm.kappa_g;     % PO4 geological restore const.[s^-1];
kappa_p = parm.kappa_p;     % POP solubilization rate constant
sigma   = parm.sigma;

kappa_d = parm.kappa_d;
alpha   = parm.alpha;
beta    = parm.beta;

R1 = parm.R1;     R2 = parm.R2;     R3 = parm.R3;
R4 = parm.R4;     R5 = parm.R5;     R6 = parm.R6;
R7 = parm.R7;     R8 = parm.R8;     R9 = parm.R9;
R10 = parm.R10;   R11 = parm.R11;   R12 = parm.R12;

PFD_r1 = parm.PFD_r1; 
PFD_r2 = parm.PFD_r2; 
PFD_r3 = parm.PFD_r3; 
PFD_r4 = parm.PFD_r4; 
PFD_r5 = parm.PFD_r5; 
PFD_r6 = parm.PFD_r6; 
PFD_r7 = parm.PFD_r7; 
PFD_r8 = parm.PFD_r8; 
PFD_r9 = parm.PFD_r9; 
PFD_r10 = parm.PFD_r10;
PFD_r11 = parm.PFD_r11;
PFD_r12 = parm.PFD_r12;

npp1 = parm.npp1;   npp2 = parm.npp2;
Lambda = parm.Lambda;
Lambda(:,:,1) = (npp1.^beta).*Lambda(:,:,1);
Lambda(:,:,2) = (npp2.^beta).*Lambda(:,:,2);
L = d0(Lambda(iwet));  % PO4 assimilation rate [s^-1];

% build Jacobian matrix
% +++++++++++++++++++++++++++++++++++++++++
% column 1 (d/dDIP)
J{1,1} = TRdiv+alpha*L+kappa_g*I;
J{2,1} = -(1-sigma)*alpha*L;
J{3,1} = -sigma*alpha*L;

% column 2 (d/dPOP)
J{1,2} = 0*I;
J{2,2} = PFD_r1*d0(R1(iwet))+PFD_r2*d0(R2(iwet))+PFD_r3*d0(R3(iwet))+...
         PFD_r4*d0(R4(iwet))+PFD_r5*d0(R5(iwet))+PFD_r6*d0(R6(iwet))+...
         PFD_r7*d0(R7(iwet))+PFD_r8*d0(R8(iwet))+PFD_r9*d0(R9(iwet))+...
         PFD_r10*d0(R10(iwet))+PFD_r11*d0(R11(iwet))+PFD_r12*d0(R12(iwet))+...
         +kappa_p*I;
% J{2,2} = PFD+kappa_p*I;
J{3,2} = -kappa_p*I;

% column 3 (d/dDOP)
J{1,3} = -kappa_d*I;
J{2,3} = 0*I;
J{3,3} = TRdiv+kappa_d*I;
% ++++++++++++++++++++++++++++++++++++++++

% right hand side of phosphate equations 
RHS{1,1} = DIPbar*kappa_g;
RHS{2,1} = sparse(nwet,1);
RHS{3,1} = sparse(nwet,1);

% factoring Jacobian matrix
Jac = cell2mat(J);
FJ = mfactor(Jac); 
% solve solutions 
P = mfactor(FJ,cell2mat(RHS));

function [Th] = eqThcycle(parm,grd,M3d,TRdiv);

iwet = find(M3d(:));  nwet = length(iwet);
I    = speye(nwet);
POP  = parm.POP;

R1 = parm.R1;   R2 = parm.R2;   R3 = parm.R3;
R4 = parm.R4;   R5 = parm.R5;   R6 = parm.R6;
R7 = parm.R7;   R8 = parm.R8;   R9 = parm.R9;
R10 = parm.R10; R11 = parm.R11; R12 = parm.R12;

PFD_r1  = parm.PFD_r1;
PFD_r2  = parm.PFD_r2;
PFD_r3  = parm.PFD_r3;
PFD_r4  = parm.PFD_r4;
PFD_r5  = parm.PFD_r5;
PFD_r6  = parm.PFD_r6;
PFD_r7  = parm.PFD_r7;
PFD_r8  = parm.PFD_r8;
PFD_r9  = parm.PFD_r9;
PFD_r10 = parm.PFD_r10;
PFD_r11 = parm.PFD_r11;
PFD_r12 = parm.PFD_r12;

u238 = parm.u238(iwet);

lambda4 = parm.lambda4;
kappa_p = parm.kappa_p;
kappa_d = parm.kappa_d;
kappa1  = parm.kappa1; %kappa_1 in note, adsorption.
kappa2  = parm.kappa2; %kappa_(-1) in note, desorption.
aa      = parm.aa;
bb      = parm.bb;
sigma   = parm.sigma;  
POP     = parm.POP(iwet);
DIP     = parm.DIP(iwet);
m2p     = d0(aa+bb*(1-tanh(DIP)));

%dTh_ddt = TRdiv*Th_d-lambda4*(u238-Th_d)-(kappa2+kappa_p)*Th_p+kappa1*POPARM.*Th_d;
dTh_ddTh_d = TRdiv+lambda4*I+d0(kappa1*m2p*POP);
dTh_ddTh_p = -(kappa2+kappa_p)*I;

%dTh_pdt = PFD*Th_p-kappa1*POPARM.*Th_d+(kappa2+kappa_p+lambda4)*Th_p;
dTh_pdTh_d = -d0(kappa1*m2p*POP);
dTh_pdTh_p = ...
    PFD_r1*d0(R1(iwet))+PFD_r2*d0(R2(iwet))+PFD_r3*d0(R3(iwet))+ ...
    PFD_r4*d0(R4(iwet))+PFD_r5*d0(R5(iwet))+PFD_r6*d0(R6(iwet))+ ...
    PFD_r7*d0(R7(iwet))+PFD_r8*d0(R8(iwet))+PFD_r9*d0(R9(iwet))+ ...
    PFD_r10*d0(R10(iwet))+PFD_r11*d0(R11(iwet))+PFD_r12*d0(R12(iwet))+ ...
    (kappa2+kappa_p+lambda4)*I;

Jac = [[dTh_ddTh_d,dTh_ddTh_p];...
       [dTh_pdTh_d,dTh_pdTh_p]];

rhs = [lambda4*u238;sparse(nwet,1)];

FD = mfactor(Jac);
Th = mfactor(FD,rhs);
