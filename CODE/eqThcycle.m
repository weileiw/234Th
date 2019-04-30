function [Th,Thx,Thxx] = eqThcycle(p,grd,M3d,TRdiv,xth,ith)
x = xth; 
nth = length(ith);
nip = p.nip;   % number of P only parameters;
n = (nth+1)*nth/2;

iwet = find(M3d(:));  nwet = length(iwet);
I    = speye(nwet);
POP  = p.POP;  Px   = p.Px;  Pxx  = p.Pxx;

R1 = p.R1;   R2 = p.R2;   R3 = p.R3;
R4 = p.R4;   R5 = p.R5;   R6 = p.R6;
R7 = p.R7;   R8 = p.R8;   R9 = p.R9;
R10 = p.R10; R11 = p.R11; R12 = p.R12;

PFD_r1 = p.PFD_r1;   dPFDdb_r1  = p.dPFDdb_r1;  d2PFDdb2_r1  = p.d2PFDdb2_r1;
PFD_r2 = p.PFD_r2;   dPFDdb_r2  = p.dPFDdb_r2;  d2PFDdb2_r2  = p.d2PFDdb2_r2;
PFD_r3 = p.PFD_r3;   dPFDdb_r3  = p.dPFDdb_r3;  d2PFDdb2_r3  = p.d2PFDdb2_r3;
PFD_r4 = p.PFD_r4;   dPFDdb_r4  = p.dPFDdb_r4;  d2PFDdb2_r4  = p.d2PFDdb2_r4;
PFD_r5 = p.PFD_r5;   dPFDdb_r5  = p.dPFDdb_r5;  d2PFDdb2_r5  = p.d2PFDdb2_r5;
PFD_r6 = p.PFD_r6;   dPFDdb_r6  = p.dPFDdb_r6;  d2PFDdb2_r6  = p.d2PFDdb2_r6;
PFD_r7 = p.PFD_r7;   dPFDdb_r7  = p.dPFDdb_r7;  d2PFDdb2_r7  = p.d2PFDdb2_r7;
PFD_r8 = p.PFD_r8;   dPFDdb_r8  = p.dPFDdb_r8;  d2PFDdb2_r8  = p.d2PFDdb2_r8;
PFD_r9 = p.PFD_r9;   dPFDdb_r9  = p.dPFDdb_r9;  d2PFDdb2_r9  = p.d2PFDdb2_r9;
PFD_r10 = p.PFD_r10; dPFDdb_r10 = p.dPFDdb_r10; d2PFDdb2_r10 = p.d2PFDdb2_r10;
PFD_r11 = p.PFD_r11; dPFDdb_r11 = p.dPFDdb_r11; d2PFDdb2_r11 = p.d2PFDdb2_r11;
PFD_r12 = p.PFD_r12; dPFDdb_r12 = p.dPFDdb_r12; d2PFDdb2_r12 = p.d2PFDdb2_r12;

u238 = p.u238(iwet);

lambda4 = p.lambda4;
kappa_p = p.kappa_p;
kappa_d = p.kappa_d;
kappa1  = p.kappa1; %kappa_1 in note, adsorption.
kappa2  = p.kappa2; %kappa_(-1) in note, desorption.
aa      = p.aa;
bb      = p.bb;
sigma   = p.sigma;  
POP     = POP(iwet);
DIP     = p.po4obs(iwet);
m2p     = d0(aa+bb*(1-tanh(DIP)));
dm2pdaa = I;
dm2pdbb = d0(1-tanh(DIP));

%dTh_ddt = TRdiv*Th_d-lambda4*(u238-Th_d)-(kappa2+kappa_p)*Th_p+kappa1*POP.*Th_d;
dTh_ddTh_d = TRdiv+lambda4*I+d0(kappa1*m2p*POP);
dTh_ddTh_p = -(kappa2+kappa_p)*I;

%dTh_pdt = PFD*Th_p-kappa1*POP.*Th_d+(kappa2+kappa_p+lambda4)*Th_p;
dTh_pdTh_d = -d0(kappa1*m2p*POP);
dTh_pdTh_p = PFD_r1*d0(R1(iwet))+PFD_r2*d0(R2(iwet))+PFD_r3*d0(R3(iwet))+ ...
    PFD_r4*d0(R4(iwet))+PFD_r5*d0(R5(iwet))+PFD_r6*d0(R6(iwet))+ ...
    PFD_r7*d0(R7(iwet))+PFD_r8*d0(R8(iwet))+PFD_r9*d0(R9(iwet))+ ...
    PFD_r10*d0(R10(iwet))+PFD_r11*d0(R11(iwet))+PFD_r12*d0(R12(iwet))+ ...
    (kappa2+kappa_p+lambda4)*I;

Jac = [[dTh_ddTh_d,dTh_ddTh_p];...
       [dTh_pdTh_d,dTh_pdTh_p]];

rhs = [lambda4*u238;sparse(nwet,1)];

FD = mfactor(Jac);
Th = mfactor(FD,rhs);

if (nargout>1)
    %
    %
    % Compute the gradient of the solution wrt to parameters
    %
    %
    Th_d = Th(1:nwet);
    Th_p = Th(nwet+1:end);
    
    % get derivative 
    POPx = zeros(nwet,nth);
    POPx(:,1:nip) = Px(1+nwet:2*nwet,:);
    dFdPPx = [kappa1*m2p*(d0(Th_d)*POPx); ...
              -kappa1*m2p*(d0(Th_d)*POPx)];
    
    Z = sparse(nwet,1);
    
    for ik1 = 1:length(ith)
        switch ith(ik1)
          case 1 % b_r1
            Fx(:,ik1) = [Z; p.b1*dPFDdb_r1*d0(R1(iwet))*Th_p];
          case 2 
            Fx(:,ik1) = [Z; p.b2*dPFDdb_r2*d0(R2(iwet))*Th_p];
          case 3
            Fx(:,ik1) = [Z; p.b3*dPFDdb_r3*d0(R3(iwet))*Th_p];
          case 4
            Fx(:,ik1) = [Z; p.b4*dPFDdb_r4*d0(R4(iwet))*Th_p];
          case 5
            Fx(:,ik1) = [Z; p.b5*dPFDdb_r5*d0(R5(iwet))*Th_p];
          case 6
            Fx(:,ik1) = [Z; p.b6*dPFDdb_r6*d0(R6(iwet))*Th_p];
          case 7
            Fx(:,ik1) = [Z; p.b7*dPFDdb_r7*d0(R7(iwet))*Th_p];
          case 8
            Fx(:,ik1) = [Z; p.b8*dPFDdb_r8*d0(R8(iwet))*Th_p];
          case 9
            Fx(:,ik1) = [Z; p.b9*dPFDdb_r9*d0(R9(iwet))*Th_p];
          case 10
            Fx(:,ik1) = [Z; p.b10*dPFDdb_r10*d0(R10(iwet))*Th_p];
          case 11
            Fx(:,ik1) = [Z; p.b11*dPFDdb_r11*d0(R11(iwet))*Th_p];
          case 12
            Fx(:,ik1) = [Z; p.b12*dPFDdb_r12*d0(R12(iwet))*Th_p];
          case 13 % kappa_d
            Fx(:,ik1) = [Z;Z];
          case 14 % alpha
            Fx(:,ik1) = [Z;Z];
          case 15 % beta
            Fx(:,ik1) = [Z;Z];
          case 16 % Fdkappa1
            Fx(:,ik1) =  ([kappa1*m2p*POP.*Th_d;...
                           -kappa1*m2p*POP.*Th_d]);
          case 17 % dFdkappa2
            Fx(:,ik1) = ([-kappa2*Th_p; ...
                          kappa2*Th_p]);
          % case 18 % dFdaa
            % Fx(:,ik1) = ([kappa1*dm2pdaa*POP.*Th_d; ...
                          % -kappa1*dm2pdaa*POP.*Th_d]);
          % case 19 % dFdbb
            % Fx(:,ik1) = ([kappa1*dm2pdbb*POP.*Th_d; ...
                          % -kappa1*dm2pdbb*POP.*Th_d]);
        end
    end
    Thx = -mfactor(FD, Fx+dFdPPx);
end

% calculate 2nd derivative.
% (dFdTh)*(Thxx)+...
% (1) (dFdTh)x*(dThdx)+...
% (2) (dFdP)x*(dPdx)+...
% (3) (dFdP)*(dPdx)x+...
% (4) (dFdx)x = 0;

if (nargout>2)
    %
    %
    pos = p.pos;
    % create a 'map' for P only 2nd derivative;
    nt = 1;
    for i1 = 1:nip
        for j1 = i1:nip
            map(nt) = pos(i1,j1);
            nt = nt+1;
        end
    end
    
    Z  = zeros(nwet,1);
    ZZ = zeros(nwet,pos(end,end));
    Th_dx = Thx(1:nwet,:);
    Th_px = Thx(1+nwet:end,:);

    rhsFxx = zeros(2*nwet,pos(end,end));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (1) rhsThx = (dFdTh)x*(dThdx);
    % b1  %%%%%%%%%%%%%
    % vertical part of L
    for i1 = 1:1
        rhsFxx(:,pos(i1,1)) = rhsFxx(:,pos(i1,1))+...
            [Z;...
             p.b1*dPFDdb_r1*d0(R1(iwet))*Th_px(:,i1)];
    end
    % horizontal part of L
    r = [pos(1,1):pos(1,end)];
    rhsFxx(:,r) = rhsFxx(:,r) +...
        [ZZ(:,r);...
         p.b1*dPFDdb_r1*d0(R1(iwet))*Th_px(:,1:end)];

    %%% b2 %%%%%%%%%%
    % vertical part of L
    for i1 = 1:2
        rhsFxx(:,pos(i1,2)) = rhsFxx(:,pos(i1,2))+...
            [Z;...
             p.b2*dPFDdb_r2*d0(R2(iwet))*Th_px(:,i1)];
    end
    % horizontal part of L
    r = [pos(2,2):pos(2,end)];
    rhsFxx(:,r) = rhsFxx(:,r) +...
        [ZZ(:,r);...
         p.b2*dPFDdb_r2*d0(R2(iwet))*Th_px(:,2:end)];

    %%% b3 %%%%%%%%%%
    % vertical part of L
    for i1 = 1:3
        rhsFxx(:,pos(i1,3)) = rhsFxx(:,pos(i1,3))+...
            [Z;...
             p.b3*dPFDdb_r3*d0(R3(iwet))*Th_px(:,i1)];
    end
    % horizontal part of L
    r = [pos(3,3):pos(3,end)];
    rhsFxx(:,r) = rhsFxx(:,r) +...
        [ZZ(:,r);...
         p.b3*dPFDdb_r3*d0(R3(iwet))*Th_px(:,3:end)];

    %%% b4 %%%%%%%%%%
    % vertical part of L
    for i1 = 1:4
        rhsFxx(:,pos(i1,4)) = rhsFxx(:,pos(i1,4))+...
            [Z;...
             p.b4*dPFDdb_r4*d0(R4(iwet))*Th_px(:,i1)];
    end
    % horizontal part of L
    r = [pos(4,4):pos(4,end)];
    rhsFxx(:,r) = rhsFxx(:,r) +...
        [ZZ(:,r);...
         p.b4*dPFDdb_r4*d0(R4(iwet))*Th_px(:,4:end)];

    %%% b5 %%%%%%%%%%
    % vertical part of L
    for i1 = 1:5
        rhsFxx(:,pos(i1,5)) = rhsFxx(:,pos(i1,5))+...
            [Z;...
             p.b5*dPFDdb_r5*d0(R5(iwet))*Th_px(:,i1)];
    end
    % horizontal part of L
    r = [pos(5,5):pos(5,end)];
    rhsFxx(:,r) = rhsFxx(:,r) +...
        [ZZ(:,r);...
         p.b5*dPFDdb_r5*d0(R5(iwet))*Th_px(:,5:end)];

    %%% b6 %%%%%%%%%%
    % vertical part of L
    for i1 = 1:6
        rhsFxx(:,pos(i1,6)) = rhsFxx(:,pos(i1,6))+...
            [Z;...
             p.b6*dPFDdb_r6*d0(R6(iwet))*Th_px(:,i1)];
    end
    % horizontal part of L
    r = [pos(6,6):pos(6,end)];
    rhsFxx(:,r) = rhsFxx(:,r) +...
        [ZZ(:,r);...
         p.b6*dPFDdb_r6*d0(R6(iwet))*Th_px(:,6:end)];

    %%% b7 %%%%%%%%%%
    % vertical part of L
    for i1 = 1:7
        rhsFxx(:,pos(i1,7)) = rhsFxx(:,pos(i1,7))+...
            [Z;...
             p.b7*dPFDdb_r7*d0(R7(iwet))*Th_px(:,i1)];
    end
    % horizontal part of L
    r = [pos(7,7):pos(7,end)];
    rhsFxx(:,r) = rhsFxx(:,r) +...
        [ZZ(:,r);...
         p.b7*dPFDdb_r7*d0(R7(iwet))*Th_px(:,7:end)];

    %%% b8 %%%%%%%%%%
    % vertical part of L
    for i1 = 1:8
        rhsFxx(:,pos(i1,8)) = rhsFxx(:,pos(i1,8))+...
            [Z;...
             p.b8*dPFDdb_r8*d0(R8(iwet))*Th_px(:,i1)];
    end
    % horizontal part of L
    r = [pos(8,8):pos(8,end)];
    rhsFxx(:,r) = rhsFxx(:,r) +...
        [ZZ(:,r);...
         p.b8*dPFDdb_r8*d0(R8(iwet))*Th_px(:,8:end)];

    %%% b9 %%%%%%%%%%
    % vertical part of L
    for i1 = 1:9
        rhsFxx(:,pos(i1,9)) = rhsFxx(:,pos(i1,9))+...
            [Z;...
             p.b9*dPFDdb_r9*d0(R9(iwet))*Th_px(:,i1)];
    end
    % horizontal part of L
    r = [pos(9,9):pos(9,end)];
    rhsFxx(:,r) = rhsFxx(:,r) +...
        [ZZ(:,r);...
         p.b9*dPFDdb_r9*d0(R9(iwet))*Th_px(:,9:end)];

    %%% b10 %%%%%%%%%%
    % vertical part of L
    for i1 = 1:10
        rhsFxx(:,pos(i1,10)) = rhsFxx(:,pos(i1,10))+...
            [Z;...
             p.b10*dPFDdb_r10*d0(R10(iwet))*Th_px(:,i1)];
    end
    % horizontal part of L
    r = [pos(10,10):pos(10,end)];
    rhsFxx(:,r) = rhsFxx(:,r) +...
        [ZZ(:,r);...
         p.b10*dPFDdb_r10*d0(R10(iwet))*Th_px(:,10:end)];

    %%% b11 %%%%%%%%%%
    % vertical part of L
    for i1 = 1:11
        rhsFxx(:,pos(i1,11)) = rhsFxx(:,pos(i1,11))+...
            [Z;...
             p.b11*dPFDdb_r11*d0(R11(iwet))*Th_px(:,i1)];
    end
    % horizontal part of L
    r = [pos(11,11):pos(11,end)];
    rhsFxx(:,r) = rhsFxx(:,r) +...
        [ZZ(:,r);...
         p.b11*dPFDdb_r11*d0(R11(iwet))*Th_px(:,11:end)];

    %%% b12 %%%%%%%%%%
    % vertical part of L
    for i1 = 1:12
        rhsFxx(:,pos(i1,12)) = rhsFxx(:,pos(i1,12))+...
            [Z;...
             p.b12*dPFDdb_r12*d0(R12(iwet))*Th_px(:,i1)];
    end
    % horizontal part of L
    r = [pos(12,12):pos(12,end)];
    rhsFxx(:,r) = rhsFxx(:,r) +...
        [ZZ(:,r);...
         p.b12*dPFDdb_r12*d0(R12(iwet))*Th_px(:,12:end)];

    %%% kappa1 %%%%%%%%%%
    % rhsFxx(:,pos(16,16)) = exp(x(16))*...
    % ([d0(m2p*POP)*Th_dx(:,16);-d0(m2p*POP)*Th_dx(:,16)]);
    % vertical part of L
    for i1 = 1:16
        rhsFxx(:,pos(i1,16)) = rhsFxx(:,pos(i1,16))+...
            [kappa1*m2p*(d0(POP)*Th_dx(:,i1));
             -kappa1*m2p*(d0(POP)*Th_dx(:,i1))];
    end
    % horizontal part of L
    r = [pos(16,16):pos(16,end)];
    rhsFxx(:,r) = rhsFxx(:,r) +...
        [kappa1*m2p*(d0(POP)*Th_dx(:,16:end));
         -kappa1*m2p*(d0(POP)*Th_dx(:,16:end))];

    %%% kappa2 %%%%%%%%%%
    % rhsFxx(:,pos(17,17)) = exp(x(17))*([-Th_px(:,17); Th_px(:,17)]);
    % vertical part of L
    for i1 = 1:17
        rhsFxx(:,pos(i1,17)) = rhsFxx(:,pos(i1,17))+...
            [-kappa2*Th_px(:,i1);...
             kappa2*Th_px(:,i1)];
    end
    % horizontal part of L
    r = [pos(17,17):pos(17,end)];
    rhsFxx(:,r) = rhsFxx(:,r) +...
        [-kappa2*Th_px(:,17:end);
         kappa2*Th_px(:,17:end)];

    %%% aa %%%%%%%%%%
    % rhsFxx(:,pos(18,18)) = x(18)*...
    % ([d0(dm2pdaa*POP)*Th_dx(:,18);-d0(dm2pdaa*POP)*Th_dx(:,18)]);
    % vertical part of L
    % for i1 = 1:18
        % rhsFxx(:,pos(i1,18)) = rhsFxx(:,pos(i1,18))+...
            % [kappa1*dm2pdaa*(d0(POP)*Th_dx(:,i1));...
             % -kappa1*dm2pdaa*(d0(POP)*Th_dx(:,i1))];
    % end
    % horizontal part of L
    % r = [pos(18,18):pos(18,end)];
    % rhsFxx(:,r) = rhsFxx(:,r) +...
        % [kappa1*dm2pdaa*(d0(POP)*Th_dx(:,18:end));
         % -kappa1*dm2pdaa*(d0(POP)*Th_dx(:,18:end))];

    %%% bb %%%%%%%%%%
    % rhsFxx(:,pos(19,19)) = x(19)*...
    % ([d0(dm2pdbb*POP)*Th_dx(:,19);-d0(dm2pdbb*POP)*Th_dx(:,19)]);
    % vertical part of L
    % for i1 = 1:19
        % rhsFxx(:,pos(i1,19)) = rhsFxx(:,pos(i1,19))+...
            % [kappa1*dm2pdbb*(d0(POP)*Th_dx(:,i1));...
             % -kappa1*dm2pdbb*(d0(POP)*Th_dx(:,i1))];
    % end
    % horizontal part of L
    % r = [pos(19,19):pos(19,end)];
    % rhsFxx(:,r) = rhsFxx(:,r) +...
        % [kappa1*dm2pdbb*(d0(POP)*Th_dx(:,19:end));
         % -kappa1*dm2pdbb*(d0(POP)*Th_dx(:,19:end))];

    %%%%%%%%%%%%%% 
    for is = 1:p.nix
        r = [pos(is,is):pos(is,end)];
        rhsFxx(:,r) = rhsFxx(:,r)+...
            [kappa1*m2p*(d0(POPx(:,is))*Th_dx(:,is:end));...
             -kappa1*m2p*(d0(POPx(:,is))*Th_dx(:,is:end))];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (2) (dFdP)x*(dPdx)
    for i1 = 1:16
        rhsFxx(:,pos(i1,16)) = rhsFxx(:,pos(i1,16))+...
            [kappa1*m2p*(d0(Th_d)*POPx(:,i1));...
             -kappa1*m2p*(d0(Th_d)*POPx(:,i1))];
    end

    % for i1 = 1:18
        % rhsFxx(:,pos(i1,18)) = rhsFxx(:,pos(i1,18))+...
            % [kappa1*dm2pdaa*(d0(Th_d)*POPx(:,i1));...
             % -kappa1*dm2pdaa*(d0(Th_d)*POPx(:,i1))];
    % end

    % for i1 = 1:19
        % rhsFxx(:,pos(i1,19)) = rhsFxx(:,pos(i1,19))+...
            % [kappa1*dm2pdbb*(d0(Th_d)*POPx(:,i1));...
             % -kappa1*dm2pdbb*(d0(Th_d)*POPx(:,i1))];
    % end

    %%%%%%%%%%%%%%%%%%%%%%
    for is = 1:p.nix
        r = [pos(is,is):pos(is,end)];
        rhsFxx(:,r) = rhsFxx(:,r)+...
            [kappa1*m2p*(d0(Th_dx(:,is))*POPx(:,is:end));...
             -kappa1*m2p*(d0(Th_dx(:,is))*POPx(:,is:end))];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (3) (dFdP)*(dPdx)x = dFdP*Pxx;
    POPxx = zeros(nwet,n);
    POPxx(:,map) = Pxx(1+nwet:2*nwet,:);
    rhsFxx = rhsFxx+...
             [kappa1*m2p*(d0(Th_d)*POPxx);...
              -kappa1*m2p*(d0(Th_d)*POPxx)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (4) (dFdx)x;
    for ik1 = 1:length(ith)
        switch ith(ik1)
          case 1 % b_r1 b_r1
            rhsFxx(:,pos(ik1,ik1)) = rhsFxx(:,pos(ik1,ik1))+...
                Fx(:,ik1)+[Z; p.b1^2*d2PFDdb2_r1*d0(R1(iwet))*Th_p];
          case 2 % b_r2 b_r2
            rhsFxx(:,pos(ik1,ik1)) = rhsFxx(:,pos(ik1,ik1))+...
                Fx(:,ik1)+[Z; p.b2^2*d2PFDdb2_r2*d0(R2(iwet))*Th_p];
          case 3 % b_r3 b_r3
            rhsFxx(:,pos(ik1,ik1)) = rhsFxx(:,pos(ik1,ik1))+...
                Fx(:,ik1)+[Z; p.b3^2*d2PFDdb2_r3*d0(R3(iwet))*Th_p];
          case 4 % b_r4 b_r4
            rhsFxx(:,pos(ik1,ik1)) = rhsFxx(:,pos(ik1,ik1))+...
                Fx(:,ik1)+[Z; p.b4^2*d2PFDdb2_r4*d0(R4(iwet))*Th_p];
          case 5 % b_r5 b_r5
            rhsFxx(:,pos(ik1,ik1)) = rhsFxx(:,pos(ik1,ik1))+...
                Fx(:,ik1)+[Z; p.b5^2*d2PFDdb2_r5*d0(R5(iwet))*Th_p];
          case 6 % b_r6 b_r6
            rhsFxx(:,pos(ik1,ik1)) = rhsFxx(:,pos(ik1,ik1))+...
                Fx(:,ik1)+[Z; p.b6^2*d2PFDdb2_r6*d0(R6(iwet))*Th_p];
          case 7 % b_r7 b_r7
            rhsFxx(:,pos(ik1,ik1)) = rhsFxx(:,pos(ik1,ik1))+...
                Fx(:,ik1)+[Z; p.b7^2*d2PFDdb2_r7*d0(R7(iwet))*Th_p];
          case 8 % b_r8 b_r8
            rhsFxx(:,pos(ik1,ik1)) = rhsFxx(:,pos(ik1,ik1))+...
                Fx(:,ik1)+[Z; p.b8^2*d2PFDdb2_r8*d0(R8(iwet))*Th_p];
          case 9 % b_r9 b_r9
            rhsFxx(:,pos(ik1,ik1)) = rhsFxx(:,pos(ik1,ik1))+...
                Fx(:,ik1)+[Z; p.b9^2*d2PFDdb2_r9*d0(R9(iwet))*Th_p];
          case 10 % b_r10 b_r10
            rhsFxx(:,pos(ik1,ik1)) = rhsFxx(:,pos(ik1,ik1))+...
                Fx(:,ik1)+[Z; p.b10^2*d2PFDdb2_r10*d0(R10(iwet))*Th_p];
          case 11 % b_r11 b_r11
            rhsFxx(:,pos(ik1,ik1)) = rhsFxx(:,pos(ik1,ik1))+...
                Fx(:,ik1)+[Z; p.b11^2*d2PFDdb2_r11*d0(R11(iwet))*Th_p];
          case 12 % b_r12 b_r12
            rhsFxx(:,pos(ik1,ik1)) = rhsFxx(:,pos(ik1,ik1))+...
                Fx(:,ik1)+[Z; p.b12^2*d2PFDdb2_r12*d0(R12(iwet))*Th_p];
          case 13 % kappa_d kappa_d
            rhsFxx(:,pos(ik1,ik1)) = rhsFxx(:,pos(ik1,ik1))+[Z;Z];
          case 14 % alpha alpha
            rhsFxx(:,pos(ik1,ik1)) = rhsFxx(:,pos(ik1,ik1))+[Z;Z];
          case 15 % beta beta
            rhsFxx(:,pos(ik1,ik1)) = rhsFxx(:,pos(ik1,ik1))+[Z;Z];
          case 16 % kappa1 kappa1
            rhsFxx(:,pos(ik1,ik1)) = rhsFxx(:,pos(ik1,ik1))+...
                Fx(:,ik1);
            % kappa1 aa
            % rhsFxx(:,pos(ik1,ik1+2)) = rhsFxx(:,pos(ik1,ik1+2))+...
                % [kappa1*dm2pdaa*d0(POP)*Th_d;...
                 % -kappa1*dm2pdaa*d0(POP)*Th_d];
            % kappa1 bb
            % rhsFxx(:,pos(ik1,ik1+3)) = rhsFxx(:,pos(ik1,ik1+3))+...
                 % [kappa1*dm2pdbb*d0(POP)*Th_d;...
                  % -kappa1*dm2pdbb*d0(POP)*Th_d];
          case 17 % kappa2 kappa2
            rhsFxx(:,pos(ik1,ik1)) = rhsFxx(:,pos(ik1,ik1))+Fx(:,ik1);
          % case 18 % aa aa
            % rhsFxx(:,pos(ik1,ik1)) = rhsFxx(:,pos(ik1,ik1))+[Z;Z];
          % case 19 % bb bb
            % rhsFxx(:,pos(ik1,ik1)) = rhsFxx(:,pos(ik1,ik1))+[Z;Z];
        end
    end
    
    Thxx = -mfactor(FD, rhsFxx);
    
end