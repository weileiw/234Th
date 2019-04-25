function [P,Px,Pxx,L] = eqPcycle(parm,grd,M3d,TRdiv,ip,x)
% output: P is model prediction of DIP,POP,and DOP
% output: D is dPdp, which is used for parameter optimization
% output: D2 is d2Pdp2, which is used for parameter optimization
% M*P = r; dMdp*P+M*dPdp = drdp; dPdp = M\(drdp-dMdp*P);
    
    iwet = find(M3d(:));          % wet point index;
    nwet = length(iwet);        % number of wet points;
    I = speye(nwet);            % make an identity matrix;
    % fixed parameters
    DIPbar = M3d(iwet)*parm.DIPbar;  % gobal arerage PO4 conc.[mmol m^-3]; 
    kappa_g = parm.kappa_g;          % PO4 geological restore const.[s^-1];
    kappa_p = parm.kappa_p;          % POP solubilization rate constant
    sigma   = parm.sigma;
    
    kappa_d = parm.kappa_d;
    alpha   = parm.alpha;
    beta    = parm.beta;
    
    R1 = parm.R1;     R2 = parm.R2;     R3 = parm.R3;
    R4 = parm.R4;     R5 = parm.R5;     R6 = parm.R6;
    R7 = parm.R7;     R8 = parm.R8;     R9 = parm.R9;
    R10 = parm.R10;   R11 = parm.R11;   R12 = parm.R12;
    
    PFD_r1 = parm.PFD_r1; dPFDdb_r1 = parm.dPFDdb_r1; d2PFDdb2_r1 = parm.d2PFDdb2_r1;
    PFD_r2 = parm.PFD_r2; dPFDdb_r2 = parm.dPFDdb_r2; d2PFDdb2_r2 = parm.d2PFDdb2_r2;
    PFD_r3 = parm.PFD_r3; dPFDdb_r3 = parm.dPFDdb_r3; d2PFDdb2_r3 = parm.d2PFDdb2_r3;
    PFD_r4 = parm.PFD_r4; dPFDdb_r4 = parm.dPFDdb_r4; d2PFDdb2_r4 = parm.d2PFDdb2_r4;
    PFD_r5 = parm.PFD_r5; dPFDdb_r5 = parm.dPFDdb_r5; d2PFDdb2_r5 = parm.d2PFDdb2_r5;
    PFD_r6 = parm.PFD_r6; dPFDdb_r6 = parm.dPFDdb_r6; d2PFDdb2_r6 = parm.d2PFDdb2_r6;
    PFD_r7 = parm.PFD_r7; dPFDdb_r7 = parm.dPFDdb_r7; d2PFDdb2_r7 = parm.d2PFDdb2_r7;
    PFD_r8 = parm.PFD_r8; dPFDdb_r8 = parm.dPFDdb_r8; d2PFDdb2_r8 = parm.d2PFDdb2_r8;
    PFD_r9 = parm.PFD_r9; dPFDdb_r9 = parm.dPFDdb_r9; d2PFDdb2_r9 = parm.d2PFDdb2_r9;
    PFD_r10 = parm.PFD_r10; dPFDdb_r10 = parm.dPFDdb_r10; d2PFDdb2_r10 = parm.d2PFDdb2_r10;
    PFD_r11 = parm.PFD_r11; dPFDdb_r11 = parm.dPFDdb_r11; d2PFDdb2_r11 = parm.d2PFDdb2_r11;
    PFD_r12 = parm.PFD_r12; dPFDdb_r12 = parm.dPFDdb_r12; d2PFDdb2_r12 = parm.d2PFDdb2_r12;

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
    
    if (nargout>1)
        %
        %
        % Compute the gradient of the solution wrt the parameters
        %
        %
        Z = zeros(nwet,1);
        DIP = P(1:nwet);
        POP = P(nwet+1:2*nwet);
        DOP = P(2*nwet+1:end);
        Fx = zeros(3*nwet,length(ip));
        % the following code is used to optimize parameters
        % initialize the derivative of the Jacobian wrt to the parameters
        for ik1 = 1:length(ip)
            switch(ip(ik1))
              case 1 % b
                Fx(:,ik1) = exp(x(ik1))*[Z; ...
                                    dPFDdb_r1*d0(R1(iwet))*POP; ...
                                    Z];
              case 2 % b                
                Fx(:,ik1) = exp(x(ik1))*[Z; ...
                                    dPFDdb_r2*d0(R2(iwet))*POP; ...
                                    Z];
              case 3 % b                
                Fx(:,ik1) = exp(x(ik1))*[Z; ...
                                    dPFDdb_r3*d0(R3(iwet))*POP; ...
                                    Z];
              case 4 % b                
                Fx(:,ik1) = exp(x(ik1))*[Z; ...
                                    dPFDdb_r4*d0(R4(iwet))*POP;...
                                    Z];
              case 5 % b
                Fx(:,ik1) = exp(x(ik1))*[Z; ...
                                    dPFDdb_r5*d0(R5(iwet))*POP; ...
                                    Z];
              case 6 % b                
                Fx(:,ik1) = exp(x(ik1))*[Z; ...
                                    dPFDdb_r6*d0(R6(iwet))*POP; ...
                                    Z];
              case 7 % b                
                Fx(:,ik1) = exp(x(ik1))*[Z; ...
                                    dPFDdb_r7*d0(R7(iwet))*POP; ...
                                    Z];
              case 8 % b                
                Fx(:,ik1) = exp(x(ik1))*[Z; ...
                                    dPFDdb_r8*d0(R8(iwet))*POP; ...
                                    Z];
              case 9 % b                
                Fx(:,ik1) = exp(x(ik1))*[Z; ...
                                    dPFDdb_r9*d0(R9(iwet))*POP; ...
                                    Z];
              case 10 % b                
                Fx(:,ik1) = exp(x(ik1))*[Z; ...
                                    dPFDdb_r10*d0(R10(iwet))*POP; ...
                                    Z];
              case 11 % b                
                Fx(:,ik1) = exp(x(ik1))*[Z; ...
                                    dPFDdb_r11*d0(R11(iwet))*POP; ...
                                    Z];
              case 12 % b                
                Fx(:,ik1) = exp(x(ik1))*[Z; ...
                                    dPFDdb_r12*d0(R12(iwet))*POP; ...
                                    Z];
              case 13 % kappa_d
                Fx(:,ik1) =  exp(x(ik1))*[-DOP; ...
                                    Z; ...
                                    DOP];
              case 14 % alpha
                Fx(:,ik1) =  exp(x(ik1))*[L*DIP; ...
                                    -(1-sigma)*L*DIP; ...
                                    -sigma*L*DIP];
              case 15 %beta
                dLambdadbeta = 0*Lambda;
                dLambdadbeta(:,:,1) = log(npp1).*Lambda(:,:,1);
                dLambdadbeta(:,:,2) = log(npp2).*Lambda(:,:,2);
                iz = find(isinf(dLambdadbeta(:)));
                dLambdadbeta(iz) = 0;
                inan = find(isnan(dLambdadbeta(:)));
                dLambdadbeta(inan) = 0;
                dLdbeta = d0(dLambdadbeta(iwet)); 
                Fx(:,ik1) =  exp(x(ik1))*[alpha*dLdbeta*DIP; ...
                                    -(1-sigma)*alpha*dLdbeta*DIP; ...
                                    -sigma*alpha*dLdbeta*DIP];
            end
        end
        % D is the derivative of the solution wrt to the parameters
        Px = mfactor(FJ,-Fx);
    end
    
    if (nargout > 2)
        %
        %
        % Compute the 2nd derivative of the solution wrt the parameters
        %
        %
        
        % initialize the 2nd derivative of the Jacobian wrt the parameters
        for rr = 1:3 % 3 is the number of species DIP,POP,DOP
            for cc = 1:3 % is the number of species DIP,POP,DOP
                for i1 = 1:15 % 5 is the number of parameters, i.e. b,kappa_d,sigma,alpha,beta
                    for i2 = i1:15 % 5 is the number of parameters
                        d2J{rr,cc,i1,i2} = sparse(nwet,nwet);
                    end
                end
            end
        end
        DIPx = Px(1:nwet,:);
        POPx = Px(nwet+1:2*nwet,:);
        DOPx = Px(2*nwet+1:end,:);
        % compute only the upper triangular part of the matrix
        Z=zeros(nwet,length(ip));
        for i1 = 1:length(ip)
            switch ip(i1)
              case 1 % b
                
                FpxPx(:,:,i1) = full(exp(x(i1))*[Z;...
                                    dPFDdb_r1*d0(R1(iwet))*POPx;...
                                    Z]);
                for i2 = i1:length(ip) 
                    switch ip(i2)
                      case 1 % b,b
                        d2J{2,2,1,1} = d2PFDdb2_r1*d0(R1(iwet));
                    end
                end
              case 2 % b
                FpxPx(:,:,i1) = full(exp(x(i1))*[Z;...
                                    dPFDdb_r2*d0(R2(iwet))*POPx;...
                                    Z]);
                for i2 = i1:length(ip) 
                    switch ip(i2)
                      case 2 
                        d2J{2,2,2,2} = d2PFDdb2_r2*d0(R2(iwet));
                    end
                end
              case 3 % b
                FpxPx(:,:,i1) = full(exp(x(i1))*[Z;...
                                    dPFDdb_r3*d0(R3(iwet))*POPx;...
                                    Z]);
                for i2 = i1:length(ip) 
                    switch ip(i2)
                      case 3 
                        d2J{2,2,3,3} = d2PFDdb2_r3*d0(R3(iwet));
                    end
                end
              case 4 % b
                FpxPx(:,:,i1) = full(exp(x(i1))*[Z;...
                                    dPFDdb_r4*d0(R4(iwet))*POPx;...
                                    Z]);
                for i2 = i1:length(ip) 
                    switch ip(i2)
                      case 4
                        d2J{2,2,4,4} = d2PFDdb2_r4*d0(R4(iwet));
                    end
                end
              case 5 % b
                FpxPx(:,:,i1) = full(exp(x(i1))*[Z;...
                                    dPFDdb_r5*d0(R5(iwet))*POPx;...
                                    Z]);
                for i2 = i1:length(ip) 
                    switch ip(i2)
                      case 5
                        d2J{2,2,5,5} = d2PFDdb2_r5*d0(R5(iwet));
                    end
                end
              case 6 % b
                FpxPx(:,:,i1) = full(exp(x(i1))*[Z;...
                                    dPFDdb_r6*d0(R6(iwet))*POPx;...
                                    Z]);
                for i2 = i1:length(ip) 
                    switch ip(i2)
                      case 6
                        d2J{2,2,6,6} = d2PFDdb2_r6*d0(R6(iwet));
                    end
                end
              case 7 % b
                FpxPx(:,:,i1) = full(exp(x(i1))*[Z;...
                                    dPFDdb_r7*d0(R7(iwet))*POPx;...
                                    Z]);
                for i2 = i1:length(ip) 
                    switch ip(i2)
                      case 7
                        d2J{2,2,7,7} = d2PFDdb2_r7*d0(R7(iwet));
                    end
                end
              case 8 % b
                FpxPx(:,:,i1) = full(exp(x(i1))*[Z;...
                                    dPFDdb_r8*d0(R8(iwet))*POPx;...
                                    Z]);
                for i2 = i1:length(ip) 
                    switch ip(i2)
                      case 8
                        d2J{2,2,8,8} = d2PFDdb2_r8*d0(R8(iwet));
                    end
                end
              case 9 % b
                FpxPx(:,:,i1) = full(exp(x(i1))*[Z;...
                                    dPFDdb_r9*d0(R9(iwet))*POPx;...
                                    Z]);
                for i2 = i1:length(ip) 
                    switch ip(i2)
                      case 9
                        d2J{2,2,9,9} = d2PFDdb2_r9*d0(R9(iwet));
                    end
                end
              case 10 % b
                FpxPx(:,:,i1) = full(exp(x(i1))*[Z;...
                                    dPFDdb_r10*d0(R10(iwet))*POPx;...
                                    Z]);
                for i2 = i1:length(ip) 
                    switch ip(i2)
                      case 10
                        d2J{2,2,10,10} = d2PFDdb2_r10*d0(R10(iwet));
                    end
                end
              case 11 % b
                FpxPx(:,:,i1) = full(exp(x(i1))*[Z;...
                                    dPFDdb_r11*d0(R11(iwet))*POPx;...
                                    Z]);
                for i2 = i1:length(ip) 
                    switch ip(i2)
                      case 11
                        d2J{2,2,11,11} = d2PFDdb2_r11*d0(R11(iwet));
                    end
                end
              case 12 % b
                FpxPx(:,:,i1) = full(exp(x(i1))*[Z;...
                                    dPFDdb_r12*d0(R12(iwet))*POPx;...
                                    Z]);
                for i2 = i1:length(ip) 
                    switch ip(i2)
                      case 12
                        d2J{2,2,12,12} = d2PFDdb2_r12*d0(R12(iwet));
                    end
                end
              case 13 % kappa_d
                FpxPx(:,:,i1) = full(exp(x(i1))*[-DOPx;...
                                    Z;...
                                    DOPx]);
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 13 % kappa_d,kappa_d
                      case 14 % kappa_d,alpha
                      case 15 % kappa_d,beta
                    end
                end
              case 14 % alpha
                FpxPx(:,:,i1) = full(exp(x(i1))*[L*DIPx;...
                                    -(1-sigma)*L*DIPx;...
                                    -sigma*L*DIPx]);
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 14 % alpha,alpha
                      case 15 % alpha,beta
                        d2J{1,1,14,15} = dLdbeta;
                        d2J{2,1,14,15} = -(1-sigma)*dLdbeta;
                        d2J{3,1,14,15} =     -sigma*dLdbeta;
                    end
                end
              case 15 % beta
                FpxPx(:,:,i1) = full(exp(x(i1))*[alpha*dLdbeta*DIPx;...
                                    -(1-sigma)*alpha*dLdbeta*DIPx;...
                                    -sigma*alpha*dLdbeta*DIPx]);
                for i2 = i1:length(ip)
                    switch ip(i2)
                      case 15 % beta,beta
                        d2Lambdadbetadbeta = 0*Lambda;
                        d2Lambdadbetadbeta(:,:,1) = log(npp1).*log(npp1).*Lambda(:,:,1);
                        d2Lambdadbetadbeta(:,:,2) = log(npp2).*log(npp2).*Lambda(:,:,2);
                        iz = find(isinf(d2Lambdadbetadbeta(:)));
                        d2Lambdadbetadbeta(iz) = 0;
                        inan = find(isnan(d2Lambdadbetadbeta(:)));
                        d2Lambdadbetadbeta(inan) = 0;
                        d2Ldbetadbeta = d0(d2Lambdadbetadbeta(iwet)); 
                        d2J{1,1,15,15} = alpha*d2Ldbetadbeta;
                        d2J{2,1,15,15} = -(1-sigma)*alpha*d2Ldbetadbeta;
                        d2J{3,1,15,15} =     -sigma*alpha*d2Ldbetadbeta;
                    end
                end
            end
        end
        symcol = @(i,j) (i>=j).*i+(i<j).*j;
        symrow = @(i,j) (i>=j).*j+(i<j).*i;
        delta = @(i,j) (i==j)*1 + (i~=j)*0.0;
        k = 0;
        for i1 = 1:length(ip)
            for i2 = i1:length(ip) 
                k = k+1;
                rhs(:,k) = -(...
                    exp(x(i1)+x(i2))*cell2mat(d2J(:,:,symrow(ip(i1),ip(i2)),symcol(ip(i1),ip(i2))))*P+ ...
                             FpxPx(:,i1,i2)+FpxPx(:,i2,i1)+delta(i1,i2)*Fx(:,i1));
            end
        end
        Pxx = mfactor(FJ,rhs);
    end
end


% --------------------------------------------------------------------

