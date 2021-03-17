function [Cstar, CTEs, As] = GMC(Constits, RUC)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Copyright 2020 United States Government as represented by the Administrator of the 
% National Aeronautics and Space Administration. No copyright is claimed in the 
% United States under Title 17, U.S. Code. All Other Rights Reserved. BY DOWNLOADING 
% OR USING THIS SOFTWARE, YOU ACKNOWLEDGE THAT YOU HAVE READ THE NASA OPEN SOURCE 
% AGREEMENT V1.3, THAT YOU UNDERSTAND IT, AND THAT YOU AGREE TO BE BOUND BY ITS 
% TERMS. IF YOU DO NOT AGREE TO THE TERMS AND CONDITIONS OF THIS AGREEMENT, DO NOT 
% USE OR DOWNLOAD THE SOFTWARE. THIS SOFTWARE IS PROVIDED AS IS WITHOUT ANY WARRANTY 
% OF ANY KIND. RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST, AND INDEMNIFIES 
% AND HOLDS HARMLESS, THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND 
% SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT. This code was prepared by Drs. 
% B.A. Bednarcyk and S.M. Arnold to complement the book “Practical Micromechanics of 
% Composite Materials” during the course of their government work.
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
% Purpose: Calculate the effective stiffness, effective CTES, and the mechanical and 
%          thermal strain concentration tensors for the GMC micromechanics theory 
% Input:
% - Constits: Struct containing constituent material information
% - RUC: Struct containing repeating unit cell (RUC) information
% Output:
% - Cstar: Effective stiffness matrix
% - CTEs: Effective CTEs
% - As: Struct containing subcell concentration tensors
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -- Number of subcells in each direction
NB = RUC.NB;
NG = RUC.NG;

% -- Obtain computational time for large RUCs
if NB*NG > 25*25
    disp(['  ** Getting Effective Properties using GMC - RUC = ', num2str(NB), ' by ', num2str(NG)])
    tic
end

% -- Assign stiffness and CTE of each subcell
props(NB, NG).C = zeros(6); % -- Preallocate props
for IB = 1:NB
    for IG = 1:NG
        switch(RUC.matsCh(IB, IG))
            case 'F'
                constit = Constits.Fiber;
            case '1'
                constit = Constits.FiberDam;
            case 'M'
                constit = Constits.Matrix;
            case '2'
                constit = Constits.MatrixDam;
            case 'I'
                constit = Constits.Interface;
            case '3'
                constit = Constits.InterfaceDam;
            otherwise
                error(['Invalid RUC.matsch in GMC - RUC.matsch = ', RUC.matsch])
        end
        props(IB, IG).C = GetCFromProps(constit);
        props(IB, IG).alpha = [constit.aL; constit.aT; constit.aT; 0; 0; 0];
    end
end

% -- Subcell dimensions
h = RUC.h;
l = RUC.l;
sumH = sum(h);
sumL = sum(l);

% -- Assemble AG and J (Eq. 4.59) as SPARSE matrices
NNZ = 60*NB*NG; % -- Estimate of upper bound of non-zeros in AG
AG = zeros(NNZ,1);
J = zeros(2*(NB + NG)+ NB*NG + 1, 6);
row = 0;
ic = 0;

% -- ES11(B,G) = EG11 (Eq. 5.29)
for IB = 1:NB
    for IG = 1:NG
        row = row + 1;
        J(row, 1) = 1;
        NS = NG*(IB - 1) + IG;
        ic = ic + 1;
        NR(ic) = row; % -- Row number in AG
        NC(ic) = 6*(NS - 1) + 1; % -- Column number number in AG
        AG(ic) = 1; % --Value number in AG
    end
end

% -- SUM [H(B)*ES22(B,G)]=H*EG22 (Eq. 5.43)
for IG = 1:NG
    row = row + 1;
    J(row, 2) = sumH;
    for IB = 1:NB
        NS = NG*(IB - 1) + IG;
        ic = ic + 1;
        NR(ic) = row;
        NC(ic) = 6*(NS - 1) + 2;
        AG(ic) = h(IB);
    end
end

% -- SUM [L(G)*ES33(B,G)]=L*EG33  (Eq. 5.44)
for IB = 1:NB
    row = row + 1;
    J(row, 3) = sumL;
    for IG = 1:NG
        NS = NG*(IB - 1) + IG;
        ic = ic + 1;
        NR(ic) = row;
        NC(ic) = 6*(NS - 1) + 3;
        AG(ic) = l(IG);
    end
end

% -- SUM SUM [H(B)*L(G)*ES23(B,G)]=H*L*EG23 (Eq. 5.48)
row = row + 1;
J(row, 4) = sumH*sumL;
for IB = 1:NB
    for IG = 1:NG
        NS = NG*(IB - 1) + IG;
        ic = ic + 1;
        NR(ic) = row;
        NC(ic) = 6*(NS - 1) + 4;
        AG(ic) = h(IB)*l(IG);
    end
end 

% -- SUM [L(G)*ES13(B,G)]=L*EG13 (5.46)
for IB = 1:NB
    row = row + 1;
    J(row, 5) = sumL;
    for IG = 1:NG
        NS = NG*(IB - 1) + IG;
        ic = ic + 1;
        NR(ic) = row;
        NC(ic) = 6*(NS - 1) + 5;
        AG(ic) = l(IG);
    end
end

% -- SUM [H(B)*ES12(B,G)]=H*EG12 (Eq. 5.47)
for IG = 1:NG
    row = row + 1;
    J(row, 6) = sumH;
    for IB = 1:NB
        NS = NG*(IB - 1) + IG;
        ic = ic + 1;
        NR(ic) = row;
        NC(ic) = 6*(NS - 1) + 6;
        AG(ic) = h(IB);
    end
end

% -- Define AG in matlab sparse format
AGsparse = sparse( NR(1:ic), NC(1:ic), AG(1:ic) );

% -- Assemble Am (Eq. 5.60)   
NNZ = 6*NB*NG; % -- Estimate of upper bound of non-zeros in AM
AM = zeros(NNZ,1);

% -- Reset counters
row = 0;
ic = 0;

 % -- sig22(B,G) = sig22(B+1,G) (Eq. 5.53)
 for IG = 1:NG
    for IB = 1:NB - 1
        row = row + 1;
        for i = 1:6
            NS = NG*(IB - 1) + IG;
            ic = ic + 1;
            NR(ic) = row;
            NC(ic) = 6*(NS - 1) + i;
            AM(ic) = props(IB, IG).C(2, i);
            NS = NG*(IB) + IG;
            ic = ic + 1;
            NR(ic) = row;
            NC(ic)=6*(NS - 1) + i;
            AM(ic)= -props(IB + 1, IG).C(2, i);
        end
    end
 end

% -- sig33(B,G) = sig33(B,G+1) (Eq. 5.54)
for IB = 1:NB
    for IG = 1:NG - 1
        row = row + 1;
        for i = 1:6
            NS = NG*(IB - 1) + IG;
            ic = ic + 1;
            NR(ic) = row;
            NC(ic)=6*(NS - 1) + i;
            AM(ic) = props(IB, IG).C(3, i);
            NS = NG*(IB - 1) + IG + 1;
            ic = ic + 1;
            NR(ic) = row;
            NC(ic)=6*(NS - 1) + i;
            AM(ic) = -props(IB, IG + 1).C(3, i);
        end
    end
end

% -- sig23(B,G) = sig23(B+1,G) (Eq. 5.55)   
for IG = 1:NG
    for IB = 1:NB - 1
        row = row + 1;
        for i = 1:6
            NS = NG*(IB - 1) + IG;
            ic = ic + 1;
            NR(ic) = row;
            NC(ic) = 6*(NS - 1) + i;
            AM(ic) = props(IB, IG).C(4, i);
            NS = NG*(IB) + IG;
            ic = ic + 1;
            NR(ic) = row;
            NC(ic)=6*(NS - 1) + i;
            AM(ic)= -props(IB + 1, IG).C(4, i);
        end
    end
end

% -- sig32(B,G) = sig32(B,G+1) (Eq. 5.56)
IB = NB;
for IG = 1:NG - 1
    row = row + 1;
    for i = 1:6
        NS = NG*(IB - 1) + IG;
        ic = ic + 1;
        NR(ic) = row;
        NC(ic)=6*(NS - 1) + i;
        AM(ic) = props(IB, IG).C(4, i);
        NS = NG*(IB - 1) + IG + 1;
        ic = ic + 1;
        NR(ic) = row;
        NC(ic)=6*(NS - 1) + i;
        AM(ic) = -props(IB, IG + 1).C(4, i);
    end
end 

% -- sig31(B,G) = sig31(B,G+1) (eq. 5.57)
for IB = 1:NB
    for IG = 1:NG - 1
        row = row + 1;
        for i = 1:6
            NS = NG*(IB - 1) + IG;
            ic = ic + 1;
            NR(ic) = row;
            NC(ic)=6*(NS - 1) + i;
            AM(ic) = props(IB, IG).C(5, i); 
            NS = NG*(IB - 1) + IG + 1;
            ic = ic + 1;
            NR(ic) = row;
            NC(ic)=6*(NS - 1) + i;
            AM(ic) = -props(IB, IG + 1).C(5, i); 
        end
    end
end

% -- sig21(B,G) = sig21(B+1,G) (Eq. 5.58)       
for IB = 1:NB - 1
    for IG = 1:NG
        row = row + 1;
        for i = 1:6
            NS = NG*(IB - 1) + IG;
            ic = ic + 1;
            NR(ic) = row;
            NC(ic)=6*(NS - 1) + i;
            AM(ic) = props(IB, IG).C(6, i);
            NS = NG*(IB) + IG;
            ic = ic + 1;
            NR(ic) = row;
            NC(ic)=6*(NS - 1) + i;
            AM(ic) = -props(IB + 1, IG).C(6, i);
        end
    end
end

% -- Define AM in matlab sparse format
AMsparse = sparse( NR(1:ic), NC(1:ic), AM(1:ic) );

% -- Assemble M and K (where R = K*epsb) (Eq. 5.62)
M = [AMsparse; AGsparse];
K = [zeros(5*NB*NG - 2*(NB + NG) - 1, 6); J]; 

% ----------------------------------------------------------------------
% -- Equations can be solved directly or by applying 6 strain components
%    See Chapter 5, Section 5.6
% ----------------------------------------------------------------------

% -- Direct Invert Solution (Faster)
% -- Solve for A
A = M\K;

% -- Extract 6x6 subcell strain concentration matrix from A
As(NB, NG).A = 0.;
for IG = 1:NG
    for IB = 1:NB
        astart=6*((IB - 1)*NG + IG - 1) + 1;
        aend = astart + 5;
        As(IB, IG).A = A(astart:aend, :);
    end
end

% ----------------------------------------------------------------------
% -- Apply 6 Components Solution (Slower)
%      for i = 1:6
%         
%         % -- Apply epsb(i) = 1, all other epsb = 0
%         dtemp = 0;
%         epsb = zeros(6,1);
%         epsb(i) = 1;
%         
%         R = K*epsb; % -- (see Eq. 5.62)
%         X=M\R;    % -- Subcell strains solution
%         
%         % -- Extract appropriate strain components from X
%         for IG = 1:NG
%             for IB = 1:NB
%                 xstart=6*((IB - 1)*NG + IG - 1) + 1;
%                 xend = xstart + 5;
%                 As(IB, IG).A(1:6, i) = X(xstart:xend);
%             end
%         end
%         
%      end
% ----------------------------------------------------------------------

% -- Calculate effective stiffness (Eq. 5.65)
Cstar = zeros(6);
for IG = 1:NG
    for IB = 1:NB
        Cstar = Cstar + h(IB)*l(IG)*props(IB,IG).C*As(IB, IG).A/(sumH*sumL);
    end
end

% ------------ THERMAL --------------
%  -- Effective compliance matrix
SS = inv(Cstar);

% -- Stress concentration matrix B
B(NB, NG).B = 0;
for IG = 1:NG
    for IB = 1:NB
        B(IB, IG).B = props(IB, IG).C*As(IB, IG).A*SS; % -- Eq. 5.66
        B(IB, IG).BT = transpose(B(IB, IG).B);
    end
end

% -- Effective CTEs (Eq. 5.68)
CTEs = zeros(6,1);
for IG = 1:NG
    for IB = 1:NB
        CTEs = CTEs + h(IB)*l(IG)*B(IB, IG).BT*props(IB, IG).alpha/(sumH*sumL);
    end
end

% -- Subcell strains
alpha_s = zeros(6*NB*NG,1);
for IG = 1:NG
    for IB = 1:NB
        astart=6*((IB - 1)*NG + IG - 1) + 1;
        alpha_s(astart) = props(IB, IG).alpha(1);
        alpha_s(astart+1) = props(IB, IG).alpha(2);
        alpha_s(astart+2) = props(IB, IG).alpha(3);
    end
end   

% -- Assemble RT (Eq. 5.74)
RT = [AMsparse; zeros(2*(NB+NG)+NB*NG+1, 6*NB*NG)]*alpha_s;

% -- Solve for subcell strains with epsb = 0 and DT = 1 (Eq. 5.73)
eps_s = M\RT;

% -- Extract 6x1 subcell thermal strain concentration matrix from eps_s
for IG = 1:NG
    for IB = 1:NB
        astart=6*((IB - 1)*NG + IG - 1) + 1;
        aend = astart + 5;
        As(IB, IG).AT = eps_s(astart:aend);
        % -- Store subcell stiffness and CTE
        As(IB, IG).C = props(IB, IG).C;
        As(IB, IG).alpha = props(IB, IG).alpha;
    end
end

% -- Get computational time for large RUCs
if NB*NG > 25*25
    ElapsedTime = toc;
    disp(['     GMC homogenization time = ', char(num2str(ElapsedTime)), ' seconds']);
end

end
