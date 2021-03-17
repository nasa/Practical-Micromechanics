function [Cstar,CTEs,As,Cstar_unavg,CTEs_unavg,As_unavg] = MOC(Fconstit, Mconstit, Vf)
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
% Purpose: Calculate the effective stiffness, effective CTES, and the mechanical and 
%          thermal strain concentration tensors for the Method of Cells (MOC) 
%          micromechanics theory (both with and without averaging for transverse
%          isotropy)
% Input:
% - Fconstit: Struct containing fiber material properties
% - Mconstit: Struct containing matrix material properties
% - Vf: Fiber volume fraction
% Output:
% - Cstar: Effective stiffness matrix (averaged)
% - CTEs: Effective CTEs (averaged)
% - As: Struct containing subcell concentration tensors (averaged)
% - Cstar_unavg: Effective stiffness matrix (unaveraged)
% - CTEs_unavg: Effective CTEs (unaveraged)
% - As_unavg: Struct containing subcell concentration tensors (unaveraged)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -- Get fiber stiffness matrix
Cf = GetCFromProps(Fconstit);

ALPHAf=[Fconstit.aL; Fconstit.aT; Fconstit.aT; 0; 0; 0];  % -- Fiber CTE vector

% -- Get matrix stiffness matrix
Cm = GetCFromProps(Mconstit);

ALPHAm=[Mconstit.aL; Mconstit.aL; Mconstit.aL; 0; 0; 0];  % -- Matrix CTE vector

% -- Get subcell dimensions from vf - square RUC
h1 = sqrt(Vf);
h2 = 1 - h1;
l1 = h1;
l2 = h2;

% -- Define matrix of coefficients for local subcell strains M(24,24) in Eq. 3.109   
M = MOCM(Cm, Cf, h1, h2, l1, l2);

% -- Initialize subcell strain concentration tensors
A11 = zeros(6);
A12 = zeros(6);
A21 = zeros(6);
A22 = zeros(6);

% -- Apply consecutively each far-field average strain component 
for I = 1: 6 
    
    % -- Vector of 6 global strain components (eps_bar)
    EPSB = zeros(6,1);
    EPSB(I) = 1;

    % -- Determine the right hand side vector R(24) in Eq. 3.109   
    R = MOCR(EPSB, h1, h2, l1, l2);

    % -- Solve Eq. 3.109 for the vector of subcell strains (stored as X)
    X = M\R; % -- Equivalent to X = inv(M)*R

    % -- Unaveraged strain concentration tensors in the four subcells
    A11(1, I) = X(1 );
    A11(2, I) = X(2 );
    A11(3, I) = X(3 );
    A11(4, I) = X(4 );
    A11(5, I) = X(5 );
    A11(6, I) = X(6 ); 

    A12(1, I) = X(7 );
    A12(2, I) = X(8 );
    A12(3, I) = X(9 );
    A12(4, I) = X(10);
    A12(5, I) = X(11);
    A12(6, I) = X(12);  

    A21(1, I) = X(13);
    A21(2, I) = X(14);
    A21(3, I) = X(15);
    A21(4, I) = X(16);
    A21(5, I) = X(17);
    A21(6, I) = X(18);  

    A22(1, I) = X(19);
    A22(2, I) = X(20);
    A22(3, I) = X(21);
    A22(4, I) = X(22);
    A22(5, I) = X(23);
    A22(6, I) = X(24);   

end  % -- End of loop over 6 applied global strain components

% -- Copy into struct
As_unavg(1,1).A = A11;
As_unavg(1,2).A = A12;
As_unavg(2,1).A = A21;
As_unavg(2,2).A = A22;

% -- Unaveraged effective stiffness matrix (Eq. 3.115)
Cstar_unavg = (h1*l1*Cf*A11 + h1*l2*Cm*A12 + h2*l1*Cm*A21 + h2*l2*Cm*A22)...
              /((h1+h2)*(l1+l2));
                                                         
% Average the stiffness terms (was an original method to get trans iso props)
%       E11 = Cstar_unavg(1,1);
%       E12 = (Cstar_unavg(1,2) + Cstar_unavg(1,3))/2;
%       E13 = E12;
%       E22 = 3*(Cstar_unavg(2,2) + Cstar_unavg(3,3))/8 + Cstar_unavg(2,3)/4 + Cstar_unavg(4,4)/2;
%       E33 = E22;
%       E23 = (Cstar_unavg(2,2) + Cstar_unavg(3,3))/8 + 3*Cstar_unavg(2,3)/4 - Cstar_unavg(4,4)/2;
%       E66 = Cstar_unavg(6,6);
%       E44 = (E22 - E23)/2;
                                                  
% -- Unaveraged effective compliance matrix
SS = inv(Cstar_unavg);

% -- Unaveraged stress concentration matrices (Eq. 3.118)
B11 = Cf*A11*SS;
B12 = Cm*A12*SS;
B21 = Cm*A21*SS;
B22 = Cm*A22*SS;

% -- Transpose of unaveraged stress concentration matrices
BT11 = transpose(B11);
BT12 = transpose(B12);
BT21 = transpose(B21);
BT22 = transpose(B22);

% -- Unaveraged CTEs (Eq. 3.165)
CTEs_unavg =(h1*l1*BT11*ALPHAf + h1*l2*BT12*ALPHAm + h2*l1*BT21*ALPHAm...
             + h2*l2*BT22*ALPHAm)/((h1 + h2)*(l1 + l2));

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -- Perform averaging to obtain transversely isotropy
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%  -- Averaged subcell strain concentration matrices (Eqs. 3.122)
COF1 = 1./2 + 1/pi;
COF2 = 1./2 - 1/pi;
COF3 = 3./8 + 1/pi;
COF4 = 3./8 - 1/pi;
COF5 = 1./2 + 1/pi;
COF6 = 1./2 - 1/pi;

AH11(1,1) = A11(1,1);
AH11(2,1) = COF1*A11(2,1) + COF2*A11(3,1);
AH11(3,1) = COF2*A11(2,1) + COF1*A11(3,1);
AH11(2,2) = COF3*A11(2,2) + COF4*A11(3,3) + (A11(2,3) + A11(3,2))/8 + A11(4,4)/4;
AH11(3,2) = COF3*A11(3,2) + COF4*A11(2,3) + (A11(2,2) + A11(3,3))/8 - A11(4,4)/4;
AH11(2,3) = COF3*A11(2,3) + COF4*A11(3,2) + (A11(2,2) + A11(3,3))/8 - A11(4,4)/4;
AH11(3,3) = COF3*A11(3,3) + COF4*A11(2,2) + (A11(2,3) + A11(3,2))/8 + A11(4,4)/4;
AH11(4,4) = (A11(2,2) + A11(3,3))/4 - (A11(2,3) + A11(3,2))/4 + A11(4,4)/2;
AH11(5,5) = COF5*A11(5,5) + COF6*A11(6,6);
AH11(6,6) = COF6*A11(5,5) + COF5*A11(6,6);

AH12(1,1) = A12(1,1);
AH12(2,1) = COF1*A12(2,1) + COF2*A12(3,1);
AH12(3,1) = COF2*A12(2,1) + COF1*A12(3,1);
AH12(2,2) = COF3*A12(2,2) + COF4*A12(3,3) + (A12(2,3) + A12(3,2))/8 + A12(4,4)/4;
AH12(3,2) = COF3*A12(3,2) + COF4*A12(2,3) + (A12(2,2) + A12(3,3))/8 - A12(4,4)/4;
AH12(2,3) = COF3*A12(2,3) + COF4*A12(3,2) + (A12(2,2) + A12(3,3))/8 - A12(4,4)/4;
AH12(3,3) = COF3*A12(3,3) + COF4*A12(2,2) + (A12(2,3) + A12(3,2))/8 + A12(4,4)/4;
AH12(4,4) = (A12(2,2) + A12(3,3))/4 - (A12(2,3) + A12(3,2))/4 + A12(4,4)/2;
AH12(5,5) = COF5*A12(5,5) + COF6*A12(6,6);
AH12(6,6) = COF6*A12(5,5) + COF5*A12(6,6);

AH21(1,1) = A21(1,1);
AH21(2,1) = COF1*A21(2,1) + COF2*A21(3,1);
AH21(3,1) = COF2*A21(2,1) + COF1*A21(3,1);
AH21(2,2) = COF3*A21(2,2) + COF4*A21(3,3) + (A21(2,3) + A21(3,2))/8 + A21(4,4)/4;
AH21(3,2) = COF3*A21(3,2) + COF4*A21(2,3) + (A21(2,2) + A21(3,3))/8 - A21(4,4)/4;
AH21(2,3) = COF3*A21(2,3) + COF4*A21(3,2) + (A21(2,2) + A21(3,3))/8 - A21(4,4)/4;
AH21(3,3) = COF3*A21(3,3) + COF4*A21(2,2) + (A21(2,3) + A21(3,2))/8 + A21(4,4)/4;
AH21(4,4) = (A21(2,2) + A21(3,3))/4 - (A21(2,3) + A21(3,2))/4 + A21(4,4)/2;
AH21(5,5) = COF5*A21(5,5) + COF6*A21(6,6);
AH21(6,6) = COF6*A21(5,5) + COF5*A21(6,6);

AH22(1,1) = A22(1,1);
AH22(2,1) = COF1*A22(2,1) + COF2*A22(3,1);
AH22(3,1) = COF2*A22(2,1) + COF1*A22(3,1);
AH22(2,2) = COF3*A22(2,2) + COF4*A22(3,3) + (A22(2,3) + A22(3,2))/8 + A22(4,4)/4;
AH22(3,2) = COF3*A22(3,2) + COF4*A22(2,3) + (A22(2,2) + A22(3,3))/8 - A22(4,4)/4;
AH22(2,3) = COF3*A22(2,3) + COF4*A22(3,2) + (A22(2,2) + A22(3,3))/8 - A22(4,4)/4;
AH22(3,3) = COF3*A22(3,3) + COF4*A22(2,2) + (A22(2,3) + A22(3,2))/8 + A22(4,4)/4;
AH22(4,4) = (A22(2,2) + A22(3,3))/4 - (A22(2,3) + A22(3,2))/4 + A22(4,4)/2;
AH22(5,5) = COF5*A22(5,5) + COF6*A22(6,6);
AH22(6,6) = COF6*A22(5,5) + COF5*A22(6,6);

% -- Copy into struct
As(1,1).A = AH11;
As(1,2).A = AH12;
As(2,1).A = AH21;
As(2,2).A = AH22;

% -- Averaged effective stiffness matrix (Eq. 3.124)
Cstar = (h1*l1*Cf*AH11 + h1*l2*Cm*AH12 + h2*l1*Cm*AH21 + h2*l2*Cm*AH22) ...
        /((h1 + h2)*(l1 + l2));

% -- Averaged effective compliance matrix
SS = inv(Cstar);

% -- Averaged stress concentration matrices (Eq. 3.125)
BH11 = Cf*AH11*SS;  
BH12 = Cm*AH12*SS;  
BH21 = Cm*AH21*SS;  
BH22 = Cm*AH22*SS;  

% -- Transpose of averaged stress concentration matrices
BHT11=transpose(BH11);
BHT12=transpose(BH12);
BHT21=transpose(BH21);
BHT22=transpose(BH22);

% -- Averaged CTEs (Eq. 3.165)
CTEs = (h1*l1*BHT11*ALPHAf + h1*l2*BHT12*ALPHAm + h2*l1*BHT21*ALPHAm + h2*l2...
        *BHT22*ALPHAm)/((h1 + h2)*(l1 + l2));

% -- Thermal strain concentration tensors

% -- Fiber and matrix thermal stress tensors (Eq. 3.128)
GAMf = Cf*ALPHAf;
GAMm = Cm*ALPHAm;

% -- R thermal (Eq. 3.172)
RT = zeros(24,1);
RT(9)  = GAMf(2) - GAMm(2);
RT(11) = GAMf(2) - GAMm(2);

% -- Subcell strains solution for EPSB = 0, DT = 1
X = M\RT;    

% -- Extract subcell thermal strain concentration tensors from X
AT11 = zeros(6,1);
AT11(1) = X(1 );
AT11(2) = X(2 );
AT11(3) = X(3 );
AT11(4) = X(4 );
AT11(5) = X(5 );
AT11(6) = X(6 ); 

AT12 = zeros(6,1);
AT12(1) = X(7 );
AT12(2) = X(8 );
AT12(3) = X(9 );
AT12(4) = X(10);
AT12(5) = X(11);
AT12(6) = X(12);  

AT21 = zeros(6,1);
AT21(1) = X(13);
AT21(2) = X(14);
AT21(3) = X(15);
AT21(4) = X(16);
AT21(5) = X(17);
AT21(6) = X(18);  

AT22 = zeros(6,1);
AT22(1) = X(19);
AT22(2) = X(20);
AT22(3) = X(21);
AT22(4) = X(22);
AT22(5) = X(23);
AT22(6) = X(24);   

% -- Copy into struct
As(1,1).AT = AT11;
As(1,2).AT = AT12;
As(2,1).AT = AT21;
As(2,2).AT = AT22;

% -- Copy into struct (unaverage AT = averaged AT)
As_unavg(1,1).AT = AT11;
As_unavg(1,2).AT = AT12;
As_unavg(2,1).AT = AT21;
As_unavg(2,2).AT = AT22;

end

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

function [M] = MOCM(Cm, Cf, h1, h2, l1, l2)

% -- Define matrix of coefficients for local subcell strains M(24,24) in Eq. 3.109   
    M = zeros(24,24);

% -- See Table 3.3 & Fig. 3.6 - Axial Strains (System Eqs 1-4)
    M(1 ,1 ) = 1; 
    M(2 ,7 ) = 1;  
    M(3 ,13) = 1;  
    M(4 ,19) = 1;  

% -- See Table 3.3 & Fig. 3.6 - Transverse Normal Strains (System Eqs 5-8)
    M(5 ,2 ) = h1;
    M(5 ,14) = h2;
    
    M(6 ,8 ) = h1;
    M(6 ,20) = h2;
    
    M(7 ,3 ) = l1;
    M(7 ,9 ) = l2;
    
    M(8 ,15) = l1;
    M(8 ,21) = l2;

% -- See Table 3.3 & Fig. 3.6 - Transverse Normal Stresses (System Eqs 9-12)
    M(9 ,1 ) = Cf(1,2);
    M(9 ,2 ) = Cf(2,2);
    M(9 ,3 ) = Cf(2,3);
    M(9 ,13) =-Cm(1,2);
    M(9 ,14) =-Cm(2,2);
    M(9 ,15) =-Cm(2,3);

    M(10,7 ) = Cm(1,2);
    M(10,8 ) = Cm(2,2);
    M(10,9 ) = Cm(2,3);
    M(10,19) =-Cm(1,2);
    M(10,20) =-Cm(2,2);
    M(10,21) =-Cm(2,3);

    M(11,1 ) = Cf(1,3);
    M(11,2 ) = Cf(2,3);
    M(11,3 ) = Cf(3,3);
    M(11,7 ) =-Cm(1,3);
    M(11,8 ) =-Cm(2,3);
    M(11,9 ) =-Cm(3,3);

    M(12,13) = Cm(1,3);
    M(12,14) = Cm(2,3);
    M(12,15) = Cm(3,3);
    M(12,19) =-Cm(1,3);
    M(12,20) =-Cm(2,3);
    M(12,21) =-Cm(3,3);

% -- See Table 3.3 & Fig. 3.6 - 23 Shear Strains (System Eq 13)
    M(13,4 ) = h1*l1;
    M(13,10) = h1*l2;
    M(13,16) = h2*l1;
    M(13,22) = h2*l2;

% -- See Table 3.3 & Fig. 3.6 - 23 Shear Stresses (System Eq 14-16)
    M(14,4 ) = 2*Cf(4,4);
    M(14,10) =-2*Cm(4,4);   

    M(15,10) = 2*Cm(4,4);   
    M(15,16) =-2*Cm(4,4);   

    M(16,16) = 2*Cm(4,4);   
    M(16,22) =-2*Cm(4,4);

% -- See Table 3.3 & Fig. 3.6 - 12 Shear Strains (System Eq 17-18)
    M(17,6 ) = h1;
    M(17,18) = h2;

    M(18,12) = h1;
    M(18,24) = h2;

% -- See Table 3.3 & Fig. 3.6 - 12 Shear Stresses (System Eq 19-20)
    M(19,6 ) = 2*Cf(6,6);
    M(19,18) =-2*Cm(6,6);   

    M(20,12) = 2*Cm(6,6);   
    M(20,24) =-2*Cm(6,6);

% -- See Table 3.3 & Fig. 3.6 - 13 Shear Strains (System Eq 21-22)
    M(21,5 ) = l1;
    M(21,11) = l2;

    M(22,17) = l1;
    M(22,23) = l2;

% -- See Table 3.3 & Fig. 3.6 - 13 Shear Stresses (System Eq 23-24)
    M(23,5 ) = 2*Cf(5,5);
    M(23,11) =-2*Cm(5,5);   

    M(24,17) = 2*Cm(5,5);   
    M(24,23) =-2*Cm(5,5);  

end

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

function [R] = MOCR(EPSB, h1, h2, l1, l2)

% -- Determine the right hand side vector R(24) in Eq. 3.109   

% -- See Table 3.3 & Fig. 3.6
    R(1 ,1) = EPSB(1); % -- Axial Strains (System Eqs 1-4)
    R(2 ,1) = EPSB(1);
    R(3 ,1) = EPSB(1);
    R(4 ,1) = EPSB(1);
    
    R(5 ,1) = (h1 + h2)*EPSB(2); % -- Transverse Normal Strains (System Eqs 5-8) 
    R(6 ,1) = (h1 + h2)*EPSB(2); 
    R(7 ,1) = (l1 + l2)*EPSB(3); 
    R(8 ,1) = (l1 + l2)*EPSB(3);
    
    R(9 ,1) = 0; % -- Transverse Normal Stresses (System Eqs 9-12)
    R(10,1) = 0;
    R(11,1) = 0;
    R(12,1) = 0;

    R(13,1) = (h1 + h2)*(l1 + l2)*EPSB(4); % -- 23 Shear Strains (System Eq 13)
    
    R(14,1) = 0; % -- 23 Shear Stresses (System Eq 14-16)
    R(15,1) = 0;
    R(16,1) = 0;

    R(17,1) = (h1 + h2)*EPSB(6); % -- 12 Shear Strains (System Eq 17-18)
    R(18,1) = (h1 + h2)*EPSB(6); 
    
    R(19,1) = 0; % -- 12 Shear Stresses (System Eq 19-20)
    R(20,1) = 0;

    R(21,1) = (l1 + l2)*EPSB(5); % -- 13 Shear Strains (System Eq 21-22) 
    R(22,1) = (l1 + l2)*EPSB(5);  
    
    R(23,1) = 0; % -- 13 Shear Stresses (System Eq 23-24)
    R(24,1) = 0; 

end

