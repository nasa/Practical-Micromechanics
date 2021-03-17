function [Cstar, CTEs, Af_MT, Am_MT, ATf_MT, ATm_MT] = MT(Fconstit, Mconstit, Vf)
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
%          thermal strain concentration tensors for the Mori-Tanaka (MT)
%          micromechanics theory
% Input:
% - Fconstit: Struct containing fiber material properties
% - Mconstit: Struct containing matrix material properties
% - Vf: Fiber volume fraction
% Output:
% - Cstar: Effective stiffness matrix
% - CTEs: Effective CTEs
% - Af_MT: MT fiber strain concentration tensor
% - Af_MT: MT matrix strain concentration tensor
% - Af_MT: MT fiber thermal strain concentration tensor
% - Af_MT: MT matrix therrmal strain concentration tensor
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -- Prevent divide by zero for Vf = 1
if Vf > 0.99999999
    Vf = 0.99999999;
end

% -- Get fiber stiffness matrix
Cf = GetCFromProps(Fconstit);

ALPHAf=[Fconstit.aL; Fconstit.aT; Fconstit.aT; 0; 0; 0];  % -- Fiber CTE vector

% -- Get matrix stiffness matrix
Cm= GetCFromProps(Mconstit);
Num = Mconstit.vL;

ALPHAm=[Mconstit.aL; Mconstit.aL; Mconstit.aL; 0; 0; 0];  % -- Matrix CTE vector

% -- Eshelby Tensor for cylindrical fiber (Eq. 3.43)
[P] = EshelbyCyl(Num);

% -- Eshelby Tensor for spherical inclusion (Eq. 3.42)
%[P] = EshelbySphere(Num);

C = (Cm^-1)*(Cm - Cf); % -- Eq. 3.41

% -- Dilute fiber strain concentration tensor (Eq. 3.40)
Af_dilute = inv(eye(6) - P*C);

% -- MT fiber strain concentration tensor (Eq. 3.57)
Af_MT = Af_dilute*inv((1-Vf)*eye(6) + Vf*Af_dilute);

% -- MT effective stiffness matrix (Eq. 3.46)
Cstar = Cm + Vf*(Cf - Cm)*Af_MT;

% -- MT matrix strain concentration tensor (Eq. 3.44)
Am_MT = (eye(6)- Vf*Af_MT)/(1 - Vf); 

% ------------ THERMAL --------------
%  -- Effective compliance matrix
SS=inv(Cstar);

% -- Stress concentration tensors (Eq. 3.24)
BHf=Cf*Af_MT*SS;  
BHm=Cm*Am_MT*SS;  

% -- Transposes of stress concentration tensors
BHTf=transpose(BHf);
BHTm=transpose(BHm);

% -- Effective CTEs (Eq. 3.138)
CTEs = Vf*BHTf*ALPHAf+(1 - Vf)*BHTm*ALPHAm;

% -- MT thermal strain concentration tensors (Eqs. 3.147 & 3.148)
ATf_MT = (eye(6) - Af_MT)*inv(Cf - Cm)*(Cf*ALPHAf - Cm*ALPHAm);
ATm_MT = (eye(6) - Am_MT)*inv(Cf - Cm)*(Cf*ALPHAf - Cm*ALPHAm);

end

%--------------------------------------------------------------
%--------------------------------------------------------------

function [P] = EshelbyCyl(Num)
    % -- Eq. 3.43
    P(2,1) = Num/(2.*(1. - Num));
    P(3,1) = Num/(2.*(1. - Num));
    P(2,2) = (1./(2.*(1. - Num))) * (3./4. + (1. - 2.*Num)/2.);    
    P(3,3) = P(2,2);
    P(2,3) = (1./(2.*(1. - Num))) * (1./4. - (1. - 2.*Num)/2.);
    P(3,2) = P(2,3);
    % -- Engineering shear strains - multiply tensorial P by 2
    P(4,4) = (2./(2.*(1. - Num))) * (1./4. + (1. - 2.*Num)/2.);
    P(5,5) = 2./4.;
    P(6,6) = 2./4.;
end

%--------------------------------------------------------------
%--------------------------------------------------------------
