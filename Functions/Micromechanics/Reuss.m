function [Cstar, CTEs, Af_R, Am_R, ATf_R, ATm_R] = Reuss(Fconstit,Mconstit,Vf)
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
%          thermal strain concentration tensors for the Reuss micomechanics theory
% Input:
% - Fconstit: Struct containing fiber material properties
% - Mconstit: Struct containing matrix material properties
% - Vf: Fiber volume fraction
% Output:
% - Cstar: Effective stiffness matrix
% - CTEs: Effective CTEs
% - Af_R: Reuss fiber strain concentration tensor
% - Am_R: Reuss matrix strain concentration tensor
% - ATf_R: Reuss fiber thermal strain concentration tensor
% - ATm_R: Reuss matrix thermal strain concentration tensor
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -- Get fiber stiffness and compliance matrices
Cf = GetCFromProps(Fconstit);
Sf = inv(Cf);

ALPHAf = [Fconstit.aL; Fconstit.aT; Fconstit.aT; 0; 0; 0];  % -- fiber CTE vector

% -- Get matrix stiffness matrix
Cm= GetCFromProps(Mconstit);
Sm = inv(Cm);

ALPHAm = [Mconstit.aL; Mconstit.aL; Mconstit.aL; 0; 0; 0];  % -- matrix CTE vector
         
% -- Reuss stress concentration tensors (Eq. 3.27)
Bf_R = eye(6);
Bm_R = eye(6);

% -- Reuss effective compliance matrix (Eq. 3.28)
Sstar = Vf*Sf*Bf_R + (1 - Vf)*Sm*Bm_R; 

% -- Effective stiffness matrix
Cstar = inv(Sstar);

% -- Reuss strain concentration tensors (Eq. 3.22)
Af_R = inv(Cf)*Bf_R*Cstar;
Am_R = inv(Cm)*Bm_R*Cstar;

% ------------ THERMAL --------------

% -- Transposes of stress concentration tensors
BHTf=transpose(Bf_R);
BHTm=transpose(Bm_R);

% -- Effective CTEs (Eq. 3.138)
CTEs = Vf*BHTf*ALPHAf+(1 - Vf)*BHTm*ALPHAm;

% -- Reuss thermal strain concentration tensors (Eqs. 3.147 & 3.148)
ATf_R = (eye(6) - Af_R)*inv(Cf - Cm)*(Cf*ALPHAf - Cm*ALPHAm);
ATm_R = (eye(6) - Am_R)*inv(Cf - Cm)*(Cf*ALPHAf - Cm*ALPHAm);

end      