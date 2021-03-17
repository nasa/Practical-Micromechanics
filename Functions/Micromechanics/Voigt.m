function [Cstar, CTEs, Af_V, Am_V, ATf_V, ATm_V] = Voigt(Fconstit, Mconstit, Vf)
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
% Purpose: Calculate the mechanical and thermal strain concentration tensors for the
%          Voigt micomechanics theory
% Input:
% - Fconstit: Struct containing fiber material properties
% - Mconstit: Struct containing matrix material properties
% - Vf: Fiber volume fraction
% Output:
% - Cstar: Effective stiffness matrix
% - CTEs: Effective CTEs
% - Af_V: Voigt fiber strain concentration tensor
% - Am_V: Voigt matrix strain concentration tensor
% - ATf_V: Voigt fiber thermal strain concentration tensor
% - ATm_V: Voigt matrix thermal strain concentration tensor
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -- Get fiber stiffness matrix
Cf = GetCFromProps(Fconstit);

ALPHAf = [Fconstit.aL; Fconstit.aT; Fconstit.aT; 0; 0; 0];  % -- fiber CTE vector

% -- Get matrix stiffness matrix
Cm = GetCFromProps(Mconstit);

ALPHAm = [Mconstit.aL; Mconstit.aL; Mconstit.aL; 0; 0; 0];  % -- matrix CTE vector

% -- Voigt strain concentration tensors (Eq. 3.25)
Af_V = eye(6);
Am_V = eye(6);

% -- Voigt effective stiffness matrix (Eq. 3.26)
Cstar = Vf*Cf*Af_V + (1 - Vf)*Cm*Am_V;


% ------------ THERMAL --------------

% -- Voigt thermal strain concentration tensors (Eq. 3.153)
ATf_V = zeros(6,1);
ATm_V = zeros(6,1);

% -- Effective compliance matrix
SS=inv(Cstar);

% -- Stress concentration tensors (Eq. 3.24)
BHf=Cf*Af_V*SS;  
BHm=Cm*Am_V*SS;  

% -- Transposes of stress concentration tensors
BHTf=transpose(BHf);
BHTm=transpose(BHm);

% -- Effective CTEs (Eq. 3.138)
CTEs = Vf*BHTf*ALPHAf+(1 - Vf)*BHTm*ALPHAm;

end         