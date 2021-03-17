function [C] = GetCFromProps(Props)
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
% Purpose: Calculates the stiffness matrix from the elastic properties
% - Assumes material is at most transversely isotropic (Eq. 2.16) 
% Input:
% - Props: Struct containing material elastic properties
% Output:
% - C: 6x6 material stiffness matrix
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

S = zeros(6);

% -- Eq. 2.17
S(1,1) = 1./Props.EL;
S(2,2) = 1./Props.ET;
S(3,3) = 1./Props.ET;

% -- Eq. 2.18
S(4,4) = 1./(Props.ET/(2.*(1. + Props.vT)));
S(5,5) = 1./Props.GL;
S(6,6) = 1./Props.GL;

% -- Eq. 2.17
S(2,1)= -Props.vL/Props.EL;
S(1,2)=S(2,1);
S(3,1)= -Props.vL/Props.EL;
S(1,3)=S(3,1);
S(3,2)= -Props.vT/Props.ET;
S(2,3)= S(3,2);

C = inv(S);

end


      

         