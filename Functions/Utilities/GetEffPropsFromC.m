function [Effprops] = GetEffPropsFromC(C)
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
% Purpose: Calculates the elastic properties from the stiffness
% - Assumes material is at most orthotropic (Eq. 2.19) 
% Input:
% - C: 6x6 material stiffness matrix
% Output:
% - Effprops: Struct containing material elastic properties
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

S = inv(C);

% -- Eqs. 2.34 - 2.36
Effprops.E1 = 1./S(1,1);
Effprops.E2 = 1./S(2,2);
Effprops.E3 = 1./S(3,3);

Effprops.G23 = 1./S(4,4);
Effprops.G13 = 1./S(5,5);
Effprops.G12 = 1./S(6,6);

Effprops.v12 = -S(2,1)/S(1,1);
Effprops.v13 = -S(3,1)/S(1,1);
Effprops.v23 = -S(3,2)/S(2,2);

end
