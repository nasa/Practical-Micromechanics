function [T1,T2] = GetT(A)
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
% Purpose: Calculate the two transformation matrices T1 and T2
% Input:
% - A: Ply orientation angle (deg) (theta - see Eqs. 2.44 and 2.45)
% Output:
% - T1 stress transformation matrix (Eq. 2.44)
% - T2 strain transformation matrix (Eq. 2.45)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

m = cosd(A);
n = sind(A);

% -- T1 is transformation matrix for stress that maps from material to
%    global coordinate system (Eq. 2.44)
T1 = [m.^2, n.^2, 2*n.*m; ...
      n.^2, m.^2, -2*m.*n; ...
      -m.*n, m.*n, (m.^2)-(n.^2)];
  
% -- T2 is transformation matrix for strains that maps from material to
%    global coordinate system (Eq. 2.45)
T2 = [m.^2, n.^2, n.*m; ...
      n.^2, m.^2, -m.*n; ...
      -2*m.*n, 2*m.*n, (m.^2)-(n.^2)];

end