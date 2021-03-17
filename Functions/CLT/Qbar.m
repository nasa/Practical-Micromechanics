function Q_bar = Qbar(A, plyprops, matID)
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
% Purpose: Calculate the transformed reduced stiffness matrix
% Input:
% - A: Ply orientation angle theta (deg)
% - plyprops: Cell array containing ply material properties
% - matID: ply-level material id number
% Output:
% - Q_bar: transformed reduced stiffness matrix
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

E1 = plyprops{matID}.E1;
E2 = plyprops{matID}.E2;
G12 = plyprops{matID}.G12;
v12 = plyprops{matID}.v12;

% -- Define plane-stress compliance matrix (Eqs. 2.17, 2.18)
S11 = 1/E1;
S12 = -v12/E1;
S22 = 1/E2;
S66 = 1/G12;

S = [S11 S12 0;S12 S22 0;0 0 S66];

% -- Define the reduced stiffness matrix in material coordinates (1,2,3)
Q = inv(S); % -- Equivalent to Eqs. 2.42

% -- Determine transformation matrices
[T1,T2] = GetT(A);

% -- Calculate transformed reduced stiffness matrix (Eq. 2.46)
Q_bar = inv(T1)*Q*T2;

end