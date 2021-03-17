function [NMT, alphbar] = GetThermLoadsNM(Orient, plyprops, plymat, N, z)
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
% Purpose: Calculate the laminate thermal force and moment results, NT and MT **per 
%          degree** - That is, these must be multiplied by DT to obtain actual force 
%          and moment resultants in Eqs. 2.67 and 2.68
% Input:
% - Orient: Array of ply orientation angles theta (deg)
% - plyprops: Cell array containing ply material properties
% - plymat: Array of ply material id numbers
% - N: Total number of plies
% - z: Array of ply boundary z-coordinates
% Output:
% - NMT: Vector containing NT and MT (see Eq. 2.70) **per degree**
% - alphbar: Cell array containing ply CTEs in the global (x-y) coordinates
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -- Preallocate alphbar
alphbar = cell(1,N);

% -- Obtain alphbar - cell array of ply CTEs in x-y coords
for k = 1:N
    a1 = plyprops{plymat(k)}.a1;   
    a2 = plyprops{plymat(k)}.a2;
    alph =[a1; a2; 0]; % -- alph is the ply CTE vector, a12 = 0 in local 
    [~,T2] = GetT(Orient(k));
    alphbar{k} = inv(T2) * alph;  % -- Eq. 2.48                     
end

% -- Preallocate NT and MT
NT = zeros(3,1);
MT = zeros(3,1);

% -- Perform summation over plies
for k = 1:N 
    Qb = Qbar(Orient(k), plyprops, plymat(k));   
    NTply = Qb*alphbar{k}*(z(k+1) - z(k)); % -- Eq. 2.67 (divided by DT)
    NT = NT + NTply;
    MTply = Qb*alphbar{k}*(z(k+1)^2 - z(k)^2)/2; % -- Eq. 2.68 (divided by DT)
    MT = MT + MTply;
end

% -- Store NT and MT in single vector (NOTE: these are per degree)
NMT = [NT;MT]; % [Nx,Ny,Nxy,Mx,My,Mxy]

end