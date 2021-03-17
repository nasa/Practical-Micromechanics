function [Aij, Bij, Dij] = GetABD(Orient, plyprops, plymat, N, z)
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
% Purpose: Calculate the laminate ABD matrix
% Input:
% - Orient: Array of ply orientation angles theta (deg)
% - plyprops: Cell array containing ply material properties
% - plymat: Array of ply material id numbers
% - N: Total number of plies
% - z: Array of ply boundary z-coordinates
% Output:
% - Aij: Laminate extensional stiffness matrix
% - Bij: Laminate coupling stiffness matrix
% - Dij: Laminate bending stiffness matrix
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -- Preaallocate arrays
Aij = zeros(3,3);
Bij = zeros(3,3);
Dij = zeros(3,3);

Aij_vect = zeros(1,N);
Bij_vect = zeros(1,N);
Dij_vect = zeros(1,N);

% -- Calculate for A,B,D vectors by using for-loops to create summation
for i = 1:3
    for j = 1:3
        for k = 1:N  % -- k is current ply
            
            Qb = Qbar(Orient(k), plyprops, plymat(k));

            % -- .* and .^ are element by element operations
            Aij_vect(k) = Qb(i,j).*(z(k+1) - z(k));           % -- Eq. 2.64
            Bij_vect(k) = Qb(i,j).*((z(k+1)).^2 - (z(k)).^2); % -- Eq. 2.65
            Dij_vect(k) = Qb(i,j).*((z(k+1)).^3 - (z(k)).^3); % -- Eq. 2.66

        end
        
        Aij(i,j) = sum(Aij_vect);   % -- Eq. 2.64
        Bij(i,j) = sum(Bij_vect)/2; % -- Eq. 2.65
        Dij(i,j) = sum(Dij_vect)/3; % -- Eq. 2.66
        
    end   
end

end
