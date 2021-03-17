function [plyprops, LamResults] = DamageLam(plyprops, Geometry, LamResults, OutInfo)
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
%
% Purpose: For laminate, call DamageMicro.m per ply to replace subcells or materials 
%          with lowest negative MoS with damaged materials (i.e., reduced stiffnesses)
%          and run micromechanics analysis to get updated effective properties and 
%          concentration tensors
% Input:
% - plyprops: Cell/struct containing effective composite material properties 
% - Geometry: Struct containing laminate definition variables
% - LamResults: Struct containing laminate analysis results
% - OutInfo: Struct containing output information
% Output:
% - plyprops: Updated cell/struct containing effective composite material properties
% - LamResults: Updated struct containing laminate analysis results
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -- Return if no Margins have been calculated
if ~LamResults.MicroMoS
    return;
end

% -- Loop through top and bottom points in each ply
for kk = 1:2*Geometry.N
    k = round(kk/2);
    mat = Geometry.plymat(k);

    % -- No bending for damage case, so can assume both points per ply are identical
    %    so only run micromechanics only for 1st point in each ply, then copy to 2nd point
    if mod(kk,2) == 1
        [plyprops{mat}, LamResults.Micro(kk).MoS] = DamageMicro(plyprops{mat}, LamResults.Micro(kk).MoS, OutInfo, k);
    else
        LamResults.Micro(kk).MoS = LamResults.Micro(kk - 1).MoS;
    end
end

end
