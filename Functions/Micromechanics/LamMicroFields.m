function [LamResults] = LamMicroFields(Geometry, Loads, plyprops, LamResults)
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
% Purpose: Calculate the local (constituent level) stress and strain fields at the 
%          top and bottom of each ply in a composite laminate by first determining 
%          the full 3D strain at each point and then calling MicroFields.m
% Input:
% - Geometry: Struct containing laminate definition variables
% - Loads: Struct containing problem loading information
% - plyprops: Cell/struct containing effective composite material properties 
% - LamResults: Struct containing laminate analysis results
% Output:
% - LamResults: Struct containing laminate analysis results updated to include local 
%               micromechanics result ("Micro(kk)") per ply location
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

DT = Loads.DT;

% -- Calculate and plot microfields per ply top and bottom location (kk)
for kk = 1: 2*Geometry.N
    
    k = round(kk/2); % -- Ply number
    mat = Geometry.plymat(k);
        
    if isfield(plyprops{mat},'micro')
        
        % -- Calculate eps_zz = eps_33 (Eq. 2.76)
        Cstar = plyprops{mat}.Cstar; 
        alpha = [plyprops{mat}.a1; plyprops{mat}.a2; plyprops{mat}.a3];
        Ptstrain(:) = LamResults.MCstrain(:, kk);
        epsz = -(Ptstrain(1) - alpha(1)*DT)*Cstar(1,3)/Cstar(3,3) ...
               -(Ptstrain(2) - alpha(2)*DT)*Cstar(2,3)/Cstar(3,3) ...
               + alpha(3)*DT;

        FullGlobalStrain = [Ptstrain(1); Ptstrain(2); epsz; 0; 0; Ptstrain(3)];

        % -- Call MicroFields.m for this ply location
        if ~isfield(LamResults, 'Micro') 
            [LamResults.Micro(kk)] = MicroFields(FullGlobalStrain, DT, plyprops{mat});
        elseif length(LamResults.Micro) < kk
            [LamResults.Micro(kk)] = MicroFields(FullGlobalStrain, DT, plyprops{mat});
        else
            [LamResults.Micro(kk)] = MicroFields(FullGlobalStrain, DT, ...
                                                 plyprops{mat}, LamResults.Micro(kk));
        end
        
    end
        
end

end