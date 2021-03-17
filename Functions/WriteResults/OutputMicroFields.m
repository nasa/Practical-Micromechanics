function [plyprops] = OutputMicroFields(OutInfo, Geometry, plyprops, LamResults)
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
% Purpose: For micromechanics-base CLT plies, call PlotMicroFields.m to make
%          microscale plots of the stress and strain fields
% Input:
% - OutInfo: Struct containing output information
% - Geometry: Struct containing laminate definition variables
% - plyprops: Cell/struct containing effective composite material properties 
% - LamResults: Struct containing laminate analysis results
% Output:
% - plyprops: Updated cell/struct containing effective composite material 
%             properties 
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

N = Geometry.N;
T = true;
F = false;
MakeMicroPlots = OutInfo.MakePlots;
Kappa = LamResults.EK(4:6);

% -- Turn on and off plotting of specific through-thickness point here
if MakeMicroPlots
    PlotPlyMicro(1:2*N) = T;
else
    PlotPlyMicro(1:2*N) = F;
end

% -- Check for bending
if (abs(sum(Kappa)) > 1.e-10)
    bending = T;
else
    bending = F;
end

% -- For cases without bending, turn off plots for even ply points so only 
%    one plot per ply
if (~bending)
    for kk = 2:2:2*N
       PlotPlyMicro(kk) = F;
    end
end

% -- Plot microfields per top/bottom points per ply
for kk = 1:2*Geometry.N
    
    k = round(kk/2);
    mat = Geometry.plymat(k);

    if isfield(plyprops{mat},'micro') % -- Check if ply is micromech-based

        if PlotPlyMicro(kk) % -- Check if plots should be made for this pt
            PlotMicroFields(OutInfo, plyprops{mat}, LamResults.Micro(kk), ...
                            Geometry.Orient(k), k, kk, bending, ...
                            Geometry.LayerZ(kk));
        end

    end
            
end

end