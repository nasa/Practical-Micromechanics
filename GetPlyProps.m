function [plyprops] = GetPlyProps(Geometry)
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
% Purpose: Provides ply level material properties to Lamination Theory
% - Ply Properties can be defined as transversely isotropic directly or calculated 
%   via micromechanics
% Input:
% - Geometry: Cell array containing laminate definition variables
% Output:
% - plyprops: Cell array containing ply material properties and allowables
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

NmatMax = 300;
plyprops = cell(1,NmatMax);

% -- Determine which mats are used in any problem
for NP = 1: length(Geometry)
    if isfield(Geometry{NP},'plymat')
        for I = 1: length(Geometry{NP}.plymat)
            mat = Geometry{NP}.plymat(I);
            plyprops{mat}.used = true;
        end
    end
end

% -- Monolithic Transversely Isotropic Materials
plyprops{1}.name = "IM7/8552";
plyprops{1}.E1 = 146.8e3;    
plyprops{1}.E2 = 8.69e3;    
plyprops{1}.G12 = 5.16e3;   
plyprops{1}.v12 = 0.32;
plyprops{1}.a1 = -0.1e-6;  
plyprops{1}.a2 = 31.e-6;   
% -- Material allowables introduced in Chapter 4
plyprops{1}.allowables.XT = 2323;
plyprops{1}.allowables.XC = -1531;
plyprops{1}.allowables.YT = 52.3;
plyprops{1}.allowables.YC = -235;
plyprops{1}.allowables.S = 88;
plyprops{1}.allowables.XeT = plyprops{1}.allowables.XT/plyprops{1}.E1;
plyprops{1}.allowables.XeC = plyprops{1}.allowables.XC/plyprops{1}.E1;
plyprops{1}.allowables.YeT = plyprops{1}.allowables.YT/plyprops{1}.E2;
plyprops{1}.allowables.YeC = plyprops{1}.allowables.YC/plyprops{1}.E2;
plyprops{1}.allowables.Se = plyprops{1}.allowables.S/plyprops{1}.G12;

plyprops{2}.name = "Glass/Epoxy";
plyprops{2}.E1 = 43.5e3;
plyprops{2}.E2 = 11.5e3;
plyprops{2}.G12 = 3.45e3;
plyprops{2}.v12 = 0.27;
plyprops{2}.a1 = 6.84e-06;
plyprops{2}.a2 = 29.e-06;
% -- Material allowables introduced in Chapter 4
plyprops{2}.allowables.XT = 1346;
plyprops{2}.allowables.XC = -944;
plyprops{2}.allowables.YT = 50;
plyprops{2}.allowables.YC = -245;
plyprops{2}.allowables.S = 122;
plyprops{2}.allowables.XeT = plyprops{2}.allowables.XT/plyprops{2}.E1;
plyprops{2}.allowables.XeC = plyprops{2}.allowables.XC/plyprops{2}.E1;
plyprops{2}.allowables.YeT = plyprops{2}.allowables.YT/plyprops{2}.E2;
plyprops{2}.allowables.YeC = plyprops{2}.allowables.YC/plyprops{2}.E2;
plyprops{2}.allowables.Se = plyprops{2}.allowables.S/plyprops{2}.G12;

plyprops{3}.name = "SCS-6/Ti-15-3";
plyprops{3}.E1 = 221.e3;
plyprops{3}.E2 = 145.e3;
plyprops{3}.G12 = 53.2e3;
plyprops{3}.v12 = 0.27;
plyprops{3}.a1 = 6.15e-06;
plyprops{3}.a2 = 7.9e-06;

plyprops{4}.name = "SiC/SiC";
plyprops{4}.E1 = 300.e3;
plyprops{4}.E2 = 152.4e3;
plyprops{4}.G12 = 69.5e3;
plyprops{4}.v12 = 0.17;
plyprops{4}.a1 = 4.88e-06;
plyprops{4}.a2 = 4.91e-06;

plyprops{5}.name = "Aluminum";
plyprops{5}.E1 = 70.e3;
plyprops{5}.E2 = 70.e3;
plyprops{5}.G12 = 26.923e3;
plyprops{5}.v12 = 0.3;
plyprops{5}.a1 = 23.e-06;
plyprops{5}.a2 = 23.e-06;

%************************************************
% -- Calculate Ply properties from micromechanics
%************************************************

%=========================================================================
% -- Define Constituent Properties (stored in constitprops)
%=========================================================================
[constitprops] = GetConstitProps;

% Check constituent compressive allowables -- must be negative
for Mat = 1: length(constitprops)
    if isfield(constitprops{Mat},'allowables')
        constitprops{Mat}.allowables.XC = -abs(constitprops{Mat}.allowables.XC);
        constitprops{Mat}.allowables.YC = -abs(constitprops{Mat}.allowables.YC);
        constitprops{Mat}.allowables.ZC = -abs(constitprops{Mat}.allowables.ZC);
    end
end

%=========================================================================
% -- Composite materials with micromechanics
%=========================================================================
[plyprops] = GetEffProps(constitprops, plyprops);


% Check ply material compressive allowables -- must be negative
for Mat = 1: length(plyprops)
    if isfield(plyprops{Mat},'allowables')
        plyprops{Mat}.allowables.XC = -abs(plyprops{Mat}.allowables.XC);
        plyprops{Mat}.allowables.YC = -abs(plyprops{Mat}.allowables.YC);
    end
end


end

