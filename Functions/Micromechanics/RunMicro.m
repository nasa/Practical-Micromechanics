function [plyprops] = RunMicro(Theory, name, Vol, Constits, plyprops, RUCid)
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
% Purpose: Runs the chosen micromechanics theory to obtain the effective properties 
%          and strain concentration tensors
% Input:
% - Theory: Chosen micromechanics theory
% - name: Micromechanics-based material name
% - Vol: Contains either Vf or is struct containing Vf and Vi
% - Constits: Struct containing constituent material properties
% - plyprops: Struct containing effective properties
% - RUCid: RUC i.d. number or full RUC struct containing RUC info
%          Optional argument, only used for GMC and HFGMC (Ch. 5 & 6)
% Output:
% - plyprops: Updated struct containing effective properties
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -- Do not run micromechanics if this material not used in any problems
if ~isfield(plyprops, 'used')
    plyprops.used = false;
    return
elseif ~plyprops.used
    return
end

% -- Extract fiber and interface vol fractions (Vol can be Vf or struct)
if isstruct(Vol)
    Vf = Vol.Vf;
    if isfield(Vol, 'Vi')
        Vi = Vol.Vi;
    else
        Vi = 0;
    end
else
    Vf = Vol;
    Vi = 0;
end

% -- Run the specified micromechanics theory
switch(Theory)
    case 'Voigt'
        [Cstar, CTEs, Af, Am, ATf, ATm] = Voigt(Constits.Fiber, Constits.Matrix, Vf);
        plyprops.micro = "Voigt";
  
    case 'Reuss'
        [Cstar, CTEs, Af, Am, ATf, ATm] = Reuss(Constits.Fiber, Constits.Matrix, Vf);
        plyprops.micro = "Reuss";
        
    case 'MT'
        [Cstar, CTEs, Af, Am, ATf, ATm] = MT(Constits.Fiber, Constits.Matrix, Vf);
         plyprops.micro = "MT";
         
    case 'MOC'
        [Cstar,CTEs,As,Cstaru,CTEsu,Asu] = MOC(Constits.Fiber, Constits.Matrix, Vf);  
        plyprops.micro = "MOC";
        [RUC] = GetRUC(Constits, Vf, Vi, 2);

    case 'MOCu'
        [Cstar,CTEs,As,Cstaru,CTEsu,Asu] = MOC(Constits.Fiber, Constits.Matrix, Vf);  
        Cstar = Cstaru;
        CTEs = CTEsu;
        As = Asu;
        plyprops.micro = "MOCu";
        
    case 'GMC'
         if ~exist('RUCid','var')
             error('*** No RUCid specified for GMC ***');
         end
        plyprops.micro = "GMC";
        if isstruct(RUCid)
            RUC = RUCid;
        else
            [RUC] = GetRUC(Constits, Vf, Vi, RUCid);
        end
        if Vf <= RUC.Vf_max 
            [Cstar,CTEs,As] = GMC(Constits, RUC); 
        else
            Cstar = eye(6,6);
            CTEs = zeros(6,1);
            As = 0;
        end
            
    case 'HFGMC'
         if ~exist('RUCid','var')
             error('*** No RUCid specified for HFGMC ***');
         end
        plyprops.micro = "HFGMC";
        if isstruct(RUCid)
            RUC = RUCid;
        else
            [RUC] = GetRUC(Constits, Vf, Vi, RUCid);
        end
        if Vf <= RUC.Vf_max 
            [Cstar,CTEs,As] = HFGMC(Constits,RUC);        
        else
            Cstar = eye(6,6);
            CTEs = zeros(6,1);
            As = 0;
        end
        
    otherwise
        error(['Invalid Theory Specified in RunMicro - Theory = ', Theory])
        
end

% -- Store info in plyprops struct
plyprops.name = name;
plyprops.Vf = Vf;
plyprops.Vi = Vi;
plyprops.Fiber = Constits.Fiber;
plyprops.Matrix = Constits.Matrix;
if isfield(Constits,'Interface')
    plyprops.Interface = Constits.Interface;
end

[Effprops] = GetEffPropsFromC(Cstar);
plyprops.E1 = Effprops.E1;
plyprops.E2 = Effprops.E2;
plyprops.E3 = Effprops.E3;
plyprops.v12 = Effprops.v12;
plyprops.v13 = Effprops.v13;
plyprops.v23 = Effprops.v23;
plyprops.G12 = Effprops.G12;
plyprops.G13 = Effprops.G13;
plyprops.G23 = Effprops.G23;

plyprops.a1 = CTEs(1);  
plyprops.a2 = CTEs(2);  
plyprops.a3 = CTEs(3);  

plyprops.Cstar = Cstar;

switch(Theory)
    case {'Voigt', 'Reuss', 'MT'}
        plyprops.Af = Af;
        plyprops.Am = Am;
        plyprops.ATf = ATf;
        plyprops.ATm = ATm;
    case {'MOC', 'MOCu'}
        plyprops.As = As;
    otherwise
        plyprops.As = As;
        plyprops.RUC = RUC;
end

