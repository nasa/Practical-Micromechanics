function [LamResults, plyprops] = CalcLamMargins(plyprops, Geometry, Loads, LamResults, Damage)
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
% Purpose: Divides laminates thermomechanical problems into separate thermal and 
%          mechanical load cases, solves CLT, and determines the local fields (if 
%          applicable) for use in calculating margins of safety
% Input:
% - plyprops: Cell/struct containing effective composite material properties 
% - Geometry: Struct containing laminate definition variables
% - Loads: Struct containing problem loading information
% - LamResults: Struct containing laminate analysis results
% - Damage: Flag indicating if problem involves prog. damage (Chapter 7)
% Output:
% - LamResults: Updated struct containing laminate analysis results
% - plyprops: Updated cell/struct containing effective composite material 
%             properties 
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -- If Damage, do not treat thermal as preload
if nargin == 4
    Damage = false;
end

% -- Determine if need to do MoS calcs based on presence of allowables
DoPlyMargins = false;
DoMicroMargins = false;

% -- Loop through plies
for k = Geometry.N
    mat = Geometry.plymat(k);
    
    % -- Check for ply-level allowables
    if isfield(plyprops{mat},'allowables')
        DoPlyMargins = true;
    end
    
    % -- Check for constituent allowables
    if isfield(plyprops{mat},'micro')
        if isfield(plyprops{mat}.Fiber,'allowables') || isfield(plyprops{mat}.Matrix,'allowables')
            DoMicroMargins = true;
        elseif isfield(plyprops, 'Interface')
            if isfield(plyprops.Interface, 'allowables')
                DoMicroMargins = true;
            end
        end
    end
end

% -- Return if no allowables specified
if ~DoPlyMargins && ~DoMicroMargins
   LamResults.MicroMoS = false;
   return 
end

% -- Determine if this is a thermomechanical problem
if Loads.DT ~= 0 && any(Loads.Value) && ~Damage
    
    ThermoMech = true;
    
    % -- Calculate ply-level & local fields for ** THERMAL PROBLEM ** only
    Loads_Th = Loads;
    Loads_Th.DT = Loads.DT;
    Loads_Th.Type = Loads.Type; % -- Must apply same types so loads superimpose
    Loads_Th.Value = zeros(1,6);
    
    % -- Analyze laminate with CLT
    [Geometry, LamResults_Th] = CLT(plyprops, Geometry, Loads_Th);

    % -- Calculate local fields for micromechanics
    [LamResults_Th] = LamMicroFields(Geometry, Loads_Th, plyprops, LamResults_Th);

    % -- Calculate ply-level & local fields for ** MECHANICAL PROBLEM ** problem only    
    Loads_Mech = Loads;
    Loads_Mech.DT = 0;
    
    % -- Analyze laminate with CLT
    [Geometry, LamResults_Mech] = CLT(plyprops, Geometry, Loads_Mech);

    % -- Calculate local fields for micromechanics
    [LamResults_Mech] = LamMicroFields(Geometry, Loads_Mech, plyprops, LamResults_Mech);    

% -- Pure thermal or pure mechechanical problem or prog damage
%    Set thermal case results to zero and mech case results to reg results
else    
    LamResults_Mech = LamResults;
    LamResults_Th.MCstress = zeros(6,2*Geometry.N);
    LamResults_Th.MCstrain = zeros(6,2*Geometry.N);
    LamResults_Th.Micro(1:2*Geometry.N) = 0;
    ThermoMech = false;
    
end


% -- Calculate ** PLY LEVEL ** margins
if DoPlyMargins
    MoSControllingZPt = 0;
    MinMoS = 99999;
    
    % -- Loop through points in each ply
    for kk = 1:2*Geometry.N

        k = round(kk/2);
        mat = Geometry.plymat(k);

        if isfield(plyprops{mat},'allowables')

            stress = zeros(6,1);
            stress(1) = LamResults_Mech.MCstress(1, kk);
            stress(2) = LamResults_Mech.MCstress(2, kk);
            stress(6) = LamResults_Mech.MCstress(3, kk);

            strain = zeros(6,1);
            strain(1) = LamResults_Mech.MCstrain(1, kk);
            strain(2) = LamResults_Mech.MCstrain(2, kk);
            strain(6) = LamResults_Mech.MCstrain(3, kk);

            stressthLC = zeros(6,1);
            stressthLC(1) = LamResults_Th.MCstress(1, kk);
            stressthLC(2) = LamResults_Th.MCstress(2, kk);
            stressthLC(6) = LamResults_Th.MCstress(3, kk);

            strainthLC = zeros(6,1);
            strainthLC(1) = LamResults_Th.MCstrain(1, kk);
            strainthLC(2) = LamResults_Th.MCstrain(2, kk);
            strainthLC(6) = LamResults_Th.MCstrain(3, kk);
            

            thstrain = zeros(6,1);
            thstrain(1) = plyprops{mat}.a1*Loads.DT;
            thstrain(2) = plyprops{mat}.a2*Loads.DT;

            % -- Max strain criterion based on mechanical strain, not total
            mechstrainthLC = strainthLC - thstrain;

            % -- Calculate margins    
            [MoS(kk)] = FCriteria(stress, strain, stressthLC, mechstrainthLC, plyprops{mat}.allowables, Loads.CriteriaOn, false);
            
            % -- Check if this point is controlling so far
            if MoS(kk).MinMoS < MinMoS
                MinMoS = MoS(kk).MinMoS;
                MoSControllingZPt = kk;
            end

        end

    end

    LamResults.MoS = MoS;
    LamResults.MoSControllingZPt = MoSControllingZPt;

end


% -- Calculate ** MICRO SCALE ** margins
LamResults.MicroMoS = false;
if DoMicroMargins

    % -- Loop through points in each ply
    for kk = 1:2*Geometry.N
       
        k = round(kk/2);
        mat = Geometry.plymat(k);

        if isfield(plyprops{mat},'micro')
            
            % -- Retain LamResults.Micro(kk) if it exists
            if isfield(LamResults.Micro(kk), 'MoS') 
                MoS = LamResults.Micro(kk).MoS;
            else
                MoS = struct;
            end
            
            % -- Calculate micro scale MoS for this point in laminate
            [MoS, plyprops{mat}] = MicroMargins(Loads.CriteriaOn, plyprops{mat}, ThermoMech, LamResults_Mech.Micro(kk), LamResults_Th.Micro(kk), MoS);
            LamResults.Micro(kk).MoS = MoS;

            if isfield(LamResults.Micro(kk), 'MoS')
                LamResults.MicroMoS = true;
            end

        end
    end

end
