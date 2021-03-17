function [MoS, props] = MicroMargins(CriteriaOn, props, ThermoMech, Results_Mech, Results_Th, MoS)
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
% Purpose: Sets up data and calls the function FCriteria to calculate micromechanics-
%          based MoS based on avg constituent or subcell fields
% Input:
% - CriteriaOn: Flags indicating which failure criteria are turned on
% - props: struct containing material properties
% - ThermoMech: Flag indicating if the problem is combined thermomechanical
% - Results_Mech: Struct containing micromechanics analysis results for problem (if 
%                 not ThermoMech) or mechanical problem only (if ThermoMech)
% - Results_Th: Struct containing micromechanics analysis results for the thermal 
%               problem (if ThermoMech), empty (if not ThermoMech)
% - MoS: Struct containing MoS results and info
% Output:
% - MoS: Updated struct containing MoS results and info
% - props: Updated struct containing material properties
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if isfield(props, 'RUC')
    RUC = props.RUC;
else % -- For MOC, create and store 2x2 RUC
    RUC.NB = 2;
    RUC.NG = 2;
    RUC.mats = [props.Fiber.constID, props.Matrix.constID; ...
                props.Matrix.constID, props.Matrix.constID];
    RUC.matsCh = ['F', 'M'; 'M', 'M'];
    RUC.Vf = props.Vf;
    RUC.id = 2;
    RUC.h(1) = sqrt(props.Vf);
    RUC.h(2) = 1-RUC.h(1);
    RUC.l = RUC.h; 
    props.RUC = RUC;
end

% -- If a material doesn't have allowables, set them to +/-1.E99 
if ~isfield(props.Fiber,'allowables')
    [props.Fiber] = SetAllowablesHigh(props.Fiber);
end

if ~isfield(props.Matrix,'allowables')
    [props.Matrix] = SetAllowablesHigh(props.Matrix);
end

if isfield(props,'Interface')
    if ~isfield(props.Interface,'allowables')
        [props.Interface] = SetAllowablesHigh(props.Interface);
    end
end


Micro_Mech = Results_Mech;
if ~ThermoMech % -- Pure thermal or pure mech load case (not ThermoMech)
    if Micro_Mech.Type == "FM" % -- Voigt, Reuss, or MT (avg fields)
        Micro_Th.MicroFieldsF.stress = zeros(6,1);
        Micro_Th.MicroFieldsF.strain = zeros(6,1);
        Micro_Th.MicroFieldsF.thstrain = zeros(6,1);
        Micro_Th.MicroFieldsM.stress = zeros(6,1);
        Micro_Th.MicroFieldsM.strain = zeros(6,1);
        Micro_Th.MicroFieldsM.thstrain = zeros(6,1);
        Micro_Mech.MicroFieldsF.strain = Micro_Mech.MicroFieldsF.strain - Micro_Mech.MicroFieldsF.thstrain;
        Micro_Mech.MicroFieldsM.strain = Micro_Mech.MicroFieldsM.strain - Micro_Mech.MicroFieldsM.thstrain;
    else % -- Subcell fields   
        for b = 1:RUC.NB
            for g = 1:RUC.NG
                Micro_Th.MicroFields(b,g).stress = zeros(6,1);
                Micro_Th.MicroFields(b,g).strain = zeros(6,1);
                Micro_Th.MicroFields(b,g).thstrain = zeros(6,1);
                Micro_Mech.MicroFields(b,g).strain = Micro_Mech.MicroFields(b,g).strain - Micro_Mech.MicroFields(b,g).thstrain;
            end
        end
    end                        
else % -- Combined thermomechanical case
    Micro_Th = Results_Th;
end


if Micro_Mech.Type == "FM" % -- Voigt, Reuss, or MT (avg fields)

     if isfield(props.Fiber,'allowables')
        MicroMechStrainF = Micro_Th.MicroFieldsF.strain - Micro_Th.MicroFieldsF.thstrain; % Max strain criterion based on mechanical strain, not total
        [MoS.MicroF] = FCriteria(Micro_Mech.MicroFieldsF.stress, Micro_Mech.MicroFieldsF.strain, ...
                                 Micro_Th.MicroFieldsF.stress, MicroMechStrainF, props.Fiber.allowables, CriteriaOn, true);
     end
     if isfield(props.Matrix,'allowables')
        MicroMechStrainM = Micro_Th.MicroFieldsM.strain - Micro_Th.MicroFieldsM.thstrain; % Max strain criterion based on mechanical strain, not total
        [MoS.MicroM] = FCriteria(Micro_Mech.MicroFieldsM.stress, Micro_Mech.MicroFieldsM.strain, ...
                                 Micro_Th.MicroFieldsM.stress, MicroMechStrainM, props.Matrix.allowables, CriteriaOn, true);
     end

else % -- Subcell fields
    MicroMechStrain = zeros(6,1);
    for b = 1:RUC.NB
        for g = 1:RUC.NG
            
            if isfield(MoS, 'Micro')
                if isfield(MoS.Micro, 'Failed')
                    if MoS.Micro(b,g).Failed 
                        continue
                    end
                end
            end

            Fid = props.Fiber.constID;
            Mid = props.Matrix.constID;
            if (RUC.mats(b,g) == Fid) && isfield(props.Fiber,'allowables')
                allowables = props.Fiber.allowables;
                MicroMechStrain = Micro_Th.MicroFields(b,g).strain - Micro_Th.MicroFields(b,g).thstrain; % Max strain criterion based on mechanical strain, not total
            elseif (RUC.mats(b,g) == Mid) && isfield(props.Matrix,'allowables')
                allowables = props.Matrix.allowables;
                MicroMechStrain = Micro_Th.MicroFields(b,g).strain - Micro_Th.MicroFields(b,g).thstrain; % Max strain criterion based on mechanical strain, not total
            elseif isfield(props,'Interface')
                Iid = props.Interface.constID;
                if (RUC.mats(b,g) == Iid) && isfield(props.Interface,'allowables')
                    allowables = props.Interface.allowables;
                    MicroMechStrain = Micro_Th.MicroFields(b,g).strain - Micro_Th.MicroFields(b,g).thstrain; % Max strain criterion based on mechanical strain, not total
                end
            else
                continue;
            end
            
            if ~isfield(MoS, 'Micro')
                [MoS.Micro(b,g)] = FCriteria(Micro_Mech.MicroFields(b,g).stress, Micro_Mech.MicroFields(b,g).strain, ...
                                   Micro_Th.MicroFields(b,g).stress, MicroMechStrain, allowables, CriteriaOn, true);
            else
                aa = size(MoS.Micro);
                if aa(1) < b || aa(2) < g
                    [MoS.Micro(b,g)] = FCriteria(Micro_Mech.MicroFields(b,g).stress, Micro_Mech.MicroFields(b,g).strain, ...
                                       Micro_Th.MicroFields(b,g).stress, MicroMechStrain, allowables, CriteriaOn, true);
                else
                    [MoS.Micro(b,g)] = FCriteria(Micro_Mech.MicroFields(b,g).stress, Micro_Mech.MicroFields(b,g).strain, ...
                                       Micro_Th.MicroFields(b,g).stress, MicroMechStrain, allowables, CriteriaOn, true, MoS.Micro(b,g));
                end
            end
            
        end
    end
end  


end

%**************************************************************************
%**************************************************************************

function [Mat] = SetAllowablesHigh(Mat)

% -- Set allowables to +/-1.E99 for materials with no allowabled specified

    Mat.allowables.XT = 1.E99;
    Mat.allowables.XC = -1.E99;
    Mat.allowables.YT = 1.E99;
    Mat.allowables.YC = -1.E99;
    Mat.allowables.ZT = 1.E99;
    Mat.allowables.ZC = -1.E99;
    Mat.allowables.Q = 1.E99;
    Mat.allowables.R = 1.E99;
    Mat.allowables.S = 1.E99;
    Mat.allowables.XeT = 1.E99;
    Mat.allowables.XeC = -1.E99;
    Mat.allowables.YeT = 1.E99;
    Mat.allowables.YeC = -1.E99;
    Mat.allowables.ZeT = 1.E99;
    Mat.allowables.ZeC = -1.E99;
    Mat.allowables.Qe = 1.E99;
    Mat.allowables.Re = 1.E99;
    Mat.allowables.Se = 1.E99;

end
