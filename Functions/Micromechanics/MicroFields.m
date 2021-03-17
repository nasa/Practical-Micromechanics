function [Results] = MicroFields(FullGlobalStrain, DT, props, Results)
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
% Purpose: Calculate the local (constituent level) stress and strain fields in the
%          composite given the concentration tensors and global fields and 
%          temperature change
% Input:
% - FullGlobalStrain: Vector of global total strain components
% - DT: Applied temperature change
% - props: Struct containing effective composite properties
% - Results: Struct containing micromechanics analysis results
% Output:
% - Results: Updated struct containing micromechanics analysis results
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -- Obtain fiber and matrix stiffness tensors
Cf = GetCFromProps(props.Fiber);
Cm = GetCFromProps(props.Matrix);

% -- Micromechanics theories with only average fiber and matrix fields
if (props.micro == "MT" || props.micro == "Voigt" || props.micro == "Reuss")
    ThStrainF = [props.Fiber.aL; props.Fiber.aT; props.Fiber.aT; 0; 0; 0]*DT;
    MicroStrainF = props.Af*FullGlobalStrain + props.ATf*DT;
    MicroStressF = Cf*(MicroStrainF - ThStrainF);
    ThStrainM = [props.Matrix.aL; props.Matrix.aT; props.Matrix.aT; 0; 0; 0]*DT;
    MicroStrainM = props.Am*FullGlobalStrain + props.ATm*DT;
    MicroStressM = Cm*(MicroStrainM - ThStrainM);
    Results.Type = "FM";
    Results.MicroFieldsF.thstrain = ThStrainF;
    Results.MicroFieldsF.strain = MicroStrainF;
    Results.MicroFieldsF.stress = MicroStressF;
    Results.MicroFieldsM.thstrain = ThStrainM;
    Results.MicroFieldsM.strain = MicroStrainM;
    Results.MicroFieldsM.stress = MicroStressM;

% -- MOC with 2x2 subcell RUC
elseif props.micro == "MOC" ||  props.micro == "MOCu"         
    for g = 1:2
        for b = 1:2
            if b == 1 && g == 1
                C = Cf;
                ThStrain = [props.Fiber.aL; props.Fiber.aT; props.Fiber.aT; ...
                            0; 0; 0]*DT;
            else
                C = Cm;
                ThStrain = [props.Matrix.aL; props.Matrix.aT; props.Matrix.aT; ...
                            0; 0; 0]*DT;
            end 
            MicroFields(b,g).strain = props.As(b,g).A*FullGlobalStrain + ...
                                      props.As(b,g).AT*DT;
            MicroFields(b,g).stress = C*(MicroFields(b,g).strain - ThStrain);
            MicroFields(b,g).thstrain = ThStrain;
            Results.Type = "bg";
            Results.MicroFields(b,g) = MicroFields(b,g);
        end
    end

% -- GMC and HFGMC, RUC with NBxNG subcells
elseif (props.micro == "GMC" || props.micro == "HFGMC")               
    NG = props.RUC.NG;
    NB = props.RUC.NB;
    Fid = props.Fiber.constID;
    Mid = props.Matrix.constID;
    if isfield(props,'Interface')
        if isfield(props.Interface,'constID')
            Iid = props.Interface.constID;
        else
            Iid = 0;
        end
    else
        Iid = 0;
    end
    RUC = props.RUC;
     for g = 1:NG
        for b = 1:NB
            C = props.As(b,g).C;
            ThStrain = props.As(b,g).alpha*DT;
            MicroFields(b,g).strain = props.As(b,g).A*FullGlobalStrain + ...
                                      props.As(b,g).AT*DT;
            MicroFields(b,g).stress = C*(MicroFields(b,g).strain - ThStrain);
            MicroFields(b,g).thstrain = ThStrain;
            Results.Type = "bg";
            Results.MicroFields(b,g) = MicroFields(b,g);
        end
     end

else
     txt = ['plyprops{',char(num2str(mat)),'}.micro = ', char(props.micro)];
     error(['** Incorrect micromechanics theory name specified: ', txt])

end 

end