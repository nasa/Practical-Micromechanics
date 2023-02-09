function [NP, OutInfo, MicroMat, Loads] = MicroProblemDef()
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
% Purpose: Sets up and defines stand-alone micromechanics problems.  Specifies
%          problem name, material to be analyzed, loading, and failure criteria.
% Output:
% - NP: Number of problems
% - OutInfo: Struct containing output information
% - MicroMat: Effective material number (from GetEff props) to be analyzed
%             per problem
% - Loads: Cell array containing problem loading information per problem
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

NProblemsMax = 200; % -- To preallocate cells
NP = 0; % -- Problem number counter
E = 1;  % -- Applied loading type identifier (strains)
S = 2;  % -- Applied loading type identifier (stresses)

% -- Preallocate cell
Loads = cell(1,NProblemsMax);

% -- Type of output requested
OutInfo.Format = "txt"; % -- Text format
% OutInfo.Format = "doc";  % -- MS Word format (Windows only, MUCH slower)

% -- Include or exclude plots from output
OutInfo.MakePlots = true;
% OutInfo.MakePlots = false;

% -- Turn on/off writing margins to file (can take long for doc output)
% OutInfo.WriteMargins = true;
OutInfo.WriteMargins = false;

%========================================================================
% -- Define all problems below
%========================================================================


%% -- Chapter 3 Standalone Micromechanics Examples

NP = NP + 1;
OutInfo.Name(NP) = "Sect 3.12.1&2 - Voigt IM7-8552 Vf=0.55";
MicroMat{NP} = 100;
Loads{NP}.DT = 0.0;
Loads{NP}.Type  = [S,  S,  S,  S,  S,  S];
Loads{NP}.Value = [0,  1,  0,  0,  0,  0];

NP = NP + 1;
OutInfo.Name(NP) = "Sect 3.12.1&2 - Reuss IM7-8552 Vf=0.55";
MicroMat{NP} = 101;
Loads{NP}.DT = 0.0;
Loads{NP}.Type  = [S,  S,  S,  S,  S,  S];
Loads{NP}.Value = [0,  1,  0,  0,  0,  0];

NP = NP + 1;
OutInfo.Name(NP) = "Sect 3.12.1&2 - MT IM7-8552 Vf=0.55";
MicroMat{NP} = 102;
Loads{NP}.DT = 0.0;
Loads{NP}.Type  = [S,  S,  S,  S,  S,  S];
Loads{NP}.Value = [0,  1,  0,  0,  0,  0];

NP = NP + 1;
OutInfo.Name(NP) = "Sect 3.12.1&2 - MOC IM7-8552 Vf=0.55";
MicroMat{NP} = 103;
Loads{NP}.DT = 0.0;
Loads{NP}.Type  = [S,  S,  S,  S,  S,  S];
Loads{NP}.Value = [0,  1,  0,  0,  0,  0];

NP = NP + 1;
OutInfo.Name(NP) = "Sect 3.12.1&2 - MOCu IM7-8552 Vf=0.55";
MicroMat{NP} = 104;
Loads{NP}.DT = 0.0;
Loads{NP}.Type  = [S,  S,  S,  S,  S,  S];
Loads{NP}.Value = [0,  1,  0,  0,  0,  0];

% NP = NP + 1;
% OutInfo.Name(NP) = "Sect 3.12.3 - All theories IM7-8552 all Vf";
% MicroMat{NP} = 105;
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Sect 3.12.3 - All theories glass-epoxy all Vf";
% MicroMat{NP} = 106;
% 
%% -- Chapter 4 Standalone Micromechanics Examples
%
% NP = NP + 1;
% OutInfo.Name(NP) = "Section 4.4.2.1 - MT IM7-8552 Vf=0.55 - S11";
% MicroMat{NP} = 102;
% Loads{NP}.DT = 0;
% Loads{NP}.Type  = [   S, S,  S,  S,  S,  S];
% Loads{NP}.Value = [1000, 0,  0,  0,  0,  0];
% Loads{NP}.CriteriaOn = [1, 1, 1, 1]; 
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Section 4.4.2.1 - MOC IM7-8552 Vf=0.55 - S11";
% MicroMat{NP} = 103;
% Loads{NP}.DT = 0;
% Loads{NP}.Type  = [   S, S,  S,  S,  S,  S];
% Loads{NP}.Value = [1000, 0,  0,  0,  0,  0];
% Loads{NP}.CriteriaOn = [1, 1, 1, 1]; 
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Section 4.4.2.1&2 - MT IM7-8552 Vf=0.55 - S22";
% MicroMat{NP} = 102;
% Loads{NP}.DT = 0;
% Loads{NP}.Type  = [ S,  S,  S,  S,  S,  S];
% Loads{NP}.Value = [ 0, 10,  0,  0,  0,  0];
% Loads{NP}.CriteriaOn = [1, 1, 1, 1]; 
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Section 4.4.2.1&2 - MOC IM7-8552 Vf=0.55 - S22";
% MicroMat{NP} = 103;
% Loads{NP}.DT = 0;
% Loads{NP}.Type  = [ S,  S,  S,  S,  S,  S];
% Loads{NP}.Value = [ 0, 10,  0,  0,  0,  0];
% Loads{NP}.CriteriaOn = [1, 1, 1, 1]; 
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Section 4.4.2.3 - MT IM7-8552 Vf=0.55 - S22 + DT";
% MicroMat{NP} = 102;
% Loads{NP}.DT = 100;
% %Loads{NP}.DT = -100;
% Loads{NP}.Type  = [ S,  S,  S,  S,  S,  S];
% Loads{NP}.Value = [ 0, 10,  0,  0,  0,  0];
% Loads{NP}.CriteriaOn = [1, 1, 1, 1]; 
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Section 4.4.2.3 - MOC IM7-8552 Vf=0.55 - S22 + DT";
% MicroMat{NP} = 103;
% %Loads{NP}.DT = 100;
% Loads{NP}.DT = -100;
% Loads{NP}.Type  = [ S,  S,  S,  S,  S,  S];
% Loads{NP}.Value = [ 0, 10,  0,  0,  0,  0];
% Loads{NP}.CriteriaOn = [1, 1, 1, 1]; 
% 
% % -- Example micro scale failure envelope - Section 4.4.3
% tttt = datetime(datetime,'Format','yyyy-MMM-dd HH.mm.ss');
% OutInfo.EnvFile = ['Output/', 'Composite ENVELOPE - ', char(tttt), '.EnvData'];
% ang = 0;
% anginc = 5;
% anginc1 = anginc;
% 
% while (ang <= 360.1)
% 
%     anginc = anginc1;
%     
%     % -- Extra angles for highly anisotropic cases
% %     if ang < 15 || (ang >= 165 && ang <= 190) || ang >= 345
% %       anginc = 0.2;
% %     end
% 
%     NP = NP + 1;
%     OutInfo.Name(NP) = string(['Ang = ', num2str(ang)]);
%     MicroMat{NP} = 102; % -- MT
% %    MicroMat{NP} = 103; % --- MOC
%     
%     S33 = sin(ang*pi/180);
%     S22 = cos(ang*pi/180);
%     % Define laminate applied loading and temperature change
%     Loads{NP}.S = S;
%     Loads{NP}.E = E;
%     Loads{NP}.Type  = [S,    S,   S,  S,  S, S];
%     Loads{NP}.Value = [0,  S22, S33,  0,  0,  0];
% 
%     Loads{NP}.DT = 0.0;    
% 
%     % Turn off individual criteria
%     Loads{NP}.CriteriaOn = [1,0,0,0]; % -- Max stress
%     %Loads{NP}.CriteriaOn = [0,1,0,0]; % -- Max strain
%     %Loads{NP}.CriteriaOn = [0,0,0,1]; % -- Tsai-Wu
%     Loads{NP}.ang = ang;
%     ang = ang + anginc;
% end
% 
%% -- Chapter 5 Standalone Micromechanics Examples
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Tsai-Hahn Glass/Epoxy Chap 5 Example 1";
% MicroMat{NP} = 130;
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Dean & Turner Carbon/Epoxy Chap 5 Example 2";
% MicroMat{NP} = 131;
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Contour Glass-Epoxy Chap 5 Fig 5.10";
% MicroMat{NP} = 132;
% Loads{NP}.DT = 0;
% Loads{NP}.Type  = [  S,   S,  E,   S,  S,  S];
% Loads{NP}.Value = [  0, 0.0,  0.02,   0,  0,  0];
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "SiC-SiC No Interface Chap 5 Table 5.3";
% MicroMat{NP} = 133;
% Loads{NP}.DT = 0;
% Loads{NP}.Type  = [  S,   S,    E,  S,  S,  S];
% Loads{NP}.Value = [  0, 0.0, 0.02,  0,  0,  0];
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "SiC-SiC BN Interface Chap 5 Fig 5.11";
% MicroMat{NP} = 134;
% Loads{NP}.DT = 0;
% Loads{NP}.Type  = [  S,      E,  S,   S,  S,  S];
% Loads{NP}.Value = [  0, 0.0006,  0.,  0,  0,  0];
% 
% % -- Micro scale failure envelopes - Example Fig 5.12 
% tttt = datetime(datetime,'Format','yyyy-MMM-dd HH.mm.ss');
% OutInfo.EnvFile = ['Output/', 'SiC-SiC-BN Composite ENVELOPE - ', char(tttt), '.EnvData'];
% ang = 0;
% anginc = 2;
% anginc1 = anginc;
% 
% while (ang <= 360.1)
% 
%     anginc = anginc1;
%     
%     % -- Extra angles for highly anisotropic cases
%     if ang < 15 || (ang >= 165 && ang <= 190) || ang >= 345
%       anginc = 0.2;
%     end
% 
%     NP = NP + 1;
%     OutInfo.Name(NP) = string(['Ang = ', num2str(ang)]);
%     %MicroMat{NP} = 133; % -- 2x2 GMC - NO interface
%     MicroMat{NP} = 134; % -- 5x5 GMC - RUCid = 105
%     
%     S22 = sin(ang*pi/180);
%     S11 = cos(ang*pi/180);
%     Loads{NP}.Type  = [S,   S,   S,  S,  S, S];
%     Loads{NP}.Value = [S11, S22, 0,  0,  0,  0];
% 
%     Loads{NP}.CriteriaOn = [1,1,1,1];
%     Loads{NP}.ang = ang;
%     ang = ang + anginc;
% end

%% -- Chapter 6 Standalone Micromechanics Examples

% NP = NP + 1;
% OutInfo.Name(NP) = "Example Section 6.10.1 - 26x26 GMC - Sun & Vaidya Carbon Epoxy";
% MicroMat{NP} = 140;
% Loads{NP}.DT = 0;
% Loads{NP}.Type  = [  S,  S,  S,  S,  S,  S];
% Loads{NP}.Value = [  0,  0,  0,  0,  0,  0];
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Example Section 6.10.1 - 26x26 HFGMC - Sun & Vaidya Carbon Epoxy";
% MicroMat{NP} = 141;
% Loads{NP}.DT = 0;
% Loads{NP}.Type  = [  S,  S,  S,  S,  S,  S];
% Loads{NP}.Value = [  0,  0,  0,  0,  0,  0];
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Example Section 6.10.1 - 26x26 GMC vs. HFGMC - glass-epoxy";
% MicroMat{NP} = 142;
% Loads{NP}.DT = 0;
% Loads{NP}.Type  = [  S,  S,  S,  S,  S,  S];
% Loads{NP}.Value = [  0,  0,  0,  0,  0,  0];
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Example Section 6.10.2 - 32x32 HFGMC - contour glass-epoxy";
% MicroMat{NP} = 143;
% Loads{NP}.DT = 0;
% Loads{NP}.Type  = [  S,  S,    E,  S,  S,  S];
% Loads{NP}.Value = [  0,  0, 0.02,  0,  0,  0];
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Example Section 6.10.2 - 100x100 HFGMC - contour glass-epoxy";
% MicroMat{NP} = 144;
% Loads{NP}.DT = 0;
% Loads{NP}.Type  = [  S,  S,    E,  S,  S,  S];
% Loads{NP}.Value = [  0,  0, 0.02,  0,  0,  0];
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Example Section 6.10.3 - 26x26 Sq GMC - contour glass-epoxy";
% MicroMat{NP} = 145;
% Loads{NP}.DT = 0;
% Loads{NP}.Type  = [  S,   S,  S,  S,  S,  S];
% Loads{NP}.Value = [  0, 100,  0,  0,  0,  0];
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Example Section 6.10.3 - 32x32 Sq HFGMC - contour glass-epoxy";
% MicroMat{NP} = 143;
% Loads{NP}.DT = 0;
% Loads{NP}.Type  = [  S,   S,  S,  S,  S,  S];
% Loads{NP}.Value = [  0, 100,  0,  0,  0,  0];
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Example Section 6.10.3 - 68x120 Hex GMC - contour glass-epoxy";
% MicroMat{NP} = 146;
% Loads{NP}.DT = 0;
% Loads{NP}.Type  = [  S,   S,  S,  S,  S,  S];
% Loads{NP}.Value = [  0, 100,  0,  0,  0,  0];
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Example Section 6.10.3 - 68x120 Hex HFGMC - contour glass-epoxy";
% MicroMat{NP} = 147;
% Loads{NP}.DT = 0;
% Loads{NP}.Type  = [  S,   S,  S,  S,  S,  S];
% Loads{NP}.Value = [  0, 100,  0,  0,  0,  0];
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Example Section 6.10.3 - 79x79 Random (json) GMC - contour glass-epoxy";
% MicroMat{NP} = 148;
% Loads{NP}.DT = 0;
% Loads{NP}.Type  = [  S,   S,  S,  S,  S,  S];
% Loads{NP}.Value = [  0, 100,  0,  0,  0,  0];
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Example Section 6.10.3 - 79x79 Random (json) HFGMC - contour glass-epoxy";
% MicroMat{NP} = 149;
% Loads{NP}.DT = 0;
% Loads{NP}.Type  = [  S,   S,  S,  S,  S,  S];
% Loads{NP}.Value = [  0, 100,  0,  0,  0,  0];
%
% NP = NP + 1;
% OutInfo.Name(NP) = "Example Section 6.10.4 - 136x236 Hex SiC-SiC HFGMC - S22";
% MicroMat{NP} = 150;
% Loads{NP}.DT = 0;
% Loads{NP}.Type  = [  S,   S,  S,  S,  S,  S];
% Loads{NP}.Value = [  0, 100,  0,  0,  0,  0];
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Example Section 6.10.4 - 136x236 Hex SiC-SiC HFGMC - S33";
% MicroMat{NP} = 150;
% Loads{NP}.DT = 0;
% Loads{NP}.Type  = [  S,  S,   S,  S,  S,  S];
% Loads{NP}.Value = [  0,  0, 100,  0,  0,  0];
% % 
% NP = NP + 1;
% OutInfo.Name(NP) = "Example Section 6.10.5 - HFGMC Subcell Grid Convergence";
% % -- Target_NG must be altered in GetRUC.m to change subcell grid
% MicroMat{NP} = 151;
% Loads{NP}.DT = 0;
% Loads{NP}.Type  = [  S,   S,  S,  S,  S,  S];
% Loads{NP}.Value = [  0, 100,  0,  0,  0,  0];
% Loads{NP}.CriteriaOn = [1,0,0,0];

% % -- This must be run with RandomRUCHist.m (with no other problems defined herein)
% NUMBER = 1000; % -- Note, it will take a long time to run 1000 random cases
% for NP1 = 1: NUMBER 
%     NP = NP + 1;
%     OutInfo.Name(NP) = "Example Section 6.10.6 - HFGMC Random 1000 Instances";
%     MicroMat{NP} = 152;
%     Loads{NP}.DT = 0;
%     Loads{NP}.Type  = [  S,  S,  S,  S,  S,  S];
%     Loads{NP}.Value = [  0,  0,  0,  0,  0,  0];
% end
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Example Section 6.10.7 - GMC";
% MicroMat{NP} = 153;
% Loads{NP}.DT = 0;
% Loads{NP}.Type  = [  S,  S,  S,   S,  S,  S];
% Loads{NP}.Value = [  0,  0,  0, 100,  0,  0];
% Loads{NP}.CriteriaOn = [1,0,0,0];
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Example Section 6.10.7 - HFGMC";
% MicroMat{NP} = 154;
% Loads{NP}.DT = 0;
% Loads{NP}.Type  = [  S,  S,  S,   S,  S,  S];
% Loads{NP}.Value = [  0,  0,  0, 100,  0,  0];
% Loads{NP}.CriteriaOn = [1,0,0,0];

%% -- Chapter 7 Standalone Micromechanics Examples
%
% -- 7.1.2 - GMC 2x2 progressive damage
% NP = NP + 1;
% OutInfo.Name(NP) = "Section 7.1.2 - GMC 2x2 progressive damage";
% MicroMat{NP} = 171; % -- GMC 2x2
% Loads{NP}.Type  = [ S,    E,  S,  S,  S,  S];
% Loads{NP}.Value = [ 0, 0.02,  0,  0,  0,  0];
% Loads{NP}.CriteriaOn = [1, 0, 0, 0]; 
% Loads{NP}.NINC = 400;
% 
% % -- 7.2 - Epoxy with notch
% MicroMat{NP + 1} = 172; % -- HF 25x25
% MicroMat{NP + 2} = 173; % -- HF 51x51
% MicroMat{NP + 3} = 174; % -- HF 101x101
% MicroMat{NP + 4} = 175; % -- HF 201x201
% 
% for I = 1:4
%     NP = NP + 1;
%     OutInfo.Name(NP) = "Section 7.2 - Epoxy with notch";
%     Loads{NP}.Type  = [S,    E,    S,    S,  S,   S];
%     Loads{NP}.Value = [0,    0.02, 0,    0,  0,   0];
%     Loads{NP}.NINC = 400;
%     Loads{NP}.CriteriaOn = [1,0,0,0];
% end
% 
% % -- 7.3 - Influence of Micromechanics Theory on Consistent Constituent Strength
% MicroMat{NP + 1} = 176; % -- MT with MT props
% MicroMat{NP + 2} = 177; % -- HF with MT props
% MicroMat{NP + 3} = 178; % -- MT with HF props
% MicroMat{NP + 4} = 179; % -- HF with HF props
% 
% for NP1 = 1:4
%     NP = NP + 1;
%     OutInfo.Name(NP) = "Section 7.3";
%     Loads{NP}.DT = 0.0;
%     Loads{NP}.Type  = [S,     E,    S,    S,    S,   S];
%     Loads{NP}.Value = [0,     0.04,    0,  0,    0,     0];
%     Loads{NP}.NINC = 800;
%     Loads{NP}.CriteriaOn = [1,0,0,0];
% end

% % -- 7.4.1.1 - Unidirectional PMC - Comparison of Progressive Damage using MT, GMC, and HFGMC
% MicroMat{NP + 1} = 178; % -- MT (HF props)
% MicroMat{NP + 2} = 180; % -- GMC 32x32
% MicroMat{NP + 3} = 181; % -- HFGMC 32x32
% 
% for NP1 = 1:3
%     NP = NP + 1;
%     OutInfo.Name(NP) = "Section 7.4.1.1 - Axial";
%     Loads{NP}.Type  = [E,     S,    S,    S,    S,   S];
%     Loads{NP}.Value = [0.04,  0,    0,    0,    0,   0];
%     Loads{NP}.NINC = 200;
%     Loads{NP}.CriteriaOn = [1,0,0,0];
% end
% 
% MicroMat{NP + 1} = 178; % -- MT (HF props)
% MicroMat{NP + 2} = 180; % -- GMC 32x32
% MicroMat{NP + 3} = 181; % -- HFGMC 32x32
% 
% for NP1 = 1:3
%     NP = NP + 1;
%     OutInfo.Name(NP) = "Section 7.4.1.1 - Trans";
%     Loads{NP}.Type  = [S,     E,    S,    S,    S,   S];
%     Loads{NP}.Value = [0,  0.03,    0,    0,    0,   0];
%     Loads{NP}.NINC = 600;
%     Loads{NP}.CriteriaOn = [1,0,0,0];
% end
% 
% MicroMat{NP + 1} = 178; % -- MT (HF props)
% MicroMat{NP + 2} = 180; % -- GMC 32x32
% MicroMat{NP + 3} = 181; % -- HFGMC 32x32
% 
% for NP1 = 1:3
%     NP = NP + 1;
%     OutInfo.Name(NP) = "Section 7.4.1.1 - Shear 12";
%     Loads{NP}.Type  = [S,  S,    S,    S,    S,   E];
%     Loads{NP}.Value = [0,  0,    0,    0,    0,   0.15];
%     Loads{NP}.NINC = 3000;
%     Loads{NP}.CriteriaOn = [1,0,0,0];
% 
% end
% 
% MicroMat{NP + 1} = 178; % -- MT (HF props)
% MicroMat{NP + 2} = 180; % -- GMC 32x32
% MicroMat{NP + 3} = 181; % -- HFGMC 32x32
% 
% for NP1 = 1:3
%     NP = NP + 1;
%     OutInfo.Name(NP) = "Section 7.4.1.1 - Shear 23";
%     Loads{NP}.Type  = [S,  S,    S,    E,    S,   S];
%     Loads{NP}.Value = [0,  0,    0,    0.1,  0,   0];
%     Loads{NP}.NINC = 2000;
%     Loads{NP}.CriteriaOn = [1,0,0,0];
% end
% 
% -- 7.4.1.2, Part A - Unidirectional PMC - Effect of fiber packing arrangement on progressive damage
% MicroMat{NP + 1} = 181; % -- HFGMC 32x32
% MicroMat{NP + 2} = 179; % -- HFGMC Hex 46x80
% MicroMat{NP + 3} = 179; % -- HFGMC Hex 46x80
% MicroMat{NP + 4} = 184; % -- HFGMC Random-3 from json
% 
% for NP1 = 1:4
%     NP = NP + 1;
%     OutInfo.Name(NP) = "Section 7.4.1.2 - Packing Part A";
% 
%     Loads{NP}.S = S;
%     Loads{NP}.E = E;
%     if NP1 == 3 % -- eps33 applied
%         Loads{NP}.Type  = [S,  S,    E,    S,  S,   S];
%         Loads{NP}.Value = [0,  0,    0.02, 0,  0,   0];
%     else % -- eps22 applied
%         Loads{NP}.Type  = [S,  E,    S,    S,  S,   S];
%         Loads{NP}.Value = [0,  0.02, 0,    0,  0,   0];
%     end
%     Loads{NP}.NINC = 400;
%     Loads{NP}.CriteriaOn = [1,0,0,0];
% 
% end
% 
% % -- 7.4.1.2, Part B - Unidirectional PMC - NUMBER of random fiber packing cases
% NUMBER = 3; % -- These runs are relative long
% for I = 1: NUMBER
%     
%     NP = NP + 1;
% 
%     OutInfo.Name(NP) = "Section 7.4.1.2, Part B - many random";
%     MicroMat{NP} = 183; % -- HF random
% 
%     Loads{NP}.S = S;
%     Loads{NP}.E = E;
%     Loads{NP}.Type  = [S,     E,    S,    S,    S,   S];
%     Loads{NP}.Value = [0,     0.06,    0,  0,    0,     0];
%     Loads{NP}.TerminationFactor = 5;
%     Loads{NP}.NINC = 1200;
%     Loads{NP}.CriteriaOn = [1,0,0,0];
% 
% end
% 
% -- 7.4.1.3 - Unidirectional PMC - Number of fibers in RUC
% for NP1 = 1:6 
%     NP = NP + 1;
%     % -- Random RUC, need to alter GetRUC.m - Radf = 10, alter Nfibers = 4, 8, 16 
%     OutInfo.Name(NP) = "Section 7.4.1.3 - Number of fibers";
%     MicroMat{NP} = 183; % -- HF random
%     Loads{NP}.Type  = [S,    E,  S,  S,  S,   S];
%     Loads{NP}.Value = [0, 0.04,  0,  0,  0,   0];
%     Loads{NP}.NINC = 800;
%     Loads{NP}.CriteriaOn = [1,0,0,0];
% end
% 
% 
% % -- 7.4.1.4 - Unidirectional PMC - Failure Criterion, Random-1
% Loads{NP + 1}.CriteriaOn = [1,0,0,0];
% Loads{NP + 2}.CriteriaOn = [0,1,0,0];
% Loads{NP + 3}.CriteriaOn = [0,0,1,0];
% Loads{NP + 4}.CriteriaOn = [0,0,0,1];
% 
% for NP1 = 1:4
%     
%     NP = NP + 1;
%     OutInfo.Name(NP) = "Section 7.4.1.4 - Failure Criterion, Random-1";
%     MicroMat{NP} = 182; % -- HFGMC, choose Random-1 from json
% 
%     Loads{NP}.S = S;
%     Loads{NP}.E = E;
%     Loads{NP}.Type  = [S,  E,    S,    S,  S,   S];
%     Loads{NP}.Value = [0,  0.02, 0,    0,  0,   0];
%     Loads{NP}.NINC = 400;
% 
% end
% 
% % -- 7.4.1.4 - Unidirectional PMC - Failure Criterion, Random-3
% Loads{NP + 1}.CriteriaOn = [1,0,0,0];
% Loads{NP + 2}.CriteriaOn = [0,1,0,0];
% Loads{NP + 3}.CriteriaOn = [0,0,1,0];
% Loads{NP + 4}.CriteriaOn = [0,0,0,1];
% 
% for NP1 = 1:4
%     
%     NP = NP + 1;
%     OutInfo.Name(NP) = "Section 7.4.1.4 - Failure Criterion, Random 3";
%     MicroMat{NP} = 184; % -- HFGMC, choose Random-3 from json
% 
%     Loads{NP}.S = S;
%     Loads{NP}.E = E;
%     Loads{NP}.Type  = [S,  E,    S,    S,  S,   S];
%     Loads{NP}.Value = [0,  0.02, 0,    0,  0,   0];
%     Loads{NP}.NINC = 400;
% 
% end
% 
% % -- 7.4.1.5 - Unidirectional PMC - Effect of loading increment size, max stress
% Loads{NP + 1}.NINC = 50;  % -- e_inc = 0.0003
% Loads{NP + 2}.NINC = 100; % -- e_inc = 0.00015
% Loads{NP + 3}.NINC = 150; % -- e_inc = 0.0001
% Loads{NP + 4}.NINC = 300; % -- e_inc = 0.00005
% 
% for NP1 = 1:4
%     NP = NP + 1;
%     OutInfo.Name(NP) = "Section 7.4.1.5 - Increment Size, max stress";
%     MicroMat{NP} = 181; % -- HFGMC 32x32
%     Loads{NP}.S = S;
%     Loads{NP}.E = E;
%     Loads{NP}.Type  = [S,  E,    S,    S,  S,   S];
%     Loads{NP}.Value = [0,  0.015, 0,    0,  0,   0];
%     Loads{NP}.CriteriaOn = [1,0,0,0];
% end
% 
% % -- 7.4.1.5 - Unidirectional PMC - Effect of loading increment size, Tsai-Wu
% Loads{NP + 1}.NINC = 50;  % -- e_inc = 0.0003
% Loads{NP + 2}.NINC = 100; % -- e_inc = 0.00015
% Loads{NP + 3}.NINC = 150; % -- e_inc = 0.0001
% Loads{NP + 4}.NINC = 300; % -- e_inc = 0.00005
% 
% for NP1 = 1:4
%     NP = NP + 1;
%     OutInfo.Name(NP) = "Section 7.4.1.5 - Increment Size, Tsai-Wu";
%     MicroMat{NP} = 181; % -- HFGMC 32x32
%     Loads{NP}.S = S;
%     Loads{NP}.E = E;
%     Loads{NP}.Type  = [S,  E,    S,    S,  S,   S];
%     Loads{NP}.Value = [0,  0.015, 0,    0,  0,   0];
%     Loads{NP}.CriteriaOn = [0,0,0,1];
% end
% 
% % -- 7.4.1.6 - Unidirectional PMC - Effect of geometric discretization
% NP = NP + 1;
% OutInfo.Name(NP) = "Section 7.4.1.6 - geometric discretization";
% % -- Random RUC, need to alter GetRUC.m - Nfibers = 1, alter Radf = 5, 10, 20, 40, 80 
% MicroMat{NP} = 183; % -- HF random
% Loads{NP}.Type  = [S,  E,    S,    S,  S,   S];
% Loads{NP}.Value = [0,  0.01, 0,    0,  0,   0];
% Loads{NP}.NINC = 400;
% Loads{NP}.CriteriaOn = [0,0,0,1];
% 
% % -- 7.5.1.1 - Unidirectional CMC - Effect of loading direction CMC
% % -- Alter GetRUC.m RUCid = 200 -- Target_Ni = 2, Target_NG = 58
% for NP1 = 1:6
%     NP = NP + 1;
%     OutInfo.Name(NP) = "Section 7.5.1.1 - Loading Orientation, CMC Hex";
%     MicroMat{NP} = 189; % -- HFGMC 32x58 Hex with interface
%     Loads{NP}.Type = [S, S, S, S, S, S];
%     Loads{NP}.Value = zeros(6, 1);
%     Loads{NP}.Type(NP1) = E;
%     Loads{NP}.Value(NP1) = 0.01;
%     Loads{NP}.NINC = 200;
%     Loads{NP}.CriteriaOn = [1,0,0,0];
% end
% 
% % -- 7.5.1.2 - Unidirectional CMC - Effect of order vs disorder CMC
% MicroMat{NP + 1} = 190; % -- HFGMC Random-A json with int failure
% MicroMat{NP + 2} = 191; % -- HFGMC Random-B json with int failure
% MicroMat{NP + 3} = 192; % -- HFGMC Random=C json with int failure
% MicroMat{NP + 4} = 189; % -- HFGMC 32x58 Hex with int failure
% for NP1 = 1:4
%     NP = NP + 1;
%     OutInfo.Name(NP) = "Section 7.5.1.2 - CMC Random vs. Hex";
%     Loads{NP}.Type  = [S,  E,     S,    S,  S,   S];
%     Loads{NP}.Value = [0,  0.004, 0,    0,  0,   0];
%     Loads{NP}.NINC = 80;
%     Loads{NP}.CriteriaOn = [1,0,0,0];
% end
% 
% % -- 7.5.1.3 - Unidirectional CMC - Effect of interface failure
% MicroMat{NP + 1} = 193; % -- HFGMC Random-A json, NO int failure
% MicroMat{NP + 2} = 194; % -- HFGMC Random-B json, NO int failure
% MicroMat{NP + 3} = 195; % -- HFGMC Random-C json, NO int failure
% MicroMat{NP + 4} = 196; % -- HFGMC 32x58 Hex, NO int failure
% for NP1 = 1:4
%     NP = NP + 1;
%     OutInfo.Name(NP) = "Section 7.5.1.3 - effect of interface failure";
%     Loads{NP}.Type  = [S,  E,     S,    S,  S,   S];
%     Loads{NP}.Value = [0,  0.01,  0,    0,  0,   0];
%     Loads{NP}.NINC = 200;
%     Loads{NP}.CriteriaOn = [1,0,0,0];
% end


%========================================================================
% -- End of problem definitions
%========================================================================

% -- Store applied loading type identifiers
for I = 1:NP
    Loads{I}.E = E;
    Loads{I}.S = S;
end

end

