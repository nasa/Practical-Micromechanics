function [NP, OutInfo, Geometry, Loads] = LamProblemDef()
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
% Purpose: Sets up and defines classical lamination theory (CLT) problems.  Specifies
%          problem name, ply materials, ply thicknesses, ply orientations, loading, 
%          and failure criteria.
% Input: None
% Output:
% - NP: Number of problems
% - OutInfo: Struct containing output information
% - Geometry: Cell array containing laminate definition variables per
%   problem
% - Loads: Cell array containing problem loading information per problem
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

NProblemsMax = 200; % -- To preallocate cells
NP = 0; % -- Problem number counter
NM = 2; % -- Applied loading type identifier (force/moment resultants)
EK = 1; % -- Applied loading type identifier (midplane strains/curvatures)

% -- Preallocate cells
Geometry = cell(1,NProblemsMax);
Loads = cell(1,NProblemsMax);

% -- Type of output requested
OutInfo.Format = "txt";
% OutInfo.Format = "doc";

% -- Include or exclude plots from output
OutInfo.MakePlots = true;
% OutInfo.MakePlots = false;

% -- Turn on/off writing margins to file (can take long for doc output)
% OutInfo.WriteMargins = true;
OutInfo.WriteMargins = false;

%========================================================================
% -- Define all problems below
%========================================================================

%% -- Chapter 2 Laminate Example Problems

NP = NP + 1;
OutInfo.Name(NP) = "Chapter 2 - Section 2.5.1 glass-epoxy [45]";
Geometry{NP}.Orient = [45];
Geometry{NP}.plymat = [2];
Geometry{NP}.tply = [0.15];
Loads{NP}.Type  = [NM, NM, NM, NM, NM, NM];
Loads{NP}.Value = [ 1,  0,  0,  0,  0,  0];

% NP = NP + 1;
% OutInfo.Name(NP) = "Chapter 2 - Section 2.5.2 IM7-8552 [0,90,90,0] Nx";
% Geometry{NP}.Orient = [0,90,90,0];
% nplies = length(Geometry{NP}.Orient);
% Geometry{NP}.plymat(1:nplies) = 1;
% Geometry{NP}.tply(1:nplies) = 0.15;
% Loads{NP}.Type  = [NM, NM, NM, NM, NM, NM];
% Loads{NP}.Value = [ 1,  0,  0,  0,  0,  0];
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Chapter 2 - Section 2.5.2 IM7-8552 [90,0,0,90] Nx";
% Geometry{NP}.Orient = [90,0,0,90];
% nplies = length(Geometry{NP}.Orient);
% Geometry{NP}.plymat(1:nplies) = 1;
% Geometry{NP}.tply(1:nplies) = 0.15;
% Loads{NP}.Type  = [NM, NM, NM, NM, NM, NM];
% Loads{NP}.Value = [ 1,  0,  0,  0,  0,  0];
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Chapter 2 - Section 2.5.2 IM7-8552 [0,90,90,0] Mx";
% Geometry{NP}.Orient = [0,90,90,0];
% nplies = length(Geometry{NP}.Orient);
% Geometry{NP}.plymat(1:nplies) = 1;
% Geometry{NP}.tply(1:nplies) = 0.15;
% Loads{NP}.Type  = [NM, NM, NM, NM, NM, NM];
% Loads{NP}.Value = [ 0,  0,  0,  1,  0,  0];
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Chapter 2 - Section 2.5.2 IM7-8552 [90,0,0,90] Mx";
% Geometry{NP}.Orient = [90,0,0,90];
% nplies = length(Geometry{NP}.Orient);
% Geometry{NP}.plymat(1:nplies) = 1;
% Geometry{NP}.tply(1:nplies) = 0.15;
% Loads{NP}.Type  = [NM, NM, NM, NM, NM, NM];
% Loads{NP}.Value = [ 0,  0,  0,  1,  0,  0];    
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Chapter 2 - Section 2.5.2 IM7-8552 [0,0,90,90] Nx";
% Geometry{NP}.Orient = [0,0,90,90];
% nplies = length(Geometry{NP}.Orient);
% Geometry{NP}.plymat(1:nplies) = 1;
% Geometry{NP}.tply(1:nplies) = 0.15;
% Loads{NP}.Type  = [NM, NM, NM, NM, NM, NM];
% Loads{NP}.Value = [ 1,  0,  0,  0,  0,  0];    
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Chapter 2 - Section 2.5.2 IM7-8552 [0,0,90,90] Mx";
% Geometry{NP}.Orient = [0,0,90,90];
% nplies = length(Geometry{NP}.Orient);
% Geometry{NP}.plymat(1:nplies) = 1;
% Geometry{NP}.tply(1:nplies) = 0.15;
% Loads{NP}.Type  = [NM, NM, NM, NM, NM, NM];
% Loads{NP}.Value = [ 0,  0,  0,  1,  0,  0];    
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Chapter 2 - Section 2.5.3 IM7-8552 [45,-45,-45,45]";
% Geometry{NP}.Orient = [45,-45,-45,45];
% nplies = length(Geometry{NP}.Orient);
% Geometry{NP}.plymat(1:nplies) = 1;
% Geometry{NP}.tply(1:nplies) = 0.15;
% Loads{NP}.Type  = [NM, NM, NM, NM, NM, NM];
% Loads{NP}.Value = [ 1,  0,  0,  0,  0,  0];
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Chapter 2 - Section 2.5.3 IM7-8552 [90,0,0,90]";
% Geometry{NP}.Orient = [90,0,0,90];
% nplies = length(Geometry{NP}.Orient);
% Geometry{NP}.plymat(1:nplies) = 1;
% Geometry{NP}.tply(1:nplies) = 0.15;
% Loads{NP}.Type  = [NM, NM, NM, NM, NM, NM];
% Loads{NP}.Value = [0.5,0.5,0.5, 0,  0,  0,  0];
% 
%% -- Chapter 3 Laminate Example Problem
%
% NP = NP + 1;
% OutInfo.Name(NP) = "Sect 3.12.4 - cross-ply IM7-8552 MOC";
% Geometry{NP}.Orient = [90,0,0,90];
% nplies = length(Geometry{NP}.Orient);
% Geometry{NP}.plymat(1:nplies) = 103;
% Geometry{NP}.tply(1:nplies) = 0.15;
% Loads{NP}.Type  = [NM,   EK, NM, NM, NM, NM];
% Loads{NP}.Value = [ 0, 0.02,  0,  0,  0,  0];

%% -- Chapter 4 Laminate Example Problems

% NP = NP + 1;
% OutInfo.Name(NP) = "Section 4.4.1.1 - [0] Laminate Nx = 1000";
% Geometry{NP}.Orient = [0];
% nplies = length(Geometry{NP}.Orient);
% Geometry{NP}.plymat(1:nplies) = 1;
% Geometry{NP}.tply(1:nplies) = 1;
% Loads{NP}.DT = 0.;
% Loads{NP}.Type  = [  NM,  NM, NM, NM, NM, NM];
% Loads{NP}.Value = [1000,   0,  0,  0,  0,  0];
% Loads{NP}.CriteriaOn = [1, 1, 1, 1];
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Section 4.4.1.1 - [0] Laminate Ny = 10";
% Geometry{NP}.Orient = [0];
% nplies = length(Geometry{NP}.Orient);
% Geometry{NP}.plymat(1:nplies) = 1;
% Geometry{NP}.tply(1:nplies) = 1;
% Loads{NP}.DT = 0.;
% Loads{NP}.Type  = [  NM,  NM, NM, NM, NM, NM];
% Loads{NP}.Value = [   0,  10,  0,  0,  0,  0];
% Loads{NP}.CriteriaOn = [1, 1, 1, 1];
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Section 4.4.1.1 - [45] Laminate Nx = 1000";
% Geometry{NP}.Orient = [45];
% nplies = length(Geometry{NP}.Orient);
% Geometry{NP}.plymat(1:nplies) = 1;
% Geometry{NP}.tply(1:nplies) = 1;
% Loads{NP}.DT = 0.;
% Loads{NP}.Type  = [  NM,  NM, NM, NM, NM, NM];
% Loads{NP}.Value = [1000,   0,  0,  0,  0,  0];
% Loads{NP}.CriteriaOn = [1, 1, 1, 1];
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Section 4.4.1.1 - [45] Laminate Nx = 104.6";
% Geometry{NP}.Orient = [45];
% nplies = length(Geometry{NP}.Orient);
% Geometry{NP}.plymat(1:nplies) = 1;
% Geometry{NP}.tply(1:nplies) = 1;
% Loads{NP}.DT = 0.;
% Loads{NP}.Type  = [   NM,  NM, NM, NM, NM, NM];
% Loads{NP}.Value = [104.6,   0,  0,  0,  0,  0];
% Loads{NP}.CriteriaOn = [1, 1, 1, 1];
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Section 4.4.1.1 - [45] Laminate Nx = 85.35";
% Geometry{NP}.Orient = [45];
% nplies = length(Geometry{NP}.Orient);
% Geometry{NP}.plymat(1:nplies) = 1;
% Geometry{NP}.tply(1:nplies) = 1;
% Loads{NP}.DT = 0.;
% Loads{NP}.Type  = [   NM,  NM, NM, NM, NM, NM];
% Loads{NP}.Value = [85.35,   0,  0,  0,  0,  0];
% Loads{NP}.CriteriaOn = [1, 1, 1, 1];
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Section 4.4.1.2 - [60,-60,0]s Laminate Nx = 300";
% Geometry{NP}.Orient = [60,-60,0,0,-60,60];
% nplies = length(Geometry{NP}.Orient);
% Geometry{NP}.plymat(1:nplies) = 1;
% Geometry{NP}.tply(1:nplies) = 0.15;
% Loads{NP}.DT = 0.;
% Loads{NP}.Type  = [ NM,  NM, NM, NM, NM, NM];
% Loads{NP}.Value = [300,   0,  0,  0,  0,  0];
% Loads{NP}.CriteriaOn = [1, 1, 1, 1];
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Section 4.4.1.2 - [60,-60,0]s Laminate Ny = 300";
% Geometry{NP}.Orient = [60,-60,0,0,-60,60];
% nplies = length(Geometry{NP}.Orient);
% Geometry{NP}.plymat(1:nplies) = 1;
% Geometry{NP}.tply(1:nplies) = 0.15;
% Loads{NP}.DT = 0.;
% Loads{NP}.Type  = [ NM,  NM, NM, NM, NM, NM];
% Loads{NP}.Value = [  0, 300,  0,  0,  0,  0];
% Loads{NP}.CriteriaOn = [1, 1, 1, 1];
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Section 4.4.1.3 - [60,-60,0]s Laminate Nx = 300, DT = 100";
% Geometry{NP}.Orient = [60,-60,0,0,-60,60];
% nplies = length(Geometry{NP}.Orient);
% Geometry{NP}.plymat(1:nplies) = 1;
% Geometry{NP}.tply(1:nplies) = 0.15;
% Loads{NP}.DT = 100.;
% Loads{NP}.Type  = [ NM,  NM, NM, NM, NM, NM];
% Loads{NP}.Value = [300,   0,  0,  0,  0,  0];
% Loads{NP}.CriteriaOn = [1, 1, 1, 1];
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Section 4.4.1.3 - [60,-60,0]s Laminate Nx = 300, DT = -100";
% Geometry{NP}.Orient = [60,-60,0,0,-60,60];
% nplies = length(Geometry{NP}.Orient);
% Geometry{NP}.plymat(1:nplies) = 1;
% Geometry{NP}.tply(1:nplies) = 0.15;
% Loads{NP}.DT = -100.;
% Loads{NP}.Type  = [ NM,  NM, NM, NM, NM, NM];
% Loads{NP}.Value = [300,   0,  0,  0,  0,  0];
% Loads{NP}.CriteriaOn = [1, 1, 1, 1];
% 
% % -- Laminate Failure Envelopes (Section 4.4.1.4) --
% %    Comment and uncomment to execute this section examples
% tttt = datetime(datetime,'Format','yyyy-MMM-dd HH.mm.ss');
% OutInfo.EnvFile = ['Output/', 'Composite ENVELOPE - ', char(tttt), '.EnvData'];
% ang = 0;
% anginc = 5;
% anginc1 = anginc;
% 
% while (ang <= 360.1)
%     
%     % -- Extra angles for highly anisotropic cases
%     if ang < 9.9 || (ang >= 170 && ang <= 190) || ang >= 350
%       anginc = 0.1;
%     else
%       anginc = anginc1;
%     end
%     
%     NP = NP + 1;
%     OutInfo.Name(NP) = string(['Ang = ', num2str(ang)]);
%     Geometry{NP}.Orient = [0];
%     %Geometry{NP}.Orient = [60,-60,0,0,-60,60];
%     nplies = length(Geometry{NP}.Orient);
%     for k = 1:nplies
%         Geometry{NP}.plymat(k) = 1;
%         Geometry{NP}.tply(k) = 1.;
%         %Geometry{NP}.tply(k) = 0.15;
%     end
%     
%     NY = sin(ang*pi/180);
%     NX = cos(ang*pi/180);
%     Loads{NP}.Type  = [NM, NM, NM, NM, NM, NM];
%     Loads{NP}.Value = [NX, NY,  0,  0,  0,  0];
%     Loads{NP}.DT = 0.0;    
%     %Loads{NP}.DT = 100.0;    
%     %Loads{NP}.DT = -100.0;    
%     % -- Turn criteria on/off to generate all envelopes in Fig. 4.14
%     Loads{NP}.CriteriaOn = [1,1,1,1];
%     %Loads{NP}.CriteriaOn = [1,0,0,0];
%     Loads{NP}.ang = ang;
%     ang = ang + anginc;
% end
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Section 4.4.4 - MT [60,-60,0]s Laminate Nx = 300";
% Geometry{NP}.Orient = [60,-60,0,0,-60,60];
% nplies = length(Geometry{NP}.Orient);
% Geometry{NP}.plymat(1:nplies) = 102;
% Geometry{NP}.tply(1:nplies) = 0.15;
% Loads{NP}.DT = 0.;
% Loads{NP}.Type  = [  NM,  NM, NM, NM, NM, NM];
% Loads{NP}.Value = [ 300,   0,  0,  0,  0,  0];
% Loads{NP}.CriteriaOn = [1, 1, 1, 1];
% 
% NP = NP + 1;
% OutInfo.Name(NP) = "Section 4.4.4 - MOC [60,-60,0]s Laminate Nx = 300";
% Geometry{NP}.Orient = [60,-60,0,0,-60,60];
% nplies = length(Geometry{NP}.Orient);
% Geometry{NP}.plymat(1:nplies) = 103;
% Geometry{NP}.tply(1:nplies) = 0.15;
% Loads{NP}.DT = 0.;
% Loads{NP}.Type  = [  NM,  NM, NM, NM, NM, NM];
% Loads{NP}.Value = [ 300,   0,  0,  0,  0,  0];
% Loads{NP}.CriteriaOn = [1, 1, 1, 1];
%
%% -- Chapter 5 Laminate Example Problems
%      --- None ---
%% -- Chapter 6 Laminate Example Problems
%      --- None ---
%% -- Chapter 7 Laminate Example Problems
% 
% -- 7.1.2 - glass-epoxy quasi-iso
% NP = NP + 1;
% OutInfo.Name(NP) = "Section 7.1.2 - glass-epoxy quasi-iso";
% Geometry{NP}.Orient = [-45,0,45,90,90,45,0,-45];
% nplies = length(Geometry{NP}.Orient);
% Geometry{NP}.plymat(1:nplies) = 171;
% Geometry{NP}.tply(1:nplies)   = 0.15;
% Loads{NP}.Type  = [  EK, NM, NM, NM, NM, NM];
% Loads{NP}.Value = [0.1,  0,  0,  0,  0,  0];
% Loads{NP}.CriteriaOn = [1, 0, 0, 0]; 
% Loads{NP}.PlyDamWatchOn = [1, 1, 1, 1, 0, 0, 0, 0];
% Loads{NP}.NINC = 2000;
% 
% % -- 7.4.2.1 - glass-epoxy [0,90]s
% mat(NP + 1) = 185; % -- HF 7x7
% mat(NP + 2) = 186; % -- HF 32x32
% mat(NP + 3) = 187; % -- MT with HF allowables
% mat(NP + 4) = 188; % -- MT with MT allowables
% 
% for I = 1: 4
%     NP = NP + 1;
%     OutInfo.Name(NP) = "Section 7.4.2.1 - glass-epoxy [0,90]s";
%     Geometry{NP}.Orient = [0,90,90,0];
%     nplies = length(Geometry{NP}.Orient);
%     Geometry{NP}.plymat(1:nplies) = mat(NP);
%     Geometry{NP}.tply(1:nplies)   = 0.25;
%     Loads{NP}.Type  = [  EK, NM, NM, NM, NM, NM];
%     Loads{NP}.Value = [0.04,  0,  0,  0,  0,  0];
%     Loads{NP}.CriteriaOn = [1, 0, 0, 0]; 
%     Loads{NP}.PlyDamWatchOn = [1, 1, 0, 0];
%     Loads{NP}.NINC = 800;
% end
% 
% -- 7.4.2.2 - glass-epoxy [60/-60/0]s, epsx0 applied
% mat(NP + 1) = 185; % -- HF 7x7
% mat(NP + 2) = 186; % -- HF 32x32
% mat(NP + 3) = 187; % -- MT with HF allowables
% mat(NP + 4) = 188; % -- MT with MT allowables
% 
% for I = 1: 4
%     NP = NP + 1;
%     OutInfo.Name(NP) = "Section 7.4.2.2 - glass-epoxy [60,-60,0]s";
%     Geometry{NP}.Orient = [60,-60,0,0,-60,60];
%     nplies = length(Geometry{NP}.Orient);
%     Geometry{NP}.plymat(1:nplies) = mat(NP);
%     Geometry{NP}.tply(1:nplies) = 0.1666666667;
%     Loads{NP}.Type  = [  EK, NM, NM, NM, NM, NM];
%     Loads{NP}.Value = [0.04,  0,  0,  0,  0,  0];
%     Loads{NP}.CriteriaOn = [1, 0, 0, 0]; 
%     Loads{NP}.PlyDamWatchOn = [1, 1, 1, 0, 0, 0];
%     Loads{NP}.NINC = 800;
% end
% 
% % -- 7.4.2.2 - glass-epoxy [60/-60/0]s, epsy0 applied
% mat(NP + 1) = 185; % -- HF 7x7
% mat(NP + 2) = 186; % -- HF 32x32
% 
% for I = 1: 2
%     NP = NP + 1;
%     OutInfo.Name(NP) = "Section 7.4.2.2 - glass-epoxy [60,-60,0]s";
%     Geometry{NP}.Orient = [60,-60,0,0,-60,60];
%     nplies = length(Geometry{NP}.Orient);
%     Geometry{NP}.plymat(1:nplies) = mat(NP);
%     Geometry{NP}.tply(1:nplies) = 0.1666666667;
%     Loads{NP}.Type  = [ NM,   EK, NM, NM, NM, NM];
%     Loads{NP}.Value = [  0, 0.06,  0,  0,  0,  0];
%     Loads{NP}.CriteriaOn = [1, 0, 0, 0]; 
%     Loads{NP}.PlyDamWatchOn = [1, 1, 1, 0, 0, 0];
%     Loads{NP}.NINC = 800;
% end
% 
% % -- 7.4.2.2 - glass-epoxy [45/0/-45/90]s, epsy0 applied
% OutInfo.Name(NP) = "Section 7.4.2.2 - glass-epoxy [45,0,-45,90]s";
% Geometry{NP}.Orient = [45,0,-45,90,90,-45,0,45];
% nplies = length(Geometry{NP}.Orient);
% Geometry{NP}.plymat(1:nplies) = 185; % -- HF 7x7
% Geometry{NP}.tply(1:nplies) = 0.125;
% Loads{NP}.Type  = [ NM,   EK, NM, NM, NM, NM];
% Loads{NP}.Value = [  0, 0.15,  0,  0,  0,  0];
% Loads{NP}.CriteriaOn = [1, 0, 0, 0]; 
% Loads{NP}.PlyDamWatchOn = [1, 1, 1, 1, 0, 0, 0, 0];
% Loads{NP}.NINC = 3000;
% 
% % -- 7.4.2.2 - glass-epoxy [45/0/-45/90]s, epsx0 applied
% mat(NP + 1) = 185; % -- HF 7x7
% mat(NP + 2) = 186; % -- HF 32x32
% mat(NP + 3) = 187; % -- MT with HF allowables
% mat(NP + 4) = 188; % -- MT with MT allowables
% 
% for I = 1: 4
%     NP = NP + 1;
%     OutInfo.Name(NP) = "Section 7.4.2.2 - glass-epoxy [45,0,-45,90]s";
%     Geometry{NP}.Orient = [45,0,-45,90,90,-45,0,45];
%     nplies = length(Geometry{NP}.Orient);
%     Geometry{NP}.plymat(1:nplies) = mat(NP);
%     Geometry{NP}.tply(1:nplies) = 0.125;
%     Loads{NP}.Type  = [  EK, NM, NM, NM, NM, NM];
%     Loads{NP}.Value = [0.15,  0,  0,  0,  0,  0];
%     Loads{NP}.CriteriaOn = [1, 0, 0, 0]; 
%     Loads{NP}.PlyDamWatchOn = [1, 1, 1, 1, 0, 0, 0, 0];
%     Loads{NP}.NINC = 3000;
% end
% 
% % -- 7.4.2.2 - glass-epoxy [45/-45]s, epsx0 applied
% mat(NP + 1) = 185; % -- HF 7x7
% mat(NP + 2) = 186; % -- HF 32x32
% mat(NP + 3) = 187; % -- MT with HF allowables
% mat(NP + 4) = 188; % -- MT with MT allowables
% 
% for I = 1: 4
%     NP = NP + 1;
%     OutInfo.Name(NP) = "Section 7.4.2.2 - glass-epoxy [45,-45]s";
%     Geometry{NP}.Orient = [45,-45,-45,45];
%     nplies = length(Geometry{NP}.Orient);
%     Geometry{NP}.plymat(1:nplies) = mat(NP);
%     Geometry{NP}.tply(1:nplies) = 0.25;
%     Loads{NP}.Type  = [  EK, NM, NM, NM, NM, NM];
%     Loads{NP}.Value = [0.04,  0,  0,  0,  0,  0];
%     Loads{NP}.CriteriaOn = [1, 0, 0, 0]; 
%     Loads{NP}.PlyDamWatchOn = [1, 1, 0, 0];
%     Loads{NP}.NINC = 800;
% end
% 
% % -- Section 7.5.2 - [0/90]s, CMC
% mat(NP + 1) = 197; % -- HF 13x13 Square with int failure
% mat(NP + 2) = 189; % -- HFGMC 32x58 Hex with int failure
% mat(NP + 3) = 198; % -- HFGMC 4-Fiber json with int failure
% mat(NP + 4) = 199; % -- HFGMC 8-Fiber json with int failure
% mat(NP + 5) = 200; % -- HF 13x13 Square without int failure
% mat(NP + 6) = 196; % -- HFGMC 32x58 Hex without int failure
% mat(NP + 7) = 201; % -- HFGMC 4-Fiber json without int failure
% mat(NP + 8) = 202; % -- HFGMC 8-Fiber json without int failure
% 
% for I = 1:8
%     NP = NP + 1;
%     OutInfo.Name(NP) = "Section 7.5.2 - CMC [0,90]s";
%     Geometry{NP}.Orient = [0,90,90,0];
%     nplies = length(Geometry{NP}.Orient);
%     Geometry{NP}.plymat(1:nplies) = mat(NP);
%     Geometry{NP}.tply(1:nplies) = 0.25;
%     Loads{NP}.Type  = [  EK, NM, NM, NM, NM, NM];
%     Loads{NP}.Value = [0.01,  0,  0,  0,  0,  0];
%     Loads{NP}.CriteriaOn = [1, 0, 0, 0]; 
%     Loads{NP}.PlyDamWatchOn = [1, 1, 0, 0];
%     Loads{NP}.NINC = 200;
% end


%========================================================================
% -- End of problem definitions
%========================================================================

% -- Store applied loading type identifiers
for I = 1:NP
    Loads{I}.NM = NM;
    Loads{I}.EK = EK;
end

