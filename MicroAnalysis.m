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
% Purpose: This script runs stand-alone micromechanics problems treating the 
%          composite as a material point.  The composite material is defined in
%          GetEffProps.m and the loading and problem setup is specified in 
%          MicroProblemDef.m
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -- Clear memory and close files
clear;
close all;
fclose('all');
clc;

% -- Add needed function locations to the path
addpath('Functions/Utilities');
addpath('Functions/WriteResults');
addpath('Functions/Micromechanics');
addpath('Functions/Margins');

%-----------------------------------------------------------------
% 1) Define Micromechanics Problems
%-----------------------------------------------------------------
[NProblems, OutInfo, MicroMat, Loads] = MicroProblemDef();

%-----------------------------------------------------------------
% 2) Get constituent properties
%-----------------------------------------------------------------
[constitprops] = GetConstitProps;

%-----------------------------------------------------------------
% 3) Get effective properties from micromechanics
%-----------------------------------------------------------------

% -- Preallocate
effprops = cell(1, 300);
Results = cell(1, NProblems);

% -- Determine which mats are used in any problem
for NP = 1: NProblems
    mat = MicroMat{NP};
    effprops{mat}.used = true;
end

[effprops] = GetEffProps(constitprops, effprops);
    
% -- Loop through problems
for NP = 1: NProblems
    
    mat = MicroMat{NP};
        
    % -- Check for missing problem name
    if (ismissing(OutInfo.Name(NP)))
        OutInfo.Name(NP) = string(['Problem ', char(num2str(NP))]);
    end

    % -- Echo problem info to command window
    disp(['Micro Problem #',num2str(NP),' - ', char(OutInfo.Name(NP))]);
    disp(['   Material Number ',num2str(mat)]);    

    % -- Option to quit if just getting eff props in GetEffProps
    if isfield(effprops{mat}, 'Quit')
        if effprops{mat}.Quit
            disp('   Completed calculation of effective properties -- quitting');
            disp(['  *** Problem ',char(num2str(NP)),' Completed ***'])
            disp(' ');
            continue;
        end
    end

    % -- Check that loads are specified for this problem
    if ~isfield(Loads{NP}, 'Type') ||  ~isfield(Loads{NP}, 'Value')   
        error(strcat('Problem #', num2str(NP), ' Loads not properly defined'));
    end   
    if ~isfield(Loads{NP}, 'DT') % -- Default to zero DT if not specified
        Loads{NP}.DT = 0;
    end
    
    % -- Check that the problem's ply material has been defined
    if ~isfield(effprops{mat}, 'name')
        error(strcat('effective material #', num2str(mat), ...
                     ' undefined ... check GetEffProps'));
    end
            
    %-----------------------------------------------------------------
    % 4) Solve loading and calculate local fields for micromechanics
    %-----------------------------------------------------------------
    
    % -- Calculate global thermal strains
    epsth = Loads{NP}.DT * [effprops{mat}.a1; effprops{mat}.a2; effprops{mat}.a3; ...
                            0; 0; 0];
    % -- Solve for unknown strains and stresses in SG = C*FullGlobalStrain + B
    B = -effprops{mat}.Cstar*epsth;
    [SG, FullGlobalStrain] = SolveLoading(6, effprops{mat}.Cstar, B, Loads{NP});

    % -- Calculate micro scale (constituent level) fields
    [Results{NP}] = MicroFields(FullGlobalStrain, Loads{NP}.DT, effprops{mat});

    %-----------------------------------------------------------------
    % 5) Write micromechanics property and global load results
    %-----------------------------------------------------------------
    effprops{mat}.Mat = mat;
    [OutInfo] = OutputMicro(OutInfo, NP, effprops{mat}, Loads{NP}, SG, ...
                            FullGlobalStrain, epsth);
    
    %-----------------------------------------------------------------
    % 6) Plot micro fields
    %-----------------------------------------------------------------
    if (~isfield(OutInfo,'Format'))
        OutInfo.Format = "txt";
    end
    
    if OutInfo.MakePlots && ~isfield(Loads{NP}, 'ang') % -- Skip plotting for envelopes
        PlotMicroFields(OutInfo, effprops{mat}, Results{NP});
    end
    
    %-----------------------------------------------------------------
    % 7) Calculate Margins
    %-----------------------------------------------------------------
    % -- Check for turning off some failure criteria (Tsai-Wu only by default, 0 = off)
    if (~isfield(Loads{NP},'CriteriaOn'))
        Loads{NP}.CriteriaOn = [0,0,0,1];
    end
    [Results{NP}, effprops{mat}] = CalcMicroMargins(effprops{mat}, Loads{NP}, Results{NP}); 

    %-----------------------------------------------------------------
    % 8) Write Margins
    %-----------------------------------------------------------------
    WriteMicroMoS(mat, OutInfo, Results{NP}, Loads{NP}, effprops)

    %-----------------------------------------------------------------
    % 9) Write RUC (only for RUCid = 300 or 1000)
    %-----------------------------------------------------------------
    WriteRUC(OutInfo, OutInfo.Name(NP), effprops{mat})

  
    disp(['  *** Problem ',char(num2str(NP)),' Completed ***'])
    disp(' ');
    close all;
    
end