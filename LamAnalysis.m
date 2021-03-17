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
% Purpose: Driver script for laminate analysis problems using CLT. The problem input 
%          is defined in the function LamProblemDef.m and GetPlyProps.m
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -- Clear memory and close files
clear;
close all;
fclose('all');
clc;

% -- Add needed function locations to the path
addpath('Functions/CLT');
addpath('Functions/Utilities');
addpath('Functions/WriteResults');
addpath('Functions/Micromechanics');
addpath('Functions/Margins');

%-----------------------------------------------------------------
% 1) Define Laminate Problems
%-----------------------------------------------------------------
[NProblems, OutInfo, Geometry, Loads] = LamProblemDef();

% -- Preallocate LamResults
LamResults = cell(1, NProblems);

%-----------------------------------------------------------------
% 2) Get ply properties
%-----------------------------------------------------------------
plyprops = GetPlyProps(Geometry);

% -- Loop through problems
for NP = 1: NProblems

    % -- Check for missing problem name
    if (ismissing(OutInfo.Name(NP)))
        OutInfo.Name(NP) = string(['Problem ', char(num2str(NP))]);
    end
    
    % -- Check for missing ply angles, thicknesses, materials
    if ~isfield(Geometry{NP},'Orient') || ~isfield(Geometry{NP},'tply') || ...
       ~isfield(Geometry{NP},'plymat')
        error(['Orient, tply, or plymat missing, Problem #', num2str(NP)]);
    end
    
    % -- Echo problem info to command window
    disp(['Problem #',num2str(NP),' - ', char(OutInfo.Name(NP))]);
    disp(['   Per ply angle orientations   [',num2str(Geometry{NP}.Orient), ']']);
    disp(['   Per ply thicknesses          [',num2str(Geometry{NP}.tply), ']']);
    disp(['   Per ply material assignments [',num2str(Geometry{NP}.plymat), ']']);    
    
    % -- Check that the problem's ply materials, orientation, and thickness
    %    have consistent lengths
    if length(Geometry{NP}.tply) ~= length(Geometry{NP}.Orient) || ...
       length(Geometry{NP}.tply) ~= length(Geometry{NP}.plymat) || ...
       length(Geometry{NP}.Orient) ~= length(Geometry{NP}.plymat)
            error('Lengths of tply, Orient, and plymat are not consistent');
    end

    % -- Check that loads are specified for this problem
    if  ~isfield(Loads{NP}, 'Type') || ~isfield(Loads{NP}, 'Value')   
        error(strcat('Problem #', num2str(NP), ' Loads not properly defined'));
    end   
    if ~isfield(Loads{NP}, 'DT') % -- Default to zero DT if not specified
        Loads{NP}.DT = 0;
    end

    % -- Check that the problem's ply materials have been defined
    for k = 1: length(Geometry{NP}.tply)
        if ~isfield(plyprops{Geometry{NP}.plymat(k)}, 'name')
            error(strcat('ply material #', num2str(Geometry{NP}.plymat(k)), ...
                         ' undefined ... check GetPlyProps'));
        end
    end

    % -- Check for slash character in OutInfo.Name
    k = strfind(OutInfo.Name(NP), '/');
    j = strfind(OutInfo.Name(NP), '\');
    if ~isempty(k) || ~isempty(j)
        error('Problem name contains a slash or backslash ... remove');
    end
        
    %-----------------------------------------------------------------
    % 3) Analyze laminate with CLT per problem
    %-----------------------------------------------------------------
    [Geometry{NP}, LamResults{NP}] = CLT(plyprops, Geometry{NP}, Loads{NP});

    %-----------------------------------------------------------------
    % 4) Plot and write laminate properties and stress/strain results
    %-----------------------------------------------------------------
    [OutInfo] = OutputLam(OutInfo, NP, Geometry{NP}, Loads{NP}, plyprops, ...
                          LamResults{NP});

    %-----------------------------------------------------------------
    % 5) Calculate local fields for micromechanics (Chapters 3 - 6)
    %-----------------------------------------------------------------
    [LamResults{NP}] = LamMicroFields(Geometry{NP}, Loads{NP}, plyprops, LamResults{NP});
    
    %-----------------------------------------------------------------
    % 6) Plot micro fields (Chapters 3 - 6)
    %-----------------------------------------------------------------
    if (~isfield(OutInfo,'Format'))
        OutInfo.Format = "txt";
    end
    if ~isfield(Loads{NP}, 'ang') % -- Skip plotting for envelopes
        [plyprops] = OutputMicroFields(OutInfo, Geometry{NP}, plyprops, LamResults{NP});
    end
    
    %-----------------------------------------------------------------
    % 7) Calculate Margins (Chapters 4 - 6)
    %-----------------------------------------------------------------
    % -- Check for turning off some failure criteria (Tsai-Wu only by default, 0 = off)
    if (~isfield(Loads{NP},'CriteriaOn'))
        Loads{NP}.CriteriaOn = [0,0,0,1];
    end
    [LamResults{NP}, plyprops] = CalcLamMargins(plyprops, Geometry{NP}, Loads{NP}, LamResults{NP}); 

    %-----------------------------------------------------------------
    % 8) Write Margins (Chapters 4 - 6)
    %-----------------------------------------------------------------
    WriteLamMoS(OutInfo, LamResults{NP}, Geometry{NP}, Loads{NP}, plyprops);
    
    %-----------------------------------------------------------------
    % 9) Write RUC(s) to .json file (only for RUCid = 300 or 1000)
    %-----------------------------------------------------------------
    WriteRUC(OutInfo, OutInfo.Name(NP), plyprops, Geometry{NP})

    % -- Display notification of file completion to command window
    disp(['  *** Problem ',char(num2str(NP)),' Completed ***'])
    disp(' ');
    close all;
    
end