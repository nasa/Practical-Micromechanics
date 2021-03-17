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
% Purpose: This script runs multiple stand-alone micromechanics problems, stores the
%          calculated effective properties, and plots a histogram for each property.  
%          The composite materials are defined in GetEffProps.m and the loading and 
%          problem setup is specified in MicroProblemDef.m
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -- Clear memory and close files
clear;
close all;
clc;

% -- Add needed function locations to the path
addpath('Functions/Utilities');
addpath('Functions/WriteResults');
addpath('Functions/Micromechanics');
addpath('Functions/Margins');

%-----------------------------------------------------------------
% I) Define Problems
%-----------------------------------------------------------------
[NProblems, OutInfo, MicroMat, Loads] = MicroProblemDef();

%-----------------------------------------------------------------
% II) Get constituent properties
%-----------------------------------------------------------------
[constitprops] = GetConstitProps;

%-----------------------------------------------------------------
% III) Get effective properties from micromechanics
%-----------------------------------------------------------------
effprops = cell(1,300);

% -- Setup output
[OutInfo] = SetupOutput(OutInfo, 1);

Filename = [OutInfo.OutFile,' - Random Results.txt'];

fid = fopen(Filename,'wt');
Fmf3e2 = '%d \t %E \t %E \t %E \t %E \t %E \t %E \t %E \t %E \t %E \t %E \t %E \t %E \t %E \n';
Fmst = '%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n';
Fmn = '\n';
TSHead = ["NP ,", "Vf ,", "E1 (MPa) ,", "E2 (MPa) ,", "E3 (MPa) ,", "v12 ,", "v13 ,", ...
          "v23 ,", "G12 (MPa) ,", "G13 (MPa) ,", "G23 (MPa) ,", ...
          "alpha1 (1/deg C) ,", "alpha2 (1/deg C) ,", "alpha3 (1/deg C) ,"];

fprintf(fid,Fmst, TSHead);
fprintf(fid,Fmn, ' ');

Vf = zeros(NProblems, 1);
Results = zeros(NProblems, 12);

% -- Loop through problems
for NP = 1: NProblems
    disp(' ');
    disp(['***** Problem Number ', char(num2str(NP)), ' *****']);
    disp(' ');
    mat = MicroMat{NP};
    effprops{mat}.used = true;
    [effprops] = GetEffProps(constitprops, effprops);
    Vf(NP) = effprops{mat}.RUC.Vf;
    Results(NP, 1) = effprops{mat}.E1;
    Results(NP, 2) = effprops{mat}.E2;
    Results(NP, 3) = effprops{mat}.E3;
    Results(NP, 4) = effprops{mat}.v12;
    Results(NP, 5) = effprops{mat}.v13;
    Results(NP, 6) = effprops{mat}.v23;
    Results(NP, 7) = effprops{mat}.G12;
    Results(NP, 8) = effprops{mat}.G13;
    Results(NP, 9) = effprops{mat}.G23;
    Results(NP, 10) = effprops{mat}.a1;
    Results(NP, 11) = effprops{mat}.a2;
    Results(NP, 12) = effprops{mat}.a3;
    
    fprintf(fid, Fmf3e2, NP, Vf(NP), Results(NP, 1),  Results(NP, 2),  Results(NP, 3), ...
                                     Results(NP, 4),  Results(NP, 5),  Results(NP, 6), ...
                                     Results(NP, 7),  Results(NP, 8),  Results(NP, 9), ...
                                     Results(NP, 10), Results(NP, 11), Results(NP, 12));
end

% -- Code to produce histograms from text output file
% -- NOTE: This can be run separately in command window or saved and run as separate
%          script, where 'filename.txt' below is the file containing the data
% Filename = 'filename.txt'; % -- Replace with actual filename
A = importdata(Filename);
Results = A.data(:, 3:14);
TSHead = split(A.textdata, ',');
for i = 1: 12
    figure
    set(gcf,'color','white')
    histogram(Results(:,i));
    txt = char(TSHead(i + 2));
    xlabel(txt(1:end - 1));
    ylabel('frequency');
    Mean = mean(Results(:,i));
    Std = std(Results(:,i));
    Title = ['mean = ', num2str(Mean), ', std = ', num2str(Std)];
    title(Title);
end
