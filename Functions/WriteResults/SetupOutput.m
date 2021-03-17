function [OutInfo] = SetupOutput(OutInfo, Problem)
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
% Purpose: Setup output filenames and paths
% Input: 
% - OutInfo: Struct containing output information
% - Problem: Current problem number
% Output:
% - OutInfo: Updated struct containing output information
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -- Output files tagged with date and time
tttt = datetime(datetime,'Format','yyyy-MMM-dd HH.mm.ss');
[~,~,~] = mkdir('Output');

disp(tttt);

% -- MS Word output
if OutInfo.Format == "doc"
   OutFile = ['Output/', 'Composite Word Report - ',char(OutInfo.Name(Problem)),' ',char(tttt),'-',char(num2str(Problem))];
   OutFileName = [OutFile, '.doc'];
   PlotFile = ' '; % Not used for Word output
   
% -- Text file output   
elseif OutInfo.Format == "txt"
   OutFile = ['Output/', 'Composite Text Report - ',char(OutInfo.Name(Problem)),' ',char(tttt),'-',char(num2str(Problem))];
   OutFileName = [OutFile, '.txt'];
   PlotFile = ['Output/', 'Composite PLOT - ',char(OutInfo.Name(Problem)),' ',char(tttt)];
end

% -- Store output information
OutInfo.PlotFile = PlotFile;
OutInfo.OutFile = OutFile;
OutInfo.OutFileName = OutFileName;
OutInfo.datetime = char(tttt);
OutInfo.Path = 'Output/';

end
