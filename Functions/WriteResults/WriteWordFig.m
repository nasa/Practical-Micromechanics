function WriteWordFig(file, SectTitle)
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
% Purpose: Write an existing figure to Word output file
% Input:
% - file: Word file, including path
% - SectTitle: Optional section title to write to word before figure
% Output: None
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -- Copy current fig to clipboard
print -dmeta

% -- Word ActiveX session, catch error
try
    ActX = actxserver('Word.Application');
catch ME
    if (strcmp(ME.identifier,'MATLAB:COM:InvalidProgid'))
        error('*** Cannot start Word in WriteWordFig.m.  Try text output. ***')
    end
end

% -- Open existing or make new word doc
if exist(file, 'file')
    WordDoc = invoke(ActX.Documents, 'Open', file);
else
    WordDoc = invoke(ActX.Documents, 'Add');
end

% -- Find end point of doc
EndPt = get(ActX.activedocument.content, 'end');

% -- Set insert location
set(ActX.application.selection, 'Start', EndPt);
set(ActX.application.selection, 'End', EndPt);

% -- Write optional section title
if nargin == 2
    ActX.Selection.TypeParagraph;
    SectTitle = ['_____________________________________', newline, SectTitle, ' '];
    ActX.Visible = true; 
    ActX.Selection.Font.ColorIndex = 2;
    ActX.Selection.Font.Size = 24;
    ActX.Selection.Font.Italic = 1;
    ActX.Selection.TypeText(SectTitle);
end

% -- Paste
invoke(ActX.selection, 'Paste');

% -- Save word doc
if exist(file, 'file')
  invoke(WordDoc, 'Save');
else
  invoke(WordDoc, 'SaveAs', file, 1);
end

% -- Quit out of Word
invoke(ActX,'Quit');
delete(ActX);


end