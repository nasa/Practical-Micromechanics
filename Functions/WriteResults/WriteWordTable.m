function WriteWordTable(ActX, NR, NC, Entries, CR_Leading, CR_Trailing) 
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
% Purpose: Write a table to Word doc output
% Input: 
% - ActX: Handle to activeX session
% - NR: Table number of rows
% - NC: Table number of columns
% - Entries: Cell array containing table entries as text
% - CR_Leading: Number of blank lines to add before table
% - CR_Trailing: Number of blank lines to add after table
% Output: None
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -- Blank lines before table
for k = 1: CR_Leading    
    ActX.Selection.TypeParagraph;
end

% -- Add table
ActX.ActiveDocument.Tables.Add(ActX.Selection.Range, NR, NC, 1, 1);

% -- Write entries in table
for IR= 1: NR
    for IC = 1: NC
        ActX.Selection.TypeText(Entries{IR, IC});
        if IR*IC == NR*NC
            ActX.Selection.MoveDown;
        else 
            ActX.Selection.MoveRight;
        end            
    end
end

% -- Blank lines after table
for k = 1: CR_Trailing    
    ActX.Selection.TypeParagraph;
end

end