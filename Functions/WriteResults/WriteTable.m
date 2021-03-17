function WriteTable(Format, fid, headerRow, ActX, NR, NC, Entries, enter)
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
% Purpose: Calls functions to write table to word or text output file
% Input: 
% - Format: 'doc' or 'txt' output file format
% - fid: File id for txt output
% - headerRow: flag to indicate if the first row of table is a header row
% - ActX: Handle to activeX session for doc output
% - NR: Table number of rows
% - NC: Table number of columns
% - Entries: Cell array containing table entries as text
% - enter: Array containing number of blank lines to add before and after text
% Output: None
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if Format == "doc"
    WriteWordTable(ActX, NR, NC, Entries, enter(1), enter(2))
else
    % -- TextTable function is below in this file
    table = TextTable(Entries, headerRow, enter(1), enter(2)); 
    if enter > 0
        Fmt = '\n%s';
    else
        Fmt = '%s';
    end
    fprintf(fid, Fmt, table);
end

return

end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function TextOut = TextTable(CellIn, Header, CR_Leading, CR_Trailing)

% -- Store data in a Cell as a table in text format for output

% -- Cell entries must be text already (not numbers)

% -- Get length of entries
[NR, NC] = size(CellIn);
lengths = zeros(NR, NC);
for IC = 1: NC
    for IR = 1:NR
        lengths(IR, IC) = length(CellIn{IR, IC});
    end
end

% -- Set column length to max of entries (min of 2)
LCol = max(max(lengths,[], 1), 2);

l = newline;

TextOut = '';

% -- Blank lines before text
for k = 1: CR_Leading    
    TextOut = [TextOut, l];
end

% -- Create table
for IR = 1: NR
    Row = '';
    for IC = 1: NC
        entry = CellIn{IR, IC};
        Row = [Row sprintf(' %s %-*s', ' ', LCol(IC), entry)];
    end
    Row = [Row sprintf(' %s', ' ')];

    TextOut = [TextOut  sprintf('%s\n', Row)];
    if Header && IR == 1
        Row = '  ';
        for IC = 1: NC
            entry = repmat('-', 1, LCol(IC) + 2);
            Row = [Row sprintf('%s ',  entry)];
        end
        
        TextOut = [TextOut  sprintf('%s\n', Row)];
    end
end

% -- Blank lines after text
for k = 1: CR_Trailing    
    TextOut = [TextOut, l];
end

end