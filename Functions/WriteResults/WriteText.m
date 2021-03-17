function WriteText(Format, fid, Fmt, ActX, text, style, enters, font, color)
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
% Purpose: Write text to word or text output files
% Input: 
% - Format: 'doc' or 'txt' output file format
% - fid: File id for txt output
% - Fmt: Format specifier for txt output
% - ActX: Handle to activeX session for doc output
% - text: Text to be written
% - style: Text style for doc output
% - enters: Array containing number of blank lines to add before and after text
% - font: Text font for doc output
% - color: Text color for doc output
% Output: None
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if Format == "doc"
    Options.CR_Leading = enters(1);
    Options.CR_Trailing = enters(2);
    Options.style = style;
    if nargin == 8
        Options.font = font;
    elseif nargin == 9
        Options.color = color;
    end    
    WriteWordText(ActX, text, Options)
    
else
    if string(style) == "Heading 1"
        text = ['<==== ',text,' ====>'];
    elseif string(style) == "Heading 2"
        text = ['-- ',text,' --'];
    elseif string(style) == "Heading 3"
        text = [' ',text];
    elseif string(style) == "Normal"
        text = ['   ',text];
    end
    fprintf(fid, Fmt, text);
end

return

