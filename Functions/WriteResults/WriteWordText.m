function WriteWordText(ActX, text, Options) 
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
% Purpose: Write text to word output files
% Input: 
% - ActX: Handle to activeX session for doc output
% - text: Text to be written
% - Options: Struct containing optional word text formatting options
% Output: None
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -- Blank lines before text
if isfield(Options, 'CR_Leading')
    for k = 1: Options.CR_Leading    
        ActX.Selection.TypeParagraph;
    end
end

% -- Text style
if isfield(Options, 'style')
    ActX.Selection.Style = Options.style;
else
    ActX.Selection.Style = 'Normal';
end

% -- Text font
if isfield(Options, 'font')
    ActX.Selection.Font.Name = Options.font;     
end

% -- Text color
if isfield(Options, 'color')
    ActX.Selection.Font.ColorIndex = Options.color;     
end

% -- Write text
ActX.Selection.TypeText(text);

% -- Return to regular color
ActX.Selection.Font.ColorIndex = 'wdAuto';

% -- Blank lines after text
if isfield(Options, 'CR_Trailing')
    for k = 1: Options.CR_Trailing    
        ActX.Selection.TypeParagraph;
    end
end

end