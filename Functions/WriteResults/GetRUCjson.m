function [json] = GetRUCjson(nRUC, RUC)

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
% Purpose: Write RUC information to a JSON file so a random RUC can be used in the 
%          future.
% Input:
% - nRUC: Number of RUCs in RUC
% - RUC: Cell/struct containing RUC information
% Output:
% - json: text of RUC information in json format
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

l = newline;
cl = [',', l];
t = sprintf('\t');

json = ['{', l, t, '"RUC": [', l];

for I = 1: nRUC
    
    hstr = sprintf('%g,' , RUC{I}.h);
    hstr = hstr(1: end - 1);
    lstr = sprintf('%g,' , RUC{I}.h);
    lstr = lstr(1: end - 1);
    
    mats = '[';
    matsCh = '[';
    for IB = 1: RUC{I}.NB
        matl = sprintf('%i,' , RUC{I}.mats(IB,:));
        matl = matl(1: end - 1);
        mats = [mats, l, t, t, t, t, '[', matl, '],'];
        matsCh = [matsCh, l, t, t, t, t, '"', RUC{I}.matsCh(IB,:), '",'];
    end
    
    % -- Remove last comma on mats and matsCh
    mats = mats(1: end - 1);
    matsCh = matsCh(1: end - 1);
    
    json = [json, t, t, '{', l, ...
            t, t, t, '"id": ', num2str(RUC{I}.id), cl, ...
            t, t, t, '"Mat": ', num2str(RUC{I}.Mat), cl, ...
            t, t, t, '"Vf": ', num2str(RUC{I}.Vf, '%g'), cl, ...
            t, t, t, '"Vf_max": ', num2str(RUC{I}.Vf_max), cl, ...
            t, t, t, '"Vi": ', num2str(RUC{I}.Vi), cl, ...
            t, t, t, '"NB": ', num2str(RUC{I}.NB), cl, ...
            t, t, t, '"NG": ', num2str(RUC{I}.NG), cl, ...
            t, t, t, '"h": ', '[', hstr, ']', cl, ...
            t, t, t, '"l": ', '[', lstr, ']', cl, ...
            t, t, t, '"mats": ', mats, l, t, t, t, ']', cl, ...
            t, t, t, '"matsCh": ', matsCh, l, t, t, t, ']', l];
        
    json = [json, t, t, '}'];

    if I < nRUC
        json = [json, cl];
    end
end

json = [json, l, t, ']', l '}', l];

end
