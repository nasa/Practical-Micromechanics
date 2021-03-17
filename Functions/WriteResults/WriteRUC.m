function WriteRUC(OutInfo, ProblemName, props_in, Geometry)
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
% - OutInfo: Struct containing output information
% - ProblemName: Name of the problem
% - props: Cell/struct containing effective composite properties
% - Geometry: Cell array containing laminate definition variables (optional)
% Output: None
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if iscell(props_in)
    props = props_in;
else
    props{1} = props_in;
end


ic = 0;
for plymat = 1:length(props)
    if isfield(props{plymat}, 'RUC')
        %props{plymat}.RUC.plymat = plymat;
        if isfield(props{plymat}, 'Mat')
            props{plymat}.RUC.Mat = props{plymat}.Mat;
        else
            props{plymat}.RUC.Mat = plymat;
        end
        % -- If laminate, skip if material is not in laminate
        if nargin == 4 && ~ismember(plymat, Geometry.plymat)
            continue
        end
        % -- Only write JSON for random and JSON RUCs
        if (props{plymat}.RUC.id == 300) || (props{plymat}.RUC.id == 1000)
            ic = ic + 1;
            RUCjson{ic} = props{plymat}.RUC;
        end
    end
end

if exist('RUCjson','var')
    RUCFileName = [OutInfo.Path, 'RUCs - ', char(ProblemName), ' - ', OutInfo.datetime, '.json'];
    fid = fopen(RUCFileName, 'w');
    [jsondata] = GetRUCjson(length(RUCjson), RUCjson);
    fprintf(fid, jsondata);
    fclose(fid);
    % -- If prog damage problem, write json to snapshot directory
    if isfield(OutInfo, 'DamDir')
        RUCFileName = [OutInfo.DamDir, 'RUCs - ', char(ProblemName), ' - ', OutInfo.datetime, '.json'];
        fid = fopen(RUCFileName, 'w');
        [jsondata] = GetRUCjson(length(RUCjson), RUCjson);
        fprintf(fid, jsondata);
        fclose(fid);
    end
end


end
