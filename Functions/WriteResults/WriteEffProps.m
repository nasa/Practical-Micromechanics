function WriteEffProps(Format, fid, ActX, props)
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
% Purpose: Write micromechanics-based effective properties to doc or txt output file
% Input: 
% - Format: 'doc' or 'txt' output file format
% - fid: File id for txt output
% - ActX: Handle to activeX session for doc output
% - props: Struct containing effective composite properties and constituent info
% Output: None
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Fmt1 = '\n %s';

style='Normal';
text=['Composite Material Number:',' ',char(num2str(props.Mat))];
WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])

style='Normal';
text=['Micromechanics Model:',' ',char(props.micro)];
WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])

style='Normal';
text=['Specified Fiber Volume Fraction:',' ',char(num2str(props.Vf))];
WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])

if isfield(props,'RUC')
    style='Normal';
    text=['Actual Fiber Volume Fraction:',' ',char(num2str(props.RUC.Vf))];
    WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])
end

if isfield(props,'Vi')
    if props.Vi > 0
        style='Normal';
        text=['Specified Interface Volume Fraction:',' ',char(num2str(props.Vi))];
        WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])

        style='Normal';
        text=['Specified Interface Volume Fraction:',' ',char(num2str(props.RUC.Vi))];
        WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])
    end
end

if Format == 'txt'
    WriteText(Format, fid, Fmt1, ActX, ' ', style, [1,0])
end

style='Heading 3';
text='Constituent Material Properties';
WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])

DataCell = {'Mat','Name','EL','ET','NuL','NuT','GL','AlphaL','AlphaT'; ...
            'F',char(props.Fiber.name),num2str(props.Fiber.EL), ...
            num2str(props.Fiber.ET),num2str(props.Fiber.vL), ...
            num2str(props.Fiber.vT),num2str(props.Fiber.GL), ...
            num2str(props.Fiber.aL), num2str(props.Fiber.aT); ...
            'M',char(props.Matrix.name),num2str(props.Matrix.EL), ...
            num2str(props.Matrix.ET), num2str(props.Matrix.vL), ...
            num2str(props.Matrix.vT),num2str(props.Matrix.GL), ...
            num2str(props.Matrix.aL),num2str(props.Matrix.aT)}; 
if isfield(props, 'Interface')
   if isfield(props, 'Vi')
       if props.Vi ~= 0
           Iinfo = {'I',char(props.Interface.name),num2str(props.Interface.EL), ...
                    num2str(props.Interface.ET), num2str(props.Interface.vL), ...
                    num2str(props.Interface.vT),num2str(props.Interface.GL), ...
                    num2str(props.Interface.aL),num2str(props.Interface.aT)};
           DataCell = [DataCell; Iinfo]; 
       end
   end
end
[NoRows,NoCols]=size(DataCell);          
enter = [1 0];
headerRow = true;
WriteTable(Format, fid, headerRow, ActX, NoRows, NoCols, DataCell, enter);

% Write allowables
MicroMats(1) = 2;
Fmt1 = '\r\n %s';
props1{1} = props;
WriteAllowables(Format, fid, Fmt1, ActX, MicroMats, props1);


if isfield(props, 'RUC')
    style='Heading 3';
    text='RUC Information';
    WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])

    if isfield(props, 'Vi')
        Vi = props.Vi;
    else
        Vi = 0;
    end
    DataCell = {'RUCid','Vf','Vi','NB','NG'; ...
                num2str(props.RUC.id),num2str(props.RUC.Vf),num2str(Vi),num2str(props.RUC.NB),num2str(props.RUC.NG)};
    [NoRows,NoCols]=size(DataCell);          
    enter = [1 0];
    headerRow = true;
    WriteTable(Format, fid, headerRow,ActX, NoRows, NoCols, DataCell, enter);

    
    style='Heading 4';
    text='RUC Constituent Arrangement';
    WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])
   
    NB = props.RUC.NB;
    NG = props.RUC.NG;
    style='Normal';
    for IB = 1: NB
        text = '';
        for IG = 1: NG
            text=[text, props.RUC.matsCh(IB,IG), ' '];
        end
        WriteText(Format, fid, Fmt1, ActX, text, style, [1,0], 'Consolas')
    end

    text = ' ';
    WriteText(Format, fid, Fmt1, ActX, text, style, [1,0], 'Consolas')
    
    style='Heading 4';
    text='RUC Subcell dimensions';
    WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])
    
    style='Normal';
    text = 'h =';
    for I = 1: NB
        text=[text, ' ', char(num2str(props.RUC.h(I)))];
    end
    WriteText(Format, fid, Fmt1, ActX, text, style, [1,0], 'Consolas')

    style='Normal';
    text = 'l =';
    for I = 1: NG
        text=[text, ' ', char(num2str(props.RUC.l(I)))];
    end
    WriteText(Format, fid, Fmt1, ActX, text, style, [1,0], 'Consolas')

    text = ' ';
    WriteText(Format, fid, Fmt1, ActX, text, style, [1,0], 'Consolas')
    
end


style='Heading 3';
text='Calculated Effective Mechanical Properties';
WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])

DataCell = {'E1','E2','E3','G12','G13','G23','Nu12','Nu13','Nu23'; ...
            num2str(round(props.E1,2)),num2str(round(props.E2,2)),num2str(round(props.E3,2)), ...
            num2str(round(props.G12,2)),num2str(round(props.G13,2)),num2str(round(props.G23,2)), ...
            num2str(round(props.v12,4)),num2str(round(props.v13,4)),num2str(round(props.v23,4))};
[NoRows,NoCols]=size(DataCell);          
enter = [1 0];
headerRow = true;
WriteTable(Format, fid, headerRow,ActX, NoRows, NoCols, DataCell, enter);

style='Heading 3';
text='Calculated Effective CTEs';
WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])

            
DataCell = {'alpha1','alpha2','alpha3'; ...
            num2str(props.a1),num2str(props.a2),num2str(props.a3)};
        
[NoRows,NoCols]=size(DataCell);          
enter = [1 0];
headerRow = true;
WriteTable(Format, fid, headerRow,ActX, NoRows, NoCols, DataCell, enter);


style='Heading 3';
text='Effective Stiffness Matrix - C*';
WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])

DataCell = cell(6,6);
for i = 1:6
    for j = 1:6
        DataCell{i,j} = num2str(props.Cstar(i,j));
    end
end 

[NoRows,NoCols]=size(DataCell);          
enter = [1 0];
headerRow = false;
WriteTable(Format, fid, headerRow, ActX, NoRows, NoCols, DataCell, enter);

end

