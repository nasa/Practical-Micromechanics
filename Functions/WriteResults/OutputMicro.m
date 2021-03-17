function [OutInfo] = OutputMicro(OutInfo, Problem, props, Loads, SG, ...
                     FullGlobalStrain, epsth, SimDamage)
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
% Purpose: Write output for stand-alone micromechanics problems
% Input:
% - OutInfo: Struct containing output information
% - Problem: Current problem number
% - props: Struct containing effective composite properties
% - Loads: Struct containing problem loading information
% - SG: Vector of global stress components
% - FullGlobalStrain: Vector of global total strain components
% - epsth: Vector of thermal strain components
% - SimDamage: Optional argument, flag for progressive damage problem (Ch7)
% Output:
% - OutInfo: Updated struct containing output information
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if nargin == 7
    SimDamage = false;
end

[OutInfo] = SetupOutput(OutInfo, Problem);

Format = OutInfo.Format;
Fmt = '%s \n';

CurDir = pwd;
if Format == "doc"
    file = fullfile(CurDir, OutInfo.OutFileName);
    fid = 0;

    % -- Word ActiveX session, catch error
    try
        ActX = actxserver('Word.Application');
    catch ME
        if (strcmp(ME.identifier,'MATLAB:COM:InvalidProgid'))
            error('*** Cannot start Word in WriteWordFig.m.  Try text output. ***')
        end
    end
    
    % ActX.Visible = false;
    ActX.Visible = true; % -- This allows you to see doc being created

    % -- Open existing or make new word doc
    if exist(file, 'file')
        WordDoc = invoke(ActX.Documents, 'Open', file);
    else
        WordDoc = invoke(ActX.Documents, 'Add');
    end

else
    fid = fopen(OutInfo.OutFileName,'wt');
    ActX = 0;
    
end

style='Title';
text='Composite Micromechanics Report';
WriteText(Format, fid, Fmt, ActX, text, style, [0,1])

text=char(OutInfo.Name(Problem));
WriteText(Format, fid, Fmt, ActX, text, style, [0,1])

style='Normal';
text=string(OutInfo.datetime);
WriteText(Format, fid, Fmt, ActX, text, style, [0,1])

WriteEffProps(OutInfo.Format, fid, ActX, props);

Fmt1 = '\n %s';
style='Heading 1';
text='Applied Loading (shear strains are engr)';
WriteText(Format, fid, Fmt1, ActX, text, style, [0,0])

comp = ['11'; '22'; '33'; '23'; '13'; '12'];
for I = 1:6
     if Loads.Type(I) == Loads.S
         Type(I) = string(['S',comp(I,1:2)]);
     else
         Type(I) = string(['E',comp(I,1:2)]);
     end
end

if SimDamage % -- For Chapter 7 (prog dam)
    Value = Loads.Final;
    DT = Loads.DTFinal;
else
    Value = Loads.Value;
    DT = Loads.DT;
end

DataCell = {char(Type(1)),char(Type(2)),char(Type(3)),char(Type(4)),char(Type(5)),char(Type(6)),'DT'; ...
            num2str(Value(1)),num2str(Value(2)),num2str(Value(3)), ...
            num2str(Value(4)),num2str(Value(5)),num2str(Value(6)),num2str(DT)};

[NoRows,NoCols]=size(DataCell);          
enter = [1 0];
headerRow = true;
WriteTable(Format, fid, headerRow, ActX, NoRows, NoCols, DataCell, enter);

% -- Used only for Chapter 7 (prog dam)
if SimDamage

    Fmt1 = '\n %s';
    style='Heading 1';
    text='Number of Loading Increments';
    WriteText(Format, fid, Fmt1, ActX, text, style, [0,0])
    
    style='Normal';
    text=['  ', char(num2str(Loads.NINC))];
    WriteText(Format, fid, Fmt1, ActX, text, style, [0,0])
    text='  ';
    WriteText(Format, fid, Fmt1, ActX, text, style, [0,0])
    
    Fmt1 = '\n %s';
    style='Heading 1';
    text='Active Failure Criteria';
    WriteText(Format, fid, Fmt1, ActX, text, style, [0,0])
    
    OnOff = {'Off'; 'Off'; 'Off'; 'Off'};
    if ~isfield(Loads, 'CriteriaOn')
        OnOff = {'On'; 'On'; 'On'; 'On'};
    else
        for I = 1:4
            if Loads.CriteriaOn(I) == 1
                OnOff{I} = 'On';
            end
        end
    end
    
    DataCell = {'Max Stress', 'Max Strain', 'Tsai-Hill', 'Tsai-Wu'; ...
                OnOff{1}, OnOff{2}, OnOff{3}, OnOff{4}};
            
    [NoRows,NoCols]=size(DataCell);          
    enter = [2 0];
    headerRow = true;
    WriteTable(Format, fid, headerRow, ActX, NoRows, NoCols, DataCell, enter);

    style='Normal';
    text='  ';
    WriteText(Format, fid, Fmt1, ActX, text, style, [0,0])
end

if ~SimDamage
    style='Heading 2';
    text='Resulting Global Stresses';
    WriteText(Format, fid, Fmt1, ActX, text, style, [0,0])

    for I = 1:6
        Type(I) = string(['S',comp(I,1:2)]);
    end
    DataCell = {char(Type(1)),char(Type(2)),char(Type(3)),char(Type(4)),char(Type(5)),char(Type(6)); ...
                 num2str(SG(1)),num2str(SG(2)),num2str(SG(3)), ...
                 num2str(SG(4)),num2str(SG(5)),num2str(SG(6))};

    [NoRows,NoCols]=size(DataCell);          
    enter = [1 0];
    headerRow = true;
    WriteTable(Format, fid, headerRow, ActX, NoRows, NoCols, DataCell, enter);

    style='Heading 2';
    text='Resulting Global Strains (shear strains are engr)';
    WriteText(Format, fid, Fmt1, ActX, text, style, [0,0])

    for I = 1:6
        Type(I) = string(['E',comp(I,1:2)]);
    end
    DataCell = {char(Type(1)),char(Type(2)),char(Type(3)),char(Type(4)),char(Type(5)),char(Type(6)); ...
                 num2str(FullGlobalStrain(1)),num2str(FullGlobalStrain(2)),num2str(FullGlobalStrain(3)), ...
                 num2str(FullGlobalStrain(4)),num2str(FullGlobalStrain(5)),num2str(FullGlobalStrain(6))};

    [NoRows,NoCols]=size(DataCell);          
    enter = [1 0];
    headerRow = true;
    WriteTable(Format, fid, headerRow, ActX, NoRows, NoCols, DataCell, enter);

    style='Heading 2';
    text='Thermal Strains (shear strains are engr)';
    WriteText(Format, fid, Fmt1, ActX, text, style, [0,0])

    for I = 1:6
        Type(I) = string(['Epsth',comp(I,1:2)]);
    end
    DataCell = {char(Type(1)),char(Type(2)),char(Type(3)),char(Type(4)),char(Type(5)),char(Type(6)); ...
                 num2str(epsth(1)),num2str(epsth(2)),num2str(epsth(3)), ...
                 num2str(epsth(4)),num2str(epsth(5)),num2str(epsth(6))};

    [NoRows,NoCols]=size(DataCell);          
    enter = [1 0];
    headerRow = true;
    WriteTable(Format, fid, headerRow, ActX, NoRows, NoCols, DataCell, enter);
end

if Format == "doc"

    % -- Save word doc
    if exist(file, 'file')
      invoke(WordDoc, 'Save');
    else
      invoke(WordDoc, 'SaveAs', file, 1);
    end

    % -- Quit out of Word
    invoke(ActX,'Quit');
    delete(ActX);  
    
else
    fclose(fid);
    
end

end