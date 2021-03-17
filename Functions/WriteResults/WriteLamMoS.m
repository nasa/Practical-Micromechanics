function WriteLamMoS(OutInfo, LamResults, Geometry, Loads, plyprops)
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
% Purpose: Write MoS to output file for laminate-based analyses.  Also write 
%          laminate-based failure envelope output
% Input:
% - OutInfo: Struct containing output information
% - LamResults: Struct containing laminate analysis results
% - Geometry: Struct containing laminate definition variables
% - Loads: Struct containing problem loading information
% - plyprops: Cell/struct containing effective composite material properties 
% Output: None
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Format = OutInfo.Format;
LayerZ = Geometry.LayerZ;
Orient = Geometry.Orient;
N = Geometry.N;
NM = LamResults.NM;

% -- Identify constituent materials existing in laminate
MicroMats = zeros(1,length(plyprops));
for k = 1:Geometry.N
    mat = Geometry.plymat(k);
    if isfield(plyprops{mat},'micro')
        MicroMats(mat) = 2;
    else
        MicroMats(mat) = 1;
    end 
end 
    
if OutInfo.Format == "doc"
    filespec = [pwd,'/', OutInfo.OutFile, '.doc'];

    [fpath,fname,fext] = fileparts(filespec);
    if isempty(fpath); fpath = pwd; end
    if isempty(fext); fext = '.doc'; end
    filespec = fullfile(fpath,[fname,fext]);

    % -- Start an ActiveX session with Word:
    word = actxserver('Word.Application');
    %word.Visible = 1;
    if ~exist(filespec,'file')
       % -- Create new doc:
      op = invoke(word.Documents,'Add');
    else
       % -- Open existing doc:
      op = invoke(word.Documents,'Open',filespec);
    end
    % -- Find end of document and make it the insertion point:
    end_of_doc = get(word.activedocument.content,'end');
    set(word.application.selection,'Start',end_of_doc);
    set(word.application.selection,'End',end_of_doc);
    word.Selection.InsertBreak; %pagebreak
    fid = 0;

elseif OutInfo.Format == "txt"
    fid = fopen([OutInfo.OutFile, '.txt'],'a');    
    word = 0;

end

% -- Text formats
Fmt = '%s \r\n';
Fmt1 = '\r\n %s';
Fmt2 = '\r\n %s \r\n';

% ----------------
% -- Ply level MoS
% ----------------
if ~isfield(OutInfo, 'WriteMargins')
   OutInfo.WriteMargins = false; 
end

if isfield(LamResults,'MoS') && OutInfo.WriteMargins

    MoS = LamResults.MoS;
    kk_ctrl = LamResults.MoSControllingZPt;

    style='Heading 1';
    text='Ply Level Margins of Safety';
    WriteText(Format, fid, Fmt1, word, text, style, [1,0])

    style='Heading 2';
    text='Specified Ply Stress and Strain Allowables';
    WriteText(Format, fid, Fmt1, word, text, style, [1,0])

    for mat = 1:length(plyprops)
        if MicroMats(mat) ~= 0
            style='Heading 3';
            text=['Material Number',' ',num2str(mat),' - ',char(plyprops{mat}.name)];
            WriteText(Format, fid, Fmt1, word, text, style, [1,0])
        end

        if MicroMats(mat) == 1 
            DataCell = {'XT','XC','YT','YC','S','XeT','XeC','YeT','YeC','Se'; ...
                        num2str(plyprops{mat}.allowables.XT),num2str(plyprops{mat}.allowables.XC),num2str(plyprops{mat}.allowables.YT), ...
                        num2str(plyprops{mat}.allowables.YC),num2str(plyprops{mat}.allowables.S), ...
                        num2str(plyprops{mat}.allowables.XeT),num2str(plyprops{mat}.allowables.XeC),num2str(plyprops{mat}.allowables.YeT), ...
                        num2str(plyprops{mat}.allowables.YeC),num2str(plyprops{mat}.allowables.Se)};
            [NoRows,NoCols]=size(DataCell);          
            enter = [1 0];
            headerRow = true;
            WriteTable(Format, fid, headerRow, word, NoRows, NoCols, DataCell, enter);
        end

    end

    style='Heading 2';
    text='Ply Margins of Safety through the thickness';
    WriteText(Format, fid, Fmt1, word, text, style, [1,0])

    DataCell = cell(2*N,5);
    DataCell{1,1} = 'Z-location';
    DataCell{1,2} = 'Angle';
    DataCell{1,3} = 'Criterion';
    DataCell{1,4} = 'MoS or MoS11';
    DataCell{1,5} = 'MoS22';
    DataCell{1,6} = 'MoS12';
    row = 1;
    for kk = 1:2*N
        DataCell{row+1, 1} = num2str(LayerZ(kk));
        DataCell{row+2, 1} = ' ';
        DataCell{row+3, 1} = ' ';
        DataCell{row+4, 1} = ' ';
        DataCell{row+1, 2} = num2str(Orient(round(kk/2)));
        DataCell{row+2, 2} = ' ';
        DataCell{row+3, 2} = ' ';
        DataCell{row+4, 2} = ' ';
        for i = 1:4
            row = row + 1;
            DataCell{row, 3} = char(MoS(kk).Crit{i}.name);
            if i < 3
                if Loads.CriteriaOn(i) == 0
                    DataCell{row, 4} = 'off';
                    DataCell{row, 5} = 'off';
                    DataCell{row, 6} = 'off';
                else
                    DataCell{row, 4} = num2str(MoS(kk).Crit{i}.MoS(1));
                    DataCell{row, 5} = num2str(MoS(kk).Crit{i}.MoS(2));
                    DataCell{row, 6} = num2str(MoS(kk).Crit{i}.MoS(6));
                end
            else
                if Loads.CriteriaOn(i) == 0
                    DataCell{row, 4} = 'off';
                    DataCell{row, 5} = ' ';
                    DataCell{row, 6} = ' ';
                else
                    DataCell{row, 4} = num2str(MoS(kk).Crit{i}.MoS);
                    DataCell{row, 5} = ' ';
                    DataCell{row, 6} = ' ';
                end
            end
        end
    end 

    [NoRows,NoCols]=size(DataCell);          
    enter = [1 0];
    headerRow = true;
    WriteTable(Format, fid, headerRow, word, NoRows, NoCols, DataCell, enter);

    if kk_ctrl > 0
        k = round(kk_ctrl/2);
        style='Normal';
        text=['Controlling Z-location is ', num2str(LayerZ(kk_ctrl)), ' within Ply # ',num2str(k), ' (Angle = ', num2str(Orient(k)),')'];
        WriteText(Format, fid, Fmt1, word, text, style, [1,0])

        CritNum = MoS(kk_ctrl).ControllingNum;
        if CritNum < 3
            text=['Controlling Criterion is ', char(MoS(kk_ctrl).Controlling), ' Component ',num2str(MoS(kk_ctrl).Crit{CritNum}.Controlling)];
        else
            text=['Controlling Criterion is ', char(MoS(kk_ctrl).Controlling)];
        end
        style='Normal';
        WriteText(Format, fid, Fmt2, word, text, style, [1,0])

        CritNum = MoS(kk_ctrl).ControllingNum;
        style='Normal';
        text=['Minimum MoS: ', num2str(MoS(kk_ctrl).MinMoS)];
        WriteText(Format, fid, Fmt2, word, text, style, [1,0])

        style='Heading 2';
        text='Allowable Loads Based on Scaling Prescribed Loads:';
        WriteText(Format, fid, Fmt1, word, text, style, [1,0])

        TypeNM = ["NXX", "NYY", "NXY", "MXX", "MYY", "MXY"];
        TypeEK = ["EXX", "EYY", "EXY", "KXX", "KYY", "KXY"];
        for I = 1:6
             if Loads.Type(I) == Loads.NM
                 Type(I) = TypeNM(I);
             else
                 Type(I) = TypeEK(I);
             end
         end
                
         Loads_All = Loads.Value*(MoS(kk_ctrl).MinMoS + 1);
         DataCell = {char(Type(1)),char(Type(2)),char(Type(3)),char(Type(4)),char(Type(5)),char(Type(6)); ...
         num2str(Loads_All(1)),num2str(Loads_All(2)),num2str(Loads_All(3)), ...
         num2str(Loads_All(4)),num2str(Loads_All(5)),num2str(Loads_All(6))};
        
        % -- Output for ply-level failure envelopes
        if isfield(OutInfo,'EnvFile') && isfield(Loads, 'ang')

            ang = Loads.ang;

            text = [num2str(ang),' ',num2str(Loads_All(1)),' ',num2str(Loads_All(2)),' ',num2str(Loads_All(3)), ...
                    ' ',num2str(Loads_All(4)),' ',num2str(Loads_All(5)),' ',num2str(Loads_All(6)),' ',num2str(Orient(k)),'deg',' ',char(MoS(kk_ctrl).Controlling),' ', ...
                    char(num2str(MoS(kk_ctrl).Crit{CritNum}.Controlling))];
            fid1 = fopen(OutInfo.EnvFile,'a');
            style='Normal';
            WriteText('txt', fid1, Fmt1, word, text, style, [1,0])
            fclose(fid1);
            
        end
        
        [NoRows,NoCols]=size(DataCell);          
        enter = [1 0];
        headerRow = true;
        WriteTable(Format, fid, headerRow, word, NoRows, NoCols, DataCell, enter);

    else
        style='Normal';
        text='No Controlling Ply Detected';
        WriteText(Format, fid, Fmt1, word, text, style, [1,0])

    end

    style='Heading 2';
    text=' ';
    WriteText(Format, fid, Fmt1, word, text, style, [1,0])

end

% ---------------------------
% -- Micromechanics-based MoS
% ---------------------------
if LamResults.MicroMoS && OutInfo.WriteMargins

    style='Heading 1';
    text='Margins of Safety';
    WriteText(Format, fid, Fmt1, word, text, style, [0,0])

    WriteAllowables(Format, fid, Fmt1, word, MicroMats, plyprops);

    style='Heading 2';
    text='Micro Scale Margins of Safety';
    WriteText(Format, fid, Fmt1, word, text, style, [1,0])

    N = Geometry.N;    
    LayerZ = Geometry.LayerZ;        
    for CN =1:4

        if Loads.CriteriaOn(CN) == 0
            continue
        end
        
        if CN == 1
            CritName = 'Max Stress';
        elseif CN == 2
            CritName = 'Max Strain';
        elseif CN == 3
            CritName = 'Tsai-Hill';
        elseif CN == 4
            CritName = 'Tsai-Wu';
        end

        style='Heading 3';
        text=CritName;
        WriteText(Format, fid, Fmt1, word, text, style, [1,0])

        if CN < 3
            DataCell = {'Z-location','Subcell or Mat','MoS11','MoS22','MoS33','MoS23','MoS13','MoS12'};
        else
            DataCell = {'Z-location','Subcell or Mat','MoS'};
        end

        for kk =1:2*N

           if isfield(LamResults.Micro(kk), 'MoS')

               MoS(kk) = LamResults.Micro(kk).MoS;

               k = round(kk/2);
               mat = Geometry.plymat(k);    

               [DataCell] = ExtractMicroMoS(CN, LamResults.Micro(kk).Type, ...
                            MoS(kk), plyprops{mat}, DataCell, LayerZ(kk));
               
           end

            [NoRows,NoCols]=size(DataCell);    
            if NoRows > 0 && NoCols > 0
                enter = [1 0];
                headerRow = true;
                WriteTable(Format, fid, headerRow, word, NoRows, NoCols, DataCell, enter);
            end
        
        end
    end
    
    % -- Find Controlling MoS
    clear MoS
    kk_ctrl = 0;
    MinMoS = 99999;
    for kk = 1:2*Geometry.N

        if isfield(LamResults.Micro(kk), 'MoS')
            k = round(kk/2);
            mat = Geometry.plymat(k);    
            MoS = LamResults.Micro(kk).MoS;

            [MoS] = ControllingMoS(LamResults.Micro(kk).Type, plyprops{mat}, MoS);

             if MoS.MinMoS < MinMoS
                MinMoS = MoS.MinMoS;
                kk_ctrl = kk;
             end
            
             LamResults.Micro(kk).MoS = MoS;
             
        end
    end

    % -- Write controlling and min MoS info
    if kk_ctrl > 0
        MoS = LamResults.Micro(kk_ctrl).MoS;
        k = round(kk_ctrl/2);
        mat =  Geometry.plymat(k);
        style='Normal';
        text=['Controlling Z-location is ', num2str(LayerZ(kk_ctrl)), ' within Ply # ',num2str(k), ' (Angle = ', num2str(Geometry.Orient(k)),')'];
        WriteText(Format, fid, Fmt1, word, text, style, [1,0])

        [CritName, Comp, b, g] = WriteControlling(Format, fid, Fmt1, Fmt2, word, MoS, plyprops{mat});

        style='Heading 2';
        text='Allowable Loads Based on Scaling Prescribed Loads:';
        WriteText(Format, fid, Fmt1, word, text, style, [1,0])

        
        TypeNM = ["NXX", "NYY", "NXY", "MXX", "MYY", "MXY"];
        TypeEK = ["EXX", "EYY", "EXY", "KXX", "KYY", "KXY"];
        for I = 1:6
             if Loads.Type(I) == Loads.NM
                 Type(I) = TypeNM(I);
             else
                 Type(I) = TypeEK(I);
             end
         end

         Loads_All = Loads.Value*(MoS.MinMoS + 1);
         DataCell = {char(Type(1)),char(Type(2)),char(Type(3)),char(Type(4)),char(Type(5)),char(Type(6)); ...
         num2str(Loads_All(1)),num2str(Loads_All(2)),num2str(Loads_All(3)), ...
         num2str(Loads_All(4)),num2str(Loads_All(5)),num2str(Loads_All(6))};

        [NoRows,NoCols]=size(DataCell);          
        enter = [1 0];
        headerRow = true;
        WriteTable(Format, fid, headerRow, word, NoRows, NoCols, DataCell, enter);

        % -- Output for micromechanics-based failure envelopes
        if isfield(OutInfo,'EnvFile') && isfield(Loads,'ang')
                       
            if ismember("EK", Type)
                error('Envelopes only available for all NM loading');
            end
            
            ang = Loads.ang;
            
            if string(MoS.ControllingType) == "FM"
                text=MoS.Controlling;
            elseif string(MoS.ControllingType) == "bg"
                text=['(',num2str(b),',',num2str(g),') ',plyprops{mat}.RUC.matsCh(b,g)];
            end

            text = [num2str(ang),' ',num2str(Loads_All(1)),' ',num2str(Loads_All(2)),' ',num2str(Loads_All(3)), ...
                    ' ',num2str(Loads_All(4)),' ',num2str(Loads_All(5)),' ',num2str(Loads_All(6)),' ', ...
                    num2str(Orient(k)),'deg',' ',char(char(CritName)),' ',num2str(Comp),' ',text];
                
            fid1 = fopen(OutInfo.EnvFile,'a');
            style='Normal';
            WriteText('txt', fid1, Fmt1, word, text, style, [1,0])
            fclose(fid1);
            
        end
        
        
    else
        style='Normal';
        text='No Controlling Ply Detected';
        WriteText(Format, fid, Fmt1, word, text, style, [1,0])

    end    

end

if OutInfo.Format == "doc"
    % -- Close the word window:
    invoke(op,'Close');
    % -- Quit MS Word
    invoke(word,'Quit');
    % -- Close word and terminate ActiveX:
    delete(word);  
elseif OutInfo.Format == "txt"
    fclose(fid);
end
        
return  
