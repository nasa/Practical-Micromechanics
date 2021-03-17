function WriteMicroMoS(mat, OutInfo, Results, Loads, props)
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
% Purpose: Write MoS to output file for stand-alone micomechanics analyses. Also 
%          write stand-alone micomechanics failure envelope output
% Input:
% - mat: Micromechanics material number (within props)
% - OutInfo: Struct containing output information
% - Results: Struct containing micromechanics analysis results
% - Loads: Struct containing problem loading information
% - props: Cell/struct containing composite material & consituent properties 
% Output: None
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Format = OutInfo.Format;

% -- Check that material is micromechanics-based
MicroMats = zeros(1,length(props));
if isfield(props{mat},'micro')
    MicroMats(mat) = 2;
else
    MicroMats(mat) = 1;
end 

if OutInfo.Format == "doc"
    filespec = [pwd,'/', OutInfo.OutFile, '.doc'];

    [fpath,fname,fext] = fileparts(filespec);
    if isempty(fpath); fpath = pwd; end
    if isempty(fext); fext = '.doc'; end
    filespec = fullfile(fpath,[fname,fext]);

    % -- Start an ActiveX session with PowerPoint:
    word = actxserver('Word.Application');
    %word.Visible = 1;
    if ~exist(filespec,'file')
      % -- Create new presentation:
      op = invoke(word.Documents,'Add');
    else
      % -- Open existing presentation:
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

% -- Make sure micro MoS results are present and MoS should be written
if Results.MicroMoS && isfield(OutInfo, 'WriteMargins')
    if OutInfo.WriteMargins

        style='Heading 1';
        text='Margins of Safety';
        WriteText(Format, fid, Fmt1, word, text, style, [0,0])

        WriteAllowables(Format, fid, Fmt1, word, MicroMats, props);

        style='Heading 2';
        text='Micro Scale Margins of Safety';
        WriteText(Format, fid, Fmt1, word, text, style, [1,0])

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
                DataCell = {' ','Subcell or Mat','MoS11','MoS22','MoS33','MoS23','MoS13','MoS12'};
            else
                DataCell = {' ','Subcell or Mat','MoS'};
            end

            if isfield(Results, 'MoS')

                MoS = Results.MoS;

                [DataCell] = ExtractMicroMoS(CN, Results.Type, MoS, props{mat}, DataCell);

            end

            [NoRows,NoCols]=size(DataCell);    
            DataCell1 = DataCell(:,2:NoCols);
            [NoRows,NoCols]=size(DataCell1);    
            if NoRows > 0 && NoCols > 0
                enter = [1 0];
                headerRow = true;
                if Format == "doc" && NoRows > 100
                    disp('** Writing large number of MoS to Word, will take a long time');
                    disp('   Consider killing job and switching to text output');
                end
                WriteTable(Format, fid, headerRow, word, NoRows, NoCols, DataCell1, enter);
            end

        end
        
        % -- Find Controlling MoS
        clear MoS
        MinMoS = 99999;

        if isfield(Results, 'MoS')
            MoS = Results.MoS;

            [MoS] = ControllingMoS(Results.Type, props{mat}, MoS);

             if MoS.MinMoS < MinMoS
                MinMoS = MoS.MinMoS;
             end

             Results.MoS = MoS;

        end

        % -- Write controlling and min MoS info
        MoS = Results.MoS;

        if MoS.MinMoS == 99999
            text='No Controlling subcell or constituent, MinMoS = 99999';
            WriteText(Format, fid, Fmt1, word, text, style, [1,0])

        else

            [CritName, Comp, b, g] = WriteControlling(Format, fid, Fmt1, Fmt2, word, MoS, props{mat});

            style='Heading 2';
            text='Allowable Loads Based on Scaling Prescribed Stresses or Strains:';
            WriteText(Format, fid, Fmt1, word, text, style, [1,0])

            TypeS = ["S11", "S22", "S33", "S23", "S13", "S12"];
            TypeE = ["E11", "E22", "E33", "E23", "E13", "E12"];
            for I = 1:6
                 if Loads.Type(I) == Loads.S
                     Type(I) = TypeS(I);
                 else
                     Type(I) = TypeE(I);
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

        end

        % -- Output for failure envelopes
        if isfield(OutInfo,'EnvFile') && isfield(Loads,'ang')

            if ismember(Loads.E, Loads.Type)
                error('Envelopes only for all Stress loading');
            end
            
            ang = Loads.ang;

            if MoS.MinMoS == 99999
                MoS.ControllingType = 'NA';
                Loads_All = [1E6 1E6 1E6 1E6 1E6 1E6]; 
                Comp = 1;
            else
                if string(MoS.ControllingType) == "FM"
                    text=MoS.Controlling;
                elseif string(MoS.ControllingType) == "bg"
                    text=['(',num2str(b),',',num2str(g),') ',props{mat}.RUC.matsCh(b,g)];
                end                
            end

            text = [num2str(ang),' ',num2str(Loads_All(1)),' ',num2str(Loads_All(2)),' ',num2str(Loads_All(3)), ...
                    ' ',num2str(Loads_All(4)),' ',num2str(Loads_All(5)),' ',num2str(Loads_All(6)),' ', ...
                    ' ',char(char(CritName)),' ',num2str(Comp),' ',text];

            fid1 = fopen(OutInfo.EnvFile,'a');
            style='Normal';
            WriteText('txt', fid1, Fmt1, word, text, style, [1,0])
            fclose(fid1);

        end    
        
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
