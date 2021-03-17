function WriteLamOutput(OutInfo, Problem, Geometry, Loads, LamResults, plyprops)
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
% Purpose: Write output for laminate problems
% Input: 
% - OutInfo: Struct containing output information
% - Problem: Current problem number
% - Geometry: Struct containing laminate definition variables
% - Loads: Struct containing problem loading information
% - LamResults: Struct containing laminate analysis results
% - plyprops: Cell array containing ply material properties
% Output: None
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    % -- Extract variables for convenience
    Format = OutInfo.Format;
    OutFileName = OutInfo.OutFileName;
    ProblemName = OutInfo.Name(Problem);
 
    LayerZ = Geometry.LayerZ;
    Orient = Geometry.Orient;
    tl = Geometry.tply;
    plymat = Geometry.plymat;
    N = Geometry.N;

    Aij = LamResults.A;
    Bij = LamResults.B;
    Dij = LamResults.D;
    EX = LamResults.EX;
    EY = LamResults.EY;
    NuXY = LamResults.NuXY;
    GXY = LamResults.GXY;
    AlphX = LamResults.AlphX;
    AlphY = LamResults.AlphY;
    AlphXY = LamResults.AlphXY;
    NM = LamResults.NM; 
    NMT = LamResults.NMT;
    EK = LamResults.EK;
    stress = LamResults.stress;
    strain = LamResults.strain;
    MCstress = LamResults.MCstress;
    MCstrain = LamResults.MCstrain;
    
    % -- Formats
    Fmt = '%s \n';
    Fmt1 = '\n %s';

    % -- Open word or text output files
    CurDir=pwd;
    if Format == 'doc'
        file = fullfile(CurDir, OutFileName);
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
            WordDoc = invoke(ActX.Documents, 'Open' ,file);
        else
            WordDoc = invoke(ActX.Documents, 'Add');
        end
        
    else
        fid = fopen(OutFileName,'wt');
        ActX = 0;
        
    end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% -- Write laminate output using WriteText.m and WriteTable.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    style='Title';
    text='Composite Laminate Report';
    WriteText(Format, fid, Fmt, ActX, text, style, [0,1])
    
    text=char(ProblemName);
    WriteText(Format, fid, Fmt, ActX, text, style, [0,1])
    
    style='Normal';
    text=string(datetime);
    WriteText(Format, fid, Fmt, ActX, text, style, [0,1])
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    style='Heading 1';
    text='Applied Loading';
    WriteText(Format, fid, Fmt1, ActX, text, style, [0,0])
        comps(Loads.NM, 1:6)  = ["Nxx", "Nyy", "Nxy", "Mxx", "Myy", "Mxy"];
    comps(Loads.EK, 1:6)  = ["Exx", "Eyy", "Exy", "Kxx", "Kyy", "Kxy"];
    for I = 1:6
        Type(I) = comps(Loads.Type(I), I);
    end

    if isfield(OutInfo, 'INC')
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    % -- Used only for Chapter 7 (Prog Dam)
    if isfield(Loads, 'NINC')

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    style='Heading 1';
    text='Specified Laminate Input Data';
    WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])
    
    DataCell = cell(N+1,5);
    DataCell{1,1} = 'Ply';
    DataCell{1,2} = 'Angle';
    DataCell{1,3} = 'Thickness';
    DataCell{1,4} = 'Material No. & Name';
    DataCell{1,5} = 'Material Type';
    MicroMats = zeros(1,length(plyprops));
    for k = 1:N
        mat = plymat(k);
        DataCell{k+1, 1} = num2str(k);
        DataCell{k+1, 2} = num2str(Orient(k));
        DataCell{k+1, 3} = num2str(tl(k));
        DataCell{k+1, 4} = strcat(num2str(mat),' - ',plyprops{mat}.name);
        DataCell{k+1, 4} = [char(num2str(mat)),' - ',char(plyprops{mat}.name)];

        if isfield(plyprops{mat},'micro')
            DataCell{k+1, 5} = char(plyprops{mat}.micro);
            MicroMats(mat) = 2;
        else
            DataCell{k+1, 5} = 'Eff Props';
            MicroMats(mat) = 1;
        end 

    end 
    
    [NoRows,NoCols]=size(DataCell);      
    enter = [1 0];
    headerRow = true;
    WriteTable(Format, fid, headerRow, ActX, NoRows, NoCols, DataCell, enter);
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

    style='Heading 1';
    text='Ply Materials';
    WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])

    for mat = 1:length(plyprops)
        if MicroMats(mat) ~= 0
            style='Heading 2';
            text=['Material Number',' ',num2str(mat),' - ',char(plyprops{mat}.name)];
            WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])
        end

        if MicroMats(mat) == 1 % Eff props
            style='Heading 3';
            text='Specified Effective Properties';
            WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])

            DataCell = {'E1','E2','Nu12','G12','Alpha_1','Alpha_2'; ...
                        num2str(plyprops{mat}.E1),num2str(plyprops{mat}.E2),num2str(plyprops{mat}.v12), ...
                        num2str(plyprops{mat}.G12),num2str(plyprops{mat}.a1),num2str(plyprops{mat}.a2)};
            [NoRows,NoCols]=size(DataCell);          
            enter = [1 0];
            headerRow = true;
            WriteTable(Format, fid, headerRow, ActX, NoRows, NoCols, DataCell, enter);

          
        elseif MicroMats(mat) == 2 % micromechanics
            plyprops{mat}.Mat = mat;
            WriteEffProps(Format, fid, ActX, plyprops{mat});
            
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

    style='Heading 1';
    text='ABD Matrices';
    WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])

    style='Heading 2';
    text='A Matrix';
    WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])
    
    DataCell = cell(3,3);
    for i = 1:3
        for j = 1:3
            DataCell{i,j} = num2str(Aij(i,j));
        end
    end 

    [NoRows,NoCols]=size(DataCell);          
    enter = [1 0];
    headerRow = false;
    WriteTable(Format, fid, headerRow, ActX, NoRows, NoCols, DataCell, enter);

    style='Heading 2';
    text='B Matrix';
    WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])
    
    DataCell = cell(3,3);
    for i = 1:3
        for j = 1:3
            DataCell{i,j} = num2str(Bij(i,j));
        end
    end 

    [NoRows,NoCols]=size(DataCell);          
    enter = [1 0];
    headerRow = false;
    WriteTable(Format, fid, headerRow, ActX, NoRows, NoCols, DataCell, enter);

    style='Heading 2';
    text='D Matrix';
    WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])
    
    DataCell = cell(3,3);
    for i = 1:3
        for j = 1:3
            DataCell{i,j} = num2str(Dij(i,j));
        end
    end 

    [NoRows,NoCols]=size(DataCell);          
    enter = [1 0];
    headerRow = false;
    WriteTable(Format, fid, headerRow, ActX, NoRows, NoCols, DataCell, enter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

    style='Heading 1';
    text='Inverse ABD Matrices (abd)';
    WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])
    
    ABD = [Aij,Bij;
           Bij,Dij];
    ABDinv = inv(ABD);
    aij=ABDinv(1:3,1:3);
    bij=ABDinv(1:3,4:6);
    dij=ABDinv(4:6,4:6);

    style='Heading 2';
    text='a Matrix';
    WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])

    DataCell = cell(3,3);
    for i = 1:3
        for j = 1:3
            DataCell{i,j} = num2str(aij(i,j));
        end
    end 

    [NoRows,NoCols]=size(DataCell);          
    enter = [1 0];
    headerRow = false;
    WriteTable(Format, fid, headerRow, ActX, NoRows, NoCols, DataCell, enter);

    style='Heading 2';
    text='b Matrix';
    WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])
    
    DataCell = cell(3,3);
    for i = 1:3
        for j = 1:3
            DataCell{i,j} = num2str(bij(i,j));
        end
    end 

    [NoRows,NoCols]=size(DataCell);          
    enter = [1 0];
    headerRow = false;
    WriteTable(Format, fid, headerRow, ActX, NoRows, NoCols, DataCell, enter);


    style='Heading 2';
    text='d Matrix';
    WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])
    
    DataCell = cell(3,3);
    for i = 1:3
        for j = 1:3
            DataCell{i,j} = num2str(dij(i,j));
        end
    end 

    [NoRows,NoCols]=size(DataCell);          
    enter = [1 0];
    headerRow = false;
    WriteTable(Format, fid, headerRow, ActX, NoRows, NoCols, DataCell, enter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    style='Heading 1';
    text='Laminate Effective Properties';
    WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])

    DataCell = {'Ex','Ey','Nuxy','Gxy','Alpha_x','Alpha_y','Alpha_xy'; ...
                num2str(EX),num2str(EY),num2str(NuXY),num2str(GXY),num2str(AlphX),num2str(AlphY),num2str(AlphXY)};
    
    [NoRows,NoCols]=size(DataCell);          
    enter = [1 0];
    headerRow = true;
    WriteTable(Format, fid, headerRow, ActX, NoRows, NoCols, DataCell, enter);

    % -- Don't print loads and through thickness fields for Chapter 7 (prog dam)
    if ~isfield(Loads, 'NINC')
        
        style='Heading 1';
        text='Laminate Force and Moment Resultants';
        WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])

        DataCell = {'Nxx','Nyy','Nxy','Mxx','Myy','Mxy'; ...
                    num2str(NM(1)),num2str(NM(2)),num2str(NM(3)), ...
                    num2str(NM(4)),num2str(NM(5)),num2str(NM(6))};

        [NoRows,NoCols]=size(DataCell);          
        enter = [1 0];
        headerRow = true;
        WriteTable(Format, fid, headerRow, ActX, NoRows, NoCols, DataCell, enter);



        style='Heading 1';
        text='Laminate Midplane Strains and Curvatures';
        WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])

        DataCell = {'Epsx0','Epsy0','Gamxy0','Kappax0','Kappay0','Kappaxy0'; ...
                    num2str(EK(1)),num2str(EK(2)),num2str(EK(3)),num2str(EK(4)),num2str(EK(5)),num2str(EK(6))};

        [NoRows,NoCols]=size(DataCell);          
        enter = [1 0];
        headerRow = true;
        WriteTable(Format, fid, headerRow, ActX, NoRows, NoCols, DataCell, enter);
        
        if any(NMT) ~= 0
            style='Heading 1';
            text='Laminate Thermal Force and Moment Resultants';
            WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])

            DataCell = {'NTxx','NTyy','NTxy','MTxx','MTyy','MTxy'; ...
                        num2str(NMT(1)),num2str(NMT(2)),num2str(NMT(3)), ...
                        num2str(NMT(4)),num2str(NMT(5)),num2str(NMT(6))};

            [NoRows,NoCols]=size(DataCell);          
            enter = [1 0];
            headerRow = true;
            WriteTable(Format, fid, headerRow, ActX, NoRows, NoCols, DataCell, enter);
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

        style='Heading 1';
        text='Ply stresses through the thickness';
        WriteText(Format, fid, Fmt1, ActX, text, style, [1,0])

        DataCell = cell(2*N,9);
        DataCell{1,1} = 'Z-location';
        DataCell{1,2} = 'Angle';
        DataCell{1,3} = 'SigmaX';
        DataCell{1,4} = 'SigmaY';
        DataCell{1,5} = 'TauXY';
        DataCell{1,6} = ' ';
        DataCell{1,7} = 'Sigma11';
        DataCell{1,8} = 'Sigma22';
        DataCell{1,9} = 'Tau12';    
        for kk = 1:2*N
            DataCell{kk+1, 1} = num2str(LayerZ(kk));
            DataCell{kk+1, 2} = num2str(Orient(round(kk/2)));
            DataCell{kk+1, 3} = num2str(stress(1,kk));
            DataCell{kk+1, 4} = num2str(stress(2,kk));
            DataCell{kk+1, 5} = num2str(stress(3,kk));
            DataCell{kk+1, 6} = ' ';
            DataCell{kk+1, 7} = num2str(MCstress(1,kk));
            DataCell{kk+1, 8} = num2str(MCstress(2,kk));
            DataCell{kk+1, 9} = num2str(MCstress(3,kk));
        end 

        [NoRows,NoCols]=size(DataCell);          
        enter = [1 0];
        headerRow = true;
        WriteTable(Format, fid, headerRow, ActX, NoRows, NoCols, DataCell, enter);

    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    % -- Close word or text file
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
    
    close all;
return



    
