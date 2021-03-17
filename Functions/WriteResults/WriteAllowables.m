function WriteAllowables(Format, fid, Fmt1, word, MicroMats, props)
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
% Purpose: Write material allowables to output file
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

WriteHeader = true;

for mat = 1:length(props)
    if MicroMats(mat) == 2

        % -- Determine which constituents have allowables specified
        if isfield(props{mat}.Fiber, 'allowables')
            WriteF = true;
        else
            WriteF = false;
        end

        if isfield(props{mat}.Matrix, 'allowables')
            WriteM = true;
        else
            WriteM = false;
        end

        if isfield(props{mat}, 'Interface')
            if isfield(props{mat}.Interface, 'allowables')  && props{mat}.Vi > 0
                WriteI = true;
            else
                WriteI = false;
            end
        else
            WriteI = false;
        end

        % -- Exit function if no consituents have allowables
        if ~WriteF && ~WriteM && ~WriteI
            return
        end

        if WriteHeader
            style='Heading 2';
            text='Specified Constituent Allowables';
            WriteText(Format, fid, Fmt1, word, text, style, [1,0])
            WriteHeader = false;
        end

        style='Heading 3';
        text=['Material Number',' ',num2str(mat),' - ',char(props{mat}.name)];
        WriteText(Format, fid, Fmt1, word, text, style, [1,0])

        style='Heading 4';
        text='Constituent Stress Allowables';
        WriteText(Format, fid, Fmt1, word, text, style, [1,0])

        DataCell = {'Mat','Name','XT','XC','YT','YC','ZT','ZC','Q','R','S'};

        if WriteF
            DataCellF = {'F', char(props{mat}.Fiber.name),   num2str(props{mat}.Fiber.allowables.XT), ...
                        num2str(props{mat}.Fiber.allowables.XC), num2str(props{mat}.Fiber.allowables.YT), ...
                        num2str(props{mat}.Fiber.allowables.YC), num2str(props{mat}.Fiber.allowables.ZT), ...
                        num2str(props{mat}.Fiber.allowables.ZC), num2str(props{mat}.Fiber.allowables.Q), ...
                        num2str(props{mat}.Fiber.allowables.R),  num2str(props{mat}.Fiber.allowables.S)};
            DataCell = [DataCell; DataCellF];
        end

        if WriteM
            DataCellM = {'M', char(props{mat}.Matrix.name),   num2str(props{mat}.Matrix.allowables.XT), ...
                        num2str(props{mat}.Matrix.allowables.XC), num2str(props{mat}.Matrix.allowables.YT), ...
                        num2str(props{mat}.Matrix.allowables.YC), num2str(props{mat}.Matrix.allowables.ZT), ...
                        num2str(props{mat}.Matrix.allowables.ZC), num2str(props{mat}.Matrix.allowables.Q), ...
                        num2str(props{mat}.Matrix.allowables.R),  num2str(props{mat}.Matrix.allowables.S)};
            DataCell = [DataCell; DataCellM];
        end

        if WriteI
            DataCellI = {'I', char(props{mat}.Interface.name),   num2str(props{mat}.Interface.allowables.XT), ...
                        num2str(props{mat}.Interface.allowables.XC), num2str(props{mat}.Interface.allowables.YT), ...
                        num2str(props{mat}.Interface.allowables.YC), num2str(props{mat}.Interface.allowables.ZT), ...
                        num2str(props{mat}.Interface.allowables.ZC), num2str(props{mat}.Interface.allowables.Q), ...
                        num2str(props{mat}.Interface.allowables.R),  num2str(props{mat}.Interface.allowables.S)};
            DataCell = [DataCell; DataCellI];
        end

        [NoRows,NoCols]=size(DataCell);          
        enter = [1 0];
        headerRow = true;
        WriteTable(Format, fid, headerRow, word, NoRows, NoCols, DataCell, enter);

        style='Heading 4';
        text='Constituent Strain Allowables';
        WriteText(Format, fid, Fmt1, word, text, style, [1,0])

        DataCell = {'Mat','Name','XeT','XeC','YeT','YeC','ZeT','ZeC','Qe','Re','Se'};

        if WriteF
            DataCellF = {'F', char(props{mat}.Fiber.name),   num2str(props{mat}.Fiber.allowables.XeT), ...
                        num2str(props{mat}.Fiber.allowables.XeC), num2str(props{mat}.Fiber.allowables.YeT), ...
                        num2str(props{mat}.Fiber.allowables.YeC), num2str(props{mat}.Fiber.allowables.ZeT), ...
                        num2str(props{mat}.Fiber.allowables.ZeC), num2str(props{mat}.Fiber.allowables.Qe), ...
                        num2str(props{mat}.Fiber.allowables.Re),  num2str(props{mat}.Fiber.allowables.Se)};
            DataCell = [DataCell; DataCellF];
        end

        if WriteM
            DataCellM = {'M', char(props{mat}.Matrix.name),   num2str(props{mat}.Matrix.allowables.XeT), ...
                        num2str(props{mat}.Matrix.allowables.XeC), num2str(props{mat}.Matrix.allowables.YeT), ...
                        num2str(props{mat}.Matrix.allowables.YeC), num2str(props{mat}.Matrix.allowables.ZeT), ...
                        num2str(props{mat}.Matrix.allowables.ZeC), num2str(props{mat}.Matrix.allowables.Qe), ...
                        num2str(props{mat}.Matrix.allowables.Re),  num2str(props{mat}.Matrix.allowables.Se)};
            DataCell = [DataCell; DataCellM];
        end

        if WriteI
            DataCellI = {'I', char(props{mat}.Interface.name),   num2str(props{mat}.Interface.allowables.XeT), ...
                        num2str(props{mat}.Interface.allowables.XeC), num2str(props{mat}.Interface.allowables.YeT), ...
                        num2str(props{mat}.Interface.allowables.YeC), num2str(props{mat}.Interface.allowables.ZeT), ...
                        num2str(props{mat}.Interface.allowables.ZeC), num2str(props{mat}.Interface.allowables.Qe), ...
                        num2str(props{mat}.Interface.allowables.Re),  num2str(props{mat}.Interface.allowables.Se)};
            DataCell = [DataCell; DataCellI];
        end

        [NoRows,NoCols]=size(DataCell);          
        enter = [1 0];
        headerRow = true;
        WriteTable(Format, fid, headerRow, word, NoRows, NoCols, DataCell, enter);

    end    
end

end