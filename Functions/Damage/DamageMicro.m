function [props, MoS] = DamageMicro(props, MoS, OutInfo, CurrentPly)
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
%
% Purpose: Replace subcells or materials with lowest negative MoS with damaged 
%          materials (i.e., reduced stiffnesses) and run micromechanics analysis to
%          get updated effective properties and concentration tensors
% Input:
% - props: Struct containing material properties 
% - MoS: Struct containing MoS results and info
% - OutInfo: Struct containing output information
% - CurrentPly: Current ply number within laminate (optional, laminate only)
% Output:
% - props: Updated struct containing material properties
% - MoS: Updated struct containing MoS results and info
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -- If not laminate, set ply = 0
if nargin == 3
    CurrentPly = 0;
end

% -- Write damage info to output file (txt based only, not for Word)
if OutInfo.INC == 0 && CurrentPly <= 1 && OutInfo.Format == "txt"
    fid = fopen([OutInfo.OutFile, '.txt'],'a');
    Fmt = '%s \n';
    ActXWord = 0;
    text = ' ';
    WriteText(OutInfo.Format, fid, Fmt, ActXWord, text, 'Normal', [0,1]);
    text = 'Incremental Progressive Damage Output';
    WriteText(OutInfo.Format, fid, Fmt, ActXWord, text, 'Heading 1', [0,1]);
end

if OutInfo.Format == "txt"
    fid = fopen([OutInfo.OutFile, '.txt'],'a');
    Fmt = '%s \n';
    ActXWord = 0;
    if CurrentPly <= 1
        if OutInfo.ITER <= 1
            text = ' ';
            WriteText(OutInfo.Format, fid, Fmt, ActXWord, text, 'Normal', [0,1]);
            text = ['*** Increment Number ', char(num2str(OutInfo.INC)), ' ***'];
            WriteText(OutInfo.Format, fid, Fmt, ActXWord, text, 'Normal', [0,1]);
        end
        text = ['  -- Iteration Number ', char(num2str(OutInfo.ITER)), ' --'];
        WriteText(OutInfo.Format, fid, Fmt, ActXWord, text, 'Normal', [0,1]);
    end
end

Theory = char(props.micro);

% --------------------------
% -- Theories withOUT RUCs
% --------------------------
if (props.micro == "MT" || props.micro == "Voigt" || props.micro == "Reuss")

    % -- Check if F or M material is already failed
    AlreadyFailedM = false;
    AlreadyFailedF = false;
    if isfield(MoS.MicroM, 'Failed')
        if MoS.MicroM.Failed
            AlreadyFailedM = true;
        end
    end
    if isfield(MoS.MicroF, 'Failed')
        if MoS.MicroF.Failed
            AlreadyFailedF = true;
        end
    end

    % -- If FIBER MoS negative or already failed, replace with damaged material
    if MoS.MicroF.MinMoS < 0 || AlreadyFailedF
        MoS.MicroF.Failed = true;
        Constits.Fiber = props.FiberDam;
        constitprops{props.FiberDam.constID} = props.FiberDam;
        if ~AlreadyFailedF
            CritNum = MoS.MicroF.ControllingNum;
            CritName = MoS.MicroF.Controlling;
            Component = MoS.MicroF.Crit{CritNum}.Controlling;
            if CurrentPly > 0
                text = ['Ply ', num2str(CurrentPly), ' - FIBER Failure - ', char(CritName), ' ', char(num2str(Component))];
            else
                text = ['FIBER Failure - ', char(CritName), ' ', char(num2str(Component))];
            end
            disp(text);
            if OutInfo.Format == "txt"
                WriteText(OutInfo.Format, fid, Fmt, ActXWord, text, 'Normal', [0,1])
            end
        end
    else
        Constits.Fiber = props.Fiber;
        constitprops{props.Fiber.constID} = props.Fiber;
    end
    
    % -- If MATRIX MoS negative or already failed, replace with damaged material
    if MoS.MicroM.MinMoS < 0 || AlreadyFailedM
        MoS.MicroM.Failed = true;
        Constits.Matrix = props.MatrixDam;
        constitprops{props.MatrixDam.constID} = props.MatrixDam;
        if ~AlreadyFailedM
            CritNum = MoS.MicroM.ControllingNum;
            CritName = MoS.MicroM.Controlling;
            Component = MoS.MicroM.Crit{CritNum}.Controlling;
            if CurrentPly > 0
                text = ['Ply ', num2str(CurrentPly), ' - MATRIX Failure - ', char(CritName), ' ', char(num2str(Component))];
            else
                text = ['MATRIX Failure - ', char(CritName), ' ', char(num2str(Component))];
            end
            disp(text);
            if OutInfo.Format == "txt"
                WriteText(OutInfo.Format, fid, Fmt, ActXWord, text, 'Normal', [0,1])
            end
        end
    else
        Constits.Matrix = props.Matrix;
        constitprops{props.Matrix.constID} = props.Matrix;
    end

    % -- Run micromechanics to get updated effective properties
    [props] = RunMicro(Theory, props.name, props.Vf, Constits, props);


% --------------------------
% -- Theories WITH RUCs
% --------------------------
else
    % -- Store damaged F, M, and I, if they exist
    Constits.Fiber = props.Fiber;
    if isfield(props,'FiberDam')
        Constits.FiberDam = props.FiberDam;              
    end
    Constits.Matrix = props.Matrix;
    if isfield(props,'MatrixDam')
        Constits.MatrixDam = props.MatrixDam;              
    end
    if isfield(props,'Interface')
        Constits.Interface = props.Interface;
        if isfield(props,'InterfaceDam')
            Constits.InterfaceDam = props.InterfaceDam;              
        end
    end

    % -- Find min MOS of the unfailed subcells in RUC 
    minMoS = 0;
    for b = 1: props.RUC.NB
        for g = 1: props.RUC.NG
            if ~isfield(MoS.Micro(b,g), 'MinMoS')
                continue
            end
            if isfield(MoS.Micro(b,g), 'Failed')
                if ~MoS.Micro(b,g).Failed 
                    minMoS = min(minMoS, MoS.Micro(b,g).MinMoS);
                end
            else
                minMoS = min(minMoS, MoS.Micro(b,g).MinMoS);
            end
        end        
    end

    % -- Tolerance for failing multiple subcells
    % -- Only fail subcells that are this close on MoS to the minimum
    tol = 0.001;
    %tol = -minMoS; % -- Turn off tol (all negative margin subcells fail at once)

    MoS.NewFailure = false; % -- Flag indicating if a new failure has occurred
    for b = 1: props.RUC.NB
        for g = 1: props.RUC.NG

            if ~isfield(MoS.Micro(b,g), 'MinMoS')
                continue
            end

            AlreadyFailed = false;
            if isfield(MoS.Micro(b,g), 'Failed')
                if MoS.Micro(b,g).Failed 
                    AlreadyFailed = true;
                end
            end

            % -- If MoS negative and within tolerance or already failed, 
            %    replace subcell mat with damaged material
             if (minMoS < 0) && (MoS.Micro(b,g).MinMoS - minMoS - tol) < 0 || AlreadyFailed
                if ~AlreadyFailed
                    MoS.Micro(b,g).Failed = true;
                    MoS.NewFailure = true;
                    CritNum = MoS.Micro(b,g).ControllingNum;
                    CritName = MoS.Micro(b,g).Controlling;
                    Component = MoS.Micro(b,g).Crit{CritNum}.Controlling;
                    if CurrentPly > 0
                        text = ['      Ply ', num2str(CurrentPly),' - Failure in subcell ',num2str(b),',',num2str(g), ' - ', char(CritName), ' ', char(num2str(Component))];
                    else
                        text = ['      Failure in subcell ',num2str(b),',',num2str(g), ' - ', char(CritName), ' ', char(num2str(Component))];
                    end
                    disp(text);
                    if OutInfo.Format == "txt"
                        WriteText(OutInfo.Format, fid, Fmt, ActXWord, text, 'Normal', [0,1])
                        if isfield(OutInfo, 'DamDir')
                              status = copyfile(['Output/*',OutInfo.datetime,'*.txt'], OutInfo.DamDir);
                        end
                    end
                end

                % -- Change character representation of damaged material from F,M,I to 1,2,3
                switch(props.RUC.matsCh(b,g))
                    case {'F', '1'}
                        props.RUC.mats(b,g) = props.FiberDam.constID;
                        props.RUC.matsCh(b,g) = '1';
                    case {'M', '2'}
                        props.RUC.mats(b,g) = props.MatrixDam.constID;
                        props.RUC.matsCh(b,g) = '2';
                    case {'I', '3'}
                        props.RUC.mats(b,g) = props.InterfaceDam.constID;
                        props.RUC.matsCh(b,g) = '3';
                    otherwise
                        error(['Invalid RUC.matsch in DamageMicro - RUC.matsch = ', RUC.matsch])
                end

             end
        end
    end

    % -- If there has been a new failure, run micromechanics to get updated effective 
    %    properties and concentration tensors
    if MoS.NewFailure

        if isfield(props, 'Vi')
            Vol.Vi = props.Vi;
            Vol.Vf = props.Vf;
        else
            Vol = props.Vf;
        end

        [props] = RunMicro(Theory, props.name, Vol, Constits, props, props.RUC);

    end

end

% -- Close output file
if OutInfo.Format == "txt"
    fclose(fid);
end

    
end