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
% Purpose: Driver script for laminate progressive damage analysis using CLT. The 
%          damage occurs at the micro scale (not the ply scale). The problem input 
%          is defined in the function LamProblemDef.m and GetEffProps.m
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% -- Clear memory and close files
clear;
close all;
fclose('all');
clc;

% -- Add needed function locations to the path
addpath('Functions/CLT');
addpath('Functions/Utilities');
addpath('Functions/WriteResults');
addpath('Functions/Micromechanics');
addpath('Functions/Margins');
addpath('Functions/Damage');

%-----------------------------------------------------------------
% 1) Define Laminate Problems
%-----------------------------------------------------------------
[NProblems, OutInfo, Geometry, Loads] = LamProblemDef();

% -- Perform some input checks
for NP = 1: NProblems
    
    if ~isfield(Loads{NP}, 'DT')
        Loads{NP}.DT = 0;
    elseif Loads{NP}.DT ~= 0
       error(['Thermal problem not implemented for progressive damage, DT ~= 0 in problem #', char(num2str(NP))]); 
    end
    
    if OutInfo.Format == "doc" 
        disp(' ---> Word output not supported for progressive damage, switching to txt');
        OutInfo.Format = "txt";
    end
    
end

% -- Preallocate LamResults
LamResults = cell(1,NProblems);

%-----------------------------------------------------------------
% 2) Get initial ply properties
%-----------------------------------------------------------------
plyprops = GetPlyProps(Geometry);

% -- Loop through problems
for NP = 1: NProblems
    
    % -- Clear storage of load/strain history for next problem
    if NP > 1
        clear NMhist;
        clear EKhist;
        clear fullhist;
    end
    
    OutInfo.PlottedInit = false;
    
    % -- Check for missing problem name
    if (ismissing(OutInfo.Name(NP)))
        OutInfo.Name(NP) = string(['Problem ', char(num2str(NP))]);
    end
    
    % -- Check for missing ply angles, thicknesses, materials
    if ~isfield(Geometry{NP},'Orient') || ~isfield(Geometry{NP},'tply') || ...
       ~isfield(Geometry{NP},'plymat')
        error(['Orient, tply, or plymat missing, Problem #', num2str(NP)]);
    end
        
    % -- Echo problem info to command window
    disp(['Problem #',num2str(NP),' - ', char(OutInfo.Name(NP))]);
    disp(['   Per ply angle orientations   [',num2str(Geometry{NP}.Orient), ']']);
    disp(['   Per ply thicknesses          [',num2str(Geometry{NP}.tply), ']']);
    disp(['   Per ply material assignments [',num2str(Geometry{NP}.plymat), ']']);    

    % -- Check that the problem's ply materials, orientation, and thickness
    %    have consistent lengths
    if length(Geometry{NP}.tply) ~= length(Geometry{NP}.Orient) || ...
       length(Geometry{NP}.tply) ~= length(Geometry{NP}.plymat) || ...
       length(Geometry{NP}.Orient) ~= length(Geometry{NP}.plymat)
            error('Lengths of tply, Orient, and plymat are not consistent');
    end
    
    % -- Check that the problem's ply materials have been defined and are micro-based
    for k = 1: length(Geometry{NP}.tply)
        if ~isfield(plyprops{Geometry{NP}.plymat(k)}, 'name')
            error(strcat('ply material #', num2str(Geometry{NP}.plymat(k)), ' undefined ... check GetPlyProps'));
        end
        
        if ~isfield(plyprops{Geometry{NP}.plymat(k)}, 'micro')
            error(strcat('ply material #', num2str(Geometry{NP}.plymat(k)), ' is not micromechanics-based', ...
                         ' and cannot be used in a progressive damage simulation'));
        end
        
    end
    
    % -- Check for slash character in OutInfo.Name
    k = strfind(OutInfo.Name(NP), '/');
    j = strfind(OutInfo.Name(NP), '\');
    if ~isempty(k) || ~isempty(j)
        error('Problem name contains a slash or backslash ... remove');
    end
    
    % -- Specifed loads are final values for incremental loading
    Loads{NP}.Final = Loads{NP}.Value;
    Loads{NP}.DTFinal = Loads{NP}.DT;
        
%-------------------------------------------------------------------
% 3) Make unique copy of each constituent per ply (so ply can damage
%      independently) and create damaged constituents
%-------------------------------------------------------------------
    DamageFactor = 0.0001; % -- Stiffness reduction factor for damaged material
    for k = 1:length(Geometry{NP}.tply)
        mat = Geometry{NP}.plymat(k);
        matdam = k + 1000;
        plyprops{matdam} = plyprops{mat};
        Geometry{NP}.plymat(k) = matdam;
        % -- Save Damaged Constituents
        [plyprops{matdam}.FiberDam] = DamageConstit(plyprops{matdam}.Fiber, DamageFactor);
        [plyprops{matdam}.MatrixDam] = DamageConstit(plyprops{matdam}.Matrix, DamageFactor);
        if isfield(plyprops{matdam},'Interface')
            [plyprops{matdam}.InterfaceDam] = DamageConstit(plyprops{matdam}.Interface, DamageFactor);
        end
        
    end
    
    % -- Default to 100 increments if NINC not specified
    if ~isfield(Loads{NP}, 'NINC')
        Loads{NP}.NINC = 100;
    end
    
    % -------------------------
    % -- Loading increment loop
    % -------------------------
    INC = -1;
    INCcount = 0;
    while INC < Loads{NP}.NINC

        INC = INC + 1;
        INCcount = INCcount + 1;
        
        OutInfo.INC = INC;
        
        % -- Determine current load level
        pct = INC/Loads{NP}.NINC;
        Loads{NP}.Value = pct*Loads{NP}.Final;
        Loads{NP}.DT = pct*Loads{NP}.DTFinal;

        % -- Display increment progress 
        disp(' ');
        disp('*****************************************************');
        disp(['* INC = ',num2str(INC), '  out of ', num2str(Loads{NP}.NINC)]);
        disp('*****************************************************');
        
        % -----------------
        % -- Iteration Loop
        % -----------------
        NITER = 1000000000; % -- Max number of iterations
        for ITER = 1:NITER
            
            % -- Display iteration progress 
            disp(' ');
            disp(['----- ITER = ', char(num2str(ITER)), ' -----']);
            OutInfo.ITER = ITER;
        
            %-----------------------------------------------------------------
            % 4) Analyze laminate with CLT
            %-----------------------------------------------------------------
            [Geometry{NP}, LamResults{NP}] = CLT(plyprops, Geometry{NP}, Loads{NP}, LamResults{NP});
            
            % -- Check to make sure no bending (see Section 7.1.1)
            TinyBending = 1.E-12;
            if ~all(abs(LamResults{NP}.EK(4:6)) < TinyBending)
                error(strcat('Problem #', num2str(NP), ' - Laminate is experiencing bending. This capability NOT implemented for progressive damage.'));
            end
            
            %-------------------------------------------------------------------------
            % 5) Write initial laminate properties and plot/write incremental results
            %-------------------------------------------------------------------------
            if INC == 0
                [OutInfo] = OutputLam(OutInfo, NP, Geometry{NP}, Loads{NP}, plyprops, LamResults{NP});
            end

            % -- Write global output
            if INC == 0
                OutInfo.DamDir = ['Output/DamSnapShots - ',char(OutInfo.Name(NP)),' ',OutInfo.datetime,'/'];
                [~,~,~] = mkdir(OutInfo.DamDir);

                filename = [OutInfo.DamDir, 'ITER Load vs eps0 - ', char(OutInfo.Name(NP)), ' ', OutInfo.datetime,'.txt']; 
                fid = fopen(filename, 'wt');

                filename3 = [OutInfo.DamDir, 'ITER props vs INC - ', char(OutInfo.Name(NP)),'.txt']; 
                fid3 = fopen(filename3, 'a');
                
                Fmt1 = '\r %s';
                text = 'INC ITER eps0_xx Nxx eps0_yy Nyy eps0_xy Nxy kappa_xx Mxx kappa_yy Myy kappa_xy Mxy';
                style='Normal';
                WriteText('txt', fid, Fmt1, 0, text, style, [0,0])

                text = 'INC ITER EX EY NuXY GXY AlphaX AlphaY AlphaXY';
                style='Normal';
                WriteText('txt', fid3, Fmt1, 0, text, style, [0,0])
                
                % -- Inverse of ABD
                AInv = inv(LamResults{NP}.A);

                % -- Store initial trace(inv Cstar) - used as global compliance measure
                InitTraceA = AInv(1,1) + AInv(2,2);
                TraceA = InitTraceA;

            end

            % -- Loads{NP}.PlotComponent can be specified in MicroProblemDef.m, or 
            %    will plot applied direction if only 1 component applied
            if isfield(Loads{NP}, 'PlotComponent')
                PlotComponent = Loads{NP}.PlotComponent;
            else
                PlotComponent = 2;
                for i = 1: 6
                    if Loads{NP}.Final(i) ~= 0
                        PlotComponent = i;
                        break;
                    end
                end
            end
            
            % -- Store global load and strain history
            EKhist(INCcount) = LamResults{NP}.EK(PlotComponent);
            NMhist(INCcount) = LamResults{NP}.NM(PlotComponent);
            fullhist(INCcount).NM = LamResults{NP}.NM;
            fullhist(INCcount).EK = LamResults{NP}.EK;
            
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % -- Write and plot progressive damage results as damage progresses
            %    (Note: this is done as a function in the case of MicroSimDamage.m)
            Fmt1 = '\r %s';
            text = [num2str(INC), ' ',num2str(ITER), ' ', ...
                    num2str(LamResults{NP}.EK(1)),' ',num2str(LamResults{NP}.NM(1)),' ',num2str(LamResults{NP}.EK(2)),' ', ...
                    num2str(LamResults{NP}.NM(2)),' ',num2str(LamResults{NP}.EK(3)),' ',num2str(LamResults{NP}.NM(3)),' ', ...
                    num2str(LamResults{NP}.EK(4)),' ',num2str(LamResults{NP}.NM(4)),' ',num2str(LamResults{NP}.EK(5)),' ', ...
                    num2str(LamResults{NP}.NM(5)),' ',num2str(LamResults{NP}.EK(6)),' ',num2str(LamResults{NP}.NM(6))];
            style='Normal';
            WriteText('txt', fid, Fmt1, 0, text, style, [0,0])

            text = [num2str(INC),' ',num2str(ITER),' ',num2str(LamResults{NP}.EX),' ',num2str(LamResults{NP}.EY),' ',num2str(LamResults{NP}.NuXY),' ',num2str(LamResults{NP}.GXY),' ', ...
                    num2str(LamResults{NP}.AlphX),' ',num2str(LamResults{NP}.AlphY),' ',num2str(LamResults{NP}.AlphXY)];
            style='Normal';
            WriteText('txt', fid3, Fmt1, 0, text, style, [0,0])
            
            Fmt1 = '\r %s';
            style='Normal';
            filename2 = [OutInfo.DamDir, 'Load vs eps0 - ', char(OutInfo.Name(NP)), ' ', OutInfo.datetime, '.txt']; 
            fid2 = fopen(filename2, 'wt');
            Fmt1 = '\r %s';
            text = 'eps0_xx Nxx eps0_yy Nyy eps0_xy Nxy kappa_xx Mxx kappa_yy Myy kappa_xy Mxy';
            style='Normal';
            WriteText('txt', fid2, Fmt1, 0, text, style, [0,0])
            for I = 1: length(fullhist)
                text = [num2str(fullhist(I).EK(1)),' ',num2str(fullhist(I).NM(1)),' ',num2str(fullhist(I).EK(2)),' ',num2str(fullhist(I).NM(2)),' ', ...
                        num2str(fullhist(I).EK(3)),' ',num2str(fullhist(I).NM(3)),' ',num2str(fullhist(I).EK(4)),' ',num2str(fullhist(I).NM(4)),' ',...
                        num2str(fullhist(I).EK(5)),' ',num2str(fullhist(I).NM(5)),' ',num2str(fullhist(I).EK(6)),' ',num2str(fullhist(I).NM(6))];
                WriteText('txt', fid2, Fmt1, 0, text, style, [0,0])
            end
            fclose(fid2);

            if INC == 0
                damfig = 0;
                LamResults{NP}.NewFailure = false;
            end

            if ~isfield(Loads{NP},'PlyDamWatchOn')
                Loads{NP}.PlyDamWatchOn(1:Geometry{NP}.N) = 1;
            end

            subplotcount = sum(Loads{NP}.PlyDamWatchOn(1:Geometry{NP}.N)) + 1;

            % -- Plot NM vs. EK
            if INC == 0
                damfig = figure('units','normalized','outerposition',[0 0.5 1 0.5]);
            else
                figure(damfig);
                subplot(1,subplotcount,1)
                plot(EKhist,NMhist,'-o')
                title('NM vs. EK');
                xlabel('EK'); 
                ylabel('NM (N/mm or N-mm/mm)','rotation',90);
                drawnow
            end

            % -- Display and save stress-strain plot, local von Mises stress and damage
            if LamResults{NP}.NewFailure || INC == 1 && (plyprops{mat}.micro == "GMC" || plyprops{mat}.micro == "HFGMC") 

                % -- Only plot first point in each ply
                for k = 1: Geometry{NP}.N

                    if Loads{NP}.PlyDamWatchOn(k) == 0
                        continue
                    end

                    mat = Geometry{NP}.plymat(k);

                    if (plyprops{mat}.micro == "GMC" || plyprops{mat}.micro == "HFGMC")

                        mats = 0;

                        % -- Plot failed subcells
                        for j = 1:2*plyprops{mat}.RUC.NG
                            g = round(j/2);
                            for i = 1:2*plyprops{mat}.RUC.NB
                                b = round(i/2);
                                switch plyprops{mat}.RUC.matsCh(b,g)
                                    case 'F'
                                        mats(i,j) = 1;
                                    case 'M'
                                        mats(i,j) = 2;
                                    case 'I'
                                        mats(i,j) = 1.6;
                                    otherwise
                                        mats(i,j) = 3;
                                end
                            end 
                        end 

                        x2 = 0;
                        x3 = 0;
                        % -- place an x2-point at bottom and top of each subcell
                        h = 0;
                        i = 1;
                        for b = 1:plyprops{mat}.RUC.NB
                            x2(i) = h;
                            h = h + plyprops{mat}.RUC.h(b);
                            i = i + 1;
                            x2(i) = h;
                            i = i + 1;
                        end

                        % -- place an x3-point at left and right of each subcell
                        l = 0;
                        i = 1;
                        for g = 1:plyprops{mat}.RUC.NG
                            x3(i) = l;
                            l = l + plyprops{mat}.RUC.l(g);
                            i = i + 1;
                            x3(i) = l;
                            i = i + 1;
                        end

                       % -- Create grid
                        X = 0;
                        Y = 0;
                        [X,Y] = meshgrid(x3,x2);        

                        % --- Plot damage ---
                        subplot(1,subplotcount,k+1)
                        colormap(jet);
                        %pcolor( X, Y, Failure), shading interp;
                        pcolor( X, Y, mats), shading faceted;
                        caxis([1 3]);
                        fontsize = 14;
                        c.FontSize=fontsize;
                        c.FontWeight='bold';
                        title(['ANG ',char(num2str(Geometry{NP}.Orient(k))),', PLY ', char(num2str(k)), ', INC ',char(num2str(INC))],'FontSize',fontsize);
                        set (gca,'FontSize',fontsize,'FontWeight','bold','LineWidth',2)
                        axis image;
                        axis off; 
                        drawnow

                        % -- Save damage snapshot figure
                        if OutInfo.Format == "doc"
                            WriteWordFig([pwd,'/', OutInfo.OutFile, '.doc'])
                        elseif OutInfo.Format == "txt"
                            if ~OutInfo.PlottedInit && k == subplotcount - 1
                                if ITER == 2
                                    Plotname = [OutInfo.DamDir,'Initiation Location Damage Snapshot  - ', 'INC = ', char(num2str(INC)), ' ITER = ', char(num2str(ITER)),'.bmp'];
                                    saveas(gcf,Plotname);
                                    OutInfo.PlottedInit = true;
                                end
                            end
                            OutInfo.DamDir = ['Output/DamSnapShots - ',char(OutInfo.Name(NP)),' ',OutInfo.datetime,'/'];
                            [~,~,~] = mkdir(OutInfo.DamDir);
                            Plotname = [OutInfo.DamDir,'Damage Snapshot  - ', 'INC = ', char(num2str(INC)),'.bmp'];
                            saveas(gcf,Plotname);
                        end   

                    end
                end

            end
            % -- End write and plot progressive damage results section
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         
            % -- Stop execution for if TraceA exceeds termination factor
            if isfield(Loads{NP}, 'TerminationFactor')
                factor = Loads{NP}.TerminationFactor;
            else
                factor = 50000; % -- Defaults to very high factor
            end

            if TraceA/InitTraceA > factor
                INC = Loads{NP}.NINC + 1;
                break
            end

            %-----------------------------------------------------------------
            % 6) Calculate local fields for micromechanics
            %-----------------------------------------------------------------
            [LamResults{NP}] = LamMicroFields(Geometry{NP}, Loads{NP}, plyprops, LamResults{NP});

            %-----------------------------------------------------------------
            % 7) Calculate Margins
            %-----------------------------------------------------------------
            % -- Check for turning off some failure criteria (all on by default, 0 = off)
            if (~isfield(Loads{NP},'CriteriaOn'))
                Loads{NP}.CriteriaOn = [1,1,1,1];
            end
            [LamResults{NP}, plyprops] = CalcLamMargins(plyprops, Geometry{NP}, Loads{NP}, LamResults{NP}, true); 

            % -- Initialize failure flags
            if INC == 0
                 for kk = 1: 2*Geometry{NP}.N
                     k = round(kk/2);
                     mat = Geometry{NP}.plymat(k);
                     if LamResults{NP}.Micro(kk).Type == "FM"
                        LamResults{NP}.Micro(kk).MoS.MicroF.Failed = false;
                        LamResults{NP}.Micro(kk).MoS.MicroM.Failed = false;
                     else
                        for b = 1: plyprops{mat}.RUC.NB
                            for g = 1: plyprops{mat}.RUC.NG
                                LamResults{NP}.Micro(kk).MoS.Micro(b,g).Failed = false;
                            end 
                        end
                     end
                end
            end
            
            %-----------------------------------------------------------------
            % 8) Reduce stiffness of any subcell with negative min margin
            %-----------------------------------------------------------------
            [plyprops, LamResults{NP}] = DamageLam(plyprops, Geometry{NP}, LamResults{NP}, OutInfo);

            % -- Check for new failure
            LamResults{NP}.NewFailure = false;
            if INC > 0
                for kk = 1: 2*Geometry{NP}.N
                    k = round(kk/2);
                    mat = Geometry{NP}.plymat(k);
                    if (plyprops{mat}.micro == "GMC" || plyprops{mat}.micro == "HFGMC")
                        if LamResults{NP}.Micro(kk).MoS.NewFailure
                            LamResults{NP}.NewFailure = true;
                            break;
                        end
                    end
                end
            end
            
            % -- Find min MoS to skip (jump) increments
            CheckForSkip = true;
            if ITER ~= 1
                CheckForSkip = false;
            elseif isfield(LamResults{NP}, 'NewFailure')
                if LamResults{NP}.NewFailure
                    CheckForSkip = false;
                end
            end
            
            if CheckForSkip
                MinMoS = 99999;
                for kk = 1: 2*Geometry{NP}.N
                    k = round(kk/2);
                    mat = Geometry{NP}.plymat(k);
                    if LamResults{NP}.Micro(kk).Type == "FM"
                        if isfield(LamResults{NP}.Micro(kk).MoS.MicroF, 'MinMoS')
                            MinMoS = min(MinMoS, LamResults{NP}.Micro(kk).MoS.MicroF.MinMoS);
                        end
                        if isfield(LamResults{NP}.Micro(kk).MoS.MicroM, 'MinMoS')
                            MinMoS = min(MinMoS, LamResults{NP}.Micro(kk).MoS.MicroM.MinMoS);
                        end
                    else
                        for b = 1: plyprops{mat}.RUC.NB
                            for g = 1: plyprops{mat}.RUC.NG
                                 
                                if ~isfield(LamResults{NP}.Micro(kk).MoS.Micro(b,g), 'MinMoS')
                                    continue;
                                end
                                 
                                if isfield(LamResults{NP}.Micro(kk).MoS.Micro(b,g), 'Failed') 
                                    if LamResults{NP}.Micro(kk).MoS.Micro(b,g).Failed
                                        continue;
                                    else
                                        MinMoS = min(MinMoS, LamResults{NP}.Micro(kk).MoS.Micro(b,g).MinMoS);
                                    end 
                                else
                                    MinMoS = min(MinMoS, LamResults{NP}.Micro(kk).MoS.Micro(b,g).MinMoS);
                                end
                                
                            end
                        end
                    end
                     
                 end
                   
                disp(['  MinMoS = ',char(num2str(MinMoS))]);

                % -- Step out to just before increment where next damage should occur
                %    if MinMoS is positive (no additional failures)
                if MinMoS > 0 && MinMoS < 99999 && INC < Loads{NP}.NINC
                    INC_JUMP = round( INC*( (MinMoS + 1) - 1) );
                    disp(['  INC_JUMP = ', char(num2str(INC_JUMP))]);
                    if INC_JUMP > 2
                        INC = INC + INC_JUMP - 2;
                        % -- Jump to final increment if will no more damage will occur
                        if INC > Loads{NP}.NINC - 1
                            INC = Loads{NP}.NINC - 1;
                        end
                        break;
                    end
                elseif MinMoS == 99999 && INC > 10 && INC < Loads{NP}.NINC
                    INC = Loads{NP}.NINC - 1;
                    break;
                end
                
             end
                

            % -- Inverse of ABD
            AInv = inv(LamResults{NP}.A);

            % -- Check for large stiffness drop (Ainv increase)
            TraceA = AInv(1,1) + AInv(2,2);
            disp(['  --> Ainv Trace Factor = ',char(num2str(TraceA/InitTraceA))]);

            % -- Exit iteration loop if no new failure
            if ~LamResults{NP}.NewFailure
                break
            end

        end % -- End Iteration Loop
        
    end % -- End Increment Loop

    fclose(fid);
    
    %-----------------------------------------------------------------
    % 9) Plot micro fields and write RUC (only for RUCid = 300 or 1000)
    %-----------------------------------------------------------------
    if (~isfield(OutInfo,'Format'))
        OutInfo.Format = "txt";
    end
    [plyprops] = OutputMicroFields(OutInfo, Geometry{NP}, plyprops, LamResults{NP});
    
    WriteRUC(OutInfo, OutInfo.Name(NP), plyprops, Geometry{NP})

    
    disp(['  *** Problem ',char(num2str(NP)),' Completed ***'])
    disp(' ');
    
    figh = findall(0,'type','figure');
    other_figures = setdiff(figh, damfig);
    delete(other_figures)
    
end