function PlotMicroFields(OutInfo, props, Results, Angle, ply, kk, bending, Zloc)
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
% Purpose: Make pseudocolor plots of the local fields in the composite based on 
%          micromechanics analyses.  For CLT, this is called per top/bottom point in
%          each ply (or 1 pt per ply if sym and no bending)
% Input:
% - OutInfo: Struct containing output information
% - props: Struct containing composite material and constituent properties 
% - Results: Struct containing micromechanics analysis results
% - Angle: Current ply angle if laminate (optional)
% - ply: Current ply number if lamainte (optional)
% - kk: Current through-thickness point number if laminate(optional)
% - bending: Flag indicating is bending is present if lamiante (optional)
% - Zloc: Through-thickness z-coordinate location in laminate (optional)
% Output: None
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if nargin == 3 % -- Nonlaminate case (will have only 3 arguments)
        Angle = 0;
        ply = 0;
        bending = false;
        Zloc = 0;
        AddtoTitle = "";
        AddPly = '';
    else % -- Laminate case
        if (bending) 
            if mod(kk,2) == 1
                btlabel = "B"; % -- Label for bottom of ply
            else
                btlabel = "T"; % -- Label for top of ply
            end
        else
            btlabel = "";
        end
        AddtoTitle = strcat(" in Ply No. ",num2str(ply),btlabel,", Z = ", ...
                            num2str(Zloc),", ply Angle ",num2str(Angle));
        AddPly = [' - Ply No.',num2str(ply)];
    end
    
    % -- Define an RUC for MOC and MT to plot the local fields on
    if ((props.micro == "MT")   || (props.micro == "MOC")  || ...
        (props.micro == "MOCu") ||(props.micro == "Voigt") || ...
        (props.micro == "Reuss"))
        RUC.NB = 2;
        RUC.NG = 2;
        RUC.Vf = props.Vf;
        RUC.id = 2;
        RUC.h(1) = sqrt(props.Vf);
        RUC.h(2) = 1-RUC.h(1);
        RUC.l = RUC.h; 
        RUC.fiber = props.Fiber.name;
        RUC.matrix = props.Matrix.name;
        RUC.mats = [props.Fiber.constID, props.Matrix.constID; ...
                    props.Matrix.constID, props.Matrix.constID];
        props.RUC = RUC;
        props.RUC.matsCh = ['F','M';'M','M'];
        if (props.micro == "MT" || props.micro == "Voigt" || ...
            props.micro == "Reuss") % -- Copy F and M fields to RUC subcells
          MicroFields(1,1).strain = Results.MicroFieldsF.strain;
          MicroFields(1,1).stress = Results.MicroFieldsF.stress;
          MicroFields(1,2).strain = Results.MicroFieldsM.strain;
          MicroFields(1,2).stress = Results.MicroFieldsM.stress;
          MicroFields(2,1).strain = Results.MicroFieldsM.strain;
          MicroFields(2,1).stress = Results.MicroFieldsM.stress;
          MicroFields(2,2).strain = Results.MicroFieldsM.strain;
          MicroFields(2,2).stress = Results.MicroFieldsM.stress;
        else
            MicroFields = Results.MicroFields;
        end  
    else
        MicroFields = Results.MicroFields;
    end
   
    RUC = props.RUC;
    nb = RUC.NB;
    ng = RUC.NG;
    
    % -- Boundaries between subcells must have 2 points so the plotted
    %    fields can change instantaneously between subcells rather than
    %    being smeared between centroid values
    
    % -- Place an x2-point at bottom and top of each subcell
    h = 0;
    i = 1;
    for b = 1:nb
        x2(i) = h;
        h = h + RUC.h(b);
        i = i + 1;
        x2(i) = h;
        i = i + 1;
    end

    % -- Place an x3-point at left and right of each subcell
    l = 0;
    i = 1;
    for g = 1:ng
        x3(i) = l;
        l = l + RUC.l(g);
        i = i + 1;
        x3(i) = l;
        i = i + 1;
    end
    
    % -- Copy subcell fields to each point in each subcell
    count = 0;
    for j = 1:2*ng
        g = round(j/2);
        for i = 1:2*nb
            b = round(i/2);
            count = count + 1;
            sig11(i,j) = MicroFields(b,g).stress(1);
            sig22(i,j) = MicroFields(b,g).stress(2);
            sig33(i,j) = MicroFields(b,g).stress(3);
            sig23(i,j) = MicroFields(b,g).stress(4);
            sig13(i,j) = MicroFields(b,g).stress(5);
            sig12(i,j) = MicroFields(b,g).stress(6);  
            eps11(i,j) = MicroFields(b,g).strain(1);
            eps22(i,j) = MicroFields(b,g).strain(2);
            eps33(i,j) = MicroFields(b,g).strain(3);
            gam23(i,j) = MicroFields(b,g).strain(4);
            gam13(i,j) = MicroFields(b,g).strain(5);
            gam12(i,j) = MicroFields(b,g).strain(6);  
            mats(i,j) = RUC.mats(b,g);
        end 
    end 
    
    % -- Outlining of material boundaries on/off
    outline = true;
    %outline = false;

    
    % ============ Added for Ch. 7 ============
    % -- Check for subcell failure
    SubcellFailure = false;
    if Results.Type == "bg"
        for b = 1:nb
            for g = 1:ng
                if isfield(Results, 'MoS')
                    if isfield(Results.MoS.Micro(b,g), 'Failed')
                        if Results.MoS.Micro(b,g).Failed
                            SubcellFailure = true;
                        else
                            Results.MoS.Micro(b,g).Failed = false;                    
                        end
                    else
                        Results.MoS.Micro(b,g).Failed = false;
                    end
                end
            end
        end
    end
    
    % -- Turn off outlining if there is subcell failure
    if SubcellFailure
        outline = false;
    end
    % ============ End added for Ch. 7 ============

    
    % -- Create grid
    [X,Y] = meshgrid(x3,x2);

    
    % -- NOTE: For all plots, wwitch shading from interp to faceted 
    %          to plot all subcell boundaries)
   
    % ----------------------------------------------
    % -- Plot RUC for GMC and HFGMC (Chapters 5 - 7) 
    % ----------------------------------------------
    if props.micro == "GMC" || props.micro == "HFGMC"
        figure
        colormap(jet);
        pcolor( X, Y, mats), shading faceted;
        %caxis([0 5]);  % -- Specify colorbar limits for this plot
        c=colorbar; 
        c.FontSize=12;
        c.FontWeight='bold';
        titleNoVi = ['RUC - ',char(num2str(nb)),' by ',char(num2str(ng)), ...
                     ' subcells', ' - Vf = ', char(num2str(RUC.Vf))];         

        if isfield(RUC, 'Vi')
            if RUC.Vi > 0
                title(['RUC - ',char(num2str(nb)),' by ',char(num2str(ng)), ...
                       ' subcells', ' - Vf = ', char(num2str(RUC.Vf)), ...
                       ' - Vi = ',char(num2str(RUC.Vi))],'FontSize',10);
            else
                title(titleNoVi,'FontSize',10);
            end
        else
            title(titleNoVi,'FontSize',10);
        end
        
        xlabel('\bfx_2'); 
        ylabel('\bfx_3','rotation',0);
        set (gca,'FontSize',12,'FontWeight','bold','LineWidth',2)
        axis image;
        axis off;
        if outline
            showOutline(RUC.mats,RUC.h,RUC.l);
        end
        
        % -- Save plot to Word output file or separate .bmp
        if OutInfo.Format == "doc"
           SectTitle = 'Local Fields';
           if ply > 0
               BotTop = ' Bottom of';
               if kk == 2
                   BotTop = ' Top of';
               end
               SectTitle = [SectTitle, ':', BotTop, ' Ply = ', num2str(ply), ', Angle = ', num2str(Angle)];           
           end
           WriteWordFig([pwd,'/', OutInfo.OutFile, '.doc'], SectTitle)
        elseif OutInfo.Format == "txt"
            Plotname = [OutInfo.PlotFile,' - RUC ',AddPly,'.bmp'];
            saveas(gcf,Plotname);
        end
        
    end
    
    % ----------------------
    % -- Plot local stresses
    % ----------------------
    fontsize = SetupSubplot(OutInfo, 'Local Stresses');
    
    %------------------- Sig11 -------------------
    hhh = subplot(2,2,1);
    colormap(jet);
    pcolor( X, Y,sig11), shading interp;
%     caxis([-85.12 173.54]);  % -- Specify colorbar limits for this plot
%     caxis([-40 80]);  % -- Specify colorbar limits for this plot
%    caxis([-50 120]);  % -- Specify colorbar limits for this plot
    c=colorbar; 
    c.FontSize=fontsize;
    c.FontWeight='bold';
    title('Stress 11','FontSize',fontsize);
    xlabel('\bfx_2'); 
    ylabel('\bfx_3','rotation',0);
    set (gca,'FontSize',fontsize,'FontWeight','bold','LineWidth',2)
    axis image;
    axis off;
    set(hhh, 'Units', 'normalized');
    pos = hhh.Position;
    pos(2) = pos(2) - 0.05;
    set(hhh, 'Position', pos);
    
    if outline
        showOutline(RUC.mats,RUC.h,RUC.l);
    end
    
    %------------------- Sig22 -------------------
    hhh = subplot(2,2,2);
    pcolor( X, Y,sig22), shading interp;
%     caxis([-93.62 169.33]); % -- Specify colorbar limits for this plot
%     caxis([0 200]); % -- Specify colorbar limits for this plot
%     caxis([0 300]); % -- Specify colorbar limits for this plot
%     caxis([0 420]); % -- Specify colorbar limits for this plot
%    caxis([-180 110]); % -- Specify colorbar limits for this plot
    c=colorbar; 
    c.FontSize=fontsize;
    c.FontWeight='bold';
    title('Stress 22','FontSize',fontsize);
    xlabel('\bfx_2'); 
    ylabel('\bfx_3','rotation',0);
    set (gca,'FontSize',fontsize,'FontWeight','bold','LineWidth',2)
    axis image;
    axis off;
    set(hhh, 'Units', 'normalized');
    pos = hhh.Position;
    pos(2) = pos(2) - 0.05;
    set(hhh, 'Position', pos);

    if outline
        showOutline(RUC.mats,RUC.h,RUC.l);
    end
    
    %------------------- Sig12 -------------------
    subplot(2,2,3)
    pcolor( X, Y,sig12), shading interp;
    %caxis([-50 80]); % -- Specify colorbar limits for this plot
    c=colorbar; 
    c.FontSize=fontsize;
    c.FontWeight='bold';
    title('Stress 12','FontSize',fontsize);
    xlabel('\bfx_2'); 
    ylabel('\bfx_3','rotation',0);
    set (gca,'FontSize',fontsize,'FontWeight','bold','LineWidth',2)
    axis image;
    axis off;
    
    if outline
        showOutline(RUC.mats,RUC.h,RUC.l);
    end

    %------------------- Sig33 -------------------
    subplot(2,2,4)
    pcolor( X, Y,sig33), shading interp;
%     caxis([-20.54 371.82]); % -- Specify colorbar limits for this plot
%     caxis([-50 70]); % -- Specify colorbar limits for this plot
%     caxis([-80 100]); % -- Specify colorbar limits for this plot
%     caxis([-180 110]); % -- Specify colorbar limits for this plot
%    caxis([0 420]); % -- Specify colorbar limits for this plot
    c=colorbar; 
    c.FontSize=fontsize;
    c.FontWeight='bold';
    title('Stress 33','FontSize',fontsize);
    xlabel('\bfx_2'); 
    ylabel('\bfx_3','rotation',0);
    set (gca,'FontSize',fontsize,'FontWeight','bold','LineWidth',2)
    axis image;
    axis off;
    
    if outline
        showOutline(RUC.mats,RUC.h,RUC.l);
    end

    if abs(Zloc) < 1e-6
        Zloc = 0;
    end
    
    set(gcf,'color','w');
    
   % -- Save plot to Word output file or separate .bmp
   if OutInfo.Format == "doc"
       SectTitle = 'Local Fields';
       if ply > 0
           BotTop = ' - Bottom';
           if kk == 2
               BotTop = ' - Top';
           end
           SectTitle = [SectTitle, ' - Ply: ', num2str(ply), ' - Angle: ', num2str(Angle), BotTop];           
       end
       if props.micro == "GMC" || props.micro == "HFGMC"
            WriteWordFig([pwd,'/', OutInfo.OutFile, '.doc'])
       else
            WriteWordFig([pwd,'/', OutInfo.OutFile, '.doc'], SectTitle)
       end 
   elseif OutInfo.Format == "txt"
       Plotname = [OutInfo.PlotFile,' - sig ',AddPly,'.bmp'];
       saveas(gcf,Plotname);
   end


    % -- Plot sig23 and sig13 for HFGMC (Chapters 6 and 7)
    if (props.micro == "HFGMC")

        fontsize = SetupSubplot(OutInfo, 'Local sig23 and sig13');
        
        if OutInfo.Format == "txt"
           fontsize = 28;
        end
        
        %------------------- Sig23 -------------------
        subplot(1,2,1)
        colormap(jet);
        pcolor( X, Y,sig23), shading interp;
%         caxis([-53.92 54.47]); % -- Specify colorbar limits for this plot
%         caxis([-30 30]); % -- Specify colorbar limits for this plot
%         caxis([-40 50]); % -- Specify colorbar limits for this plot
%        caxis([-120 120]); % -- Specify colorbar limits for this plot
%        caxis([-15 15]); % -- Specify colorbar limits for this plot
        c=colorbar; 
        c.FontSize=fontsize;
        c.FontWeight='bold';
        title('Stress 23','FontSize',fontsize);
        xlabel('\bfx_2'); 
        ylabel('\bfx_3','rotation',0);
        set (gca,'FontSize',fontsize,'FontWeight','bold','LineWidth',2)
        axis image;
        axis off;
        
        if outline
            showOutline(RUC.mats,RUC.h,RUC.l);
        end
        
        %------------------- Sig13 -------------------
        subplot(1,2,2)
        pcolor( X, Y,sig13), shading interp;
        %caxis([-93.62 169.33]); % -- Specify colorbar limits for this plot
        c=colorbar; 
        c.FontSize=fontsize;
        c.FontWeight='bold';
        title('Stress 13','FontSize',fontsize);
        xlabel('\bfx_2'); 
        ylabel('\bfx_3','rotation',0);
        set (gca,'FontSize',fontsize,'FontWeight','bold','LineWidth',2)
        axis image;
        axis off;
        
        if outline
            showOutline(RUC.mats,RUC.h,RUC.l);
        end
        
        set(gcf,'color','w');

        % -- Save plot to Word output file or separate .bmp
        if OutInfo.Format == "doc"
            WriteWordFig([pwd,'/', OutInfo.OutFile, '.doc'])
        elseif OutInfo.Format == "txt"
            Plotname = [OutInfo.PlotFile,' - S23-S13 ',AddPly,'.bmp'];
            saveas(gcf,Plotname);
        end   
        
    end   
   
   
    % ---------------------
    % -- Plot local strains
    % ---------------------
    fontsize = SetupSubplot(OutInfo, 'Local Strains');
    
    %------------------- Eps11 -------------------
    hhh = subplot(2,2,1);
    colormap(jet);
    pcolor( X, Y,eps11), shading interp;
    % -- Stops constant small value from being treated as zero
    caxis([min(min(eps11))-1e-6 max(max(eps11))+1e-6]);
    %caxis([-1e-4 0]);  % -- Specify colorbar limits for this plot
    c=colorbar; 
    c.FontSize=fontsize;
    c.FontWeight='bold';
    title('Strain 11','FontSize',fontsize);
    xlabel('\bfx_2'); 
    ylabel('\bfx_3','rotation',0);
    set (gca,'FontSize',fontsize,'FontWeight','bold','LineWidth',2)
    axis image;
    axis off;
    set(hhh, 'Units', 'normalized');
    pos = hhh.Position;
    pos(2) = pos(2) - 0.05;
    set(hhh, 'Position', pos);
    
    if outline
        showOutline(RUC.mats,RUC.h,RUC.l);
    end

    %------------------- Eps22 -------------------
    hhh = subplot(2,2,2);
    pcolor( X, Y,eps22), shading interp;
    % -- Stops constant small value from being treated as zero
    caxis([min(min(eps22))-1e-6 max(max(eps22))+1e-6]);
    %caxis([0.4e-4 2.2e-4]);  % -- Specify colorbar limits for this plot
    c=colorbar; 
    c.FontSize=fontsize;
    c.FontWeight='bold';
    title('Strain 22','FontSize',fontsize);
    xlabel('\bfx_2'); 
    ylabel('\bfx_3','rotation',0);
    set (gca,'FontSize',fontsize,'FontWeight','bold','LineWidth',2)
    axis image;
    axis off;
    set(hhh, 'Units', 'normalized');
    pos = hhh.Position;
    pos(2) = pos(2) - 0.05;
    set(hhh, 'Position', pos);

    if outline
        showOutline(RUC.mats,RUC.h,RUC.l);
    end

    %------------------- Eps12 -------------------
    subplot(2,2,3)
    pcolor( X, Y, gam12), shading interp;
    % -- Stops constant small value from being treated as zero
    caxis([min(min(gam12))-1e-6 max(max(gam12))+1e-6]);
    %caxis([-1.3e-4 0]); % -- Specify colorbar limits for this plot
    c=colorbar; 
    c.FontSize=fontsize;
    c.FontWeight='bold';
    title('Strain 12','FontSize',fontsize);
    xlabel('\bfx_2'); 
    ylabel('\bfx_3','rotation',0);
    set (gca,'FontSize',fontsize,'FontWeight','bold','LineWidth',2)
    axis image;
    axis off;

    if outline
        showOutline(RUC.mats,RUC.h,RUC.l);
    end

    %------------------- Eps33 -------------------
    subplot(2,2,4)
    pcolor( X, Y,eps33), shading interp;
    % -- Stops constant small value from being treated as zero
    caxis([min(min(eps33))-1e-6 max(max(eps33))+1e-6]);
    %caxis([-1.3e-4 0]);  % -- Specify colorbar limits for this plot
    c=colorbar; 
    c.FontSize=fontsize;
    c.FontWeight='bold';
    title('Strain 33','FontSize',fontsize);
    xlabel('\bfx_2'); 
    ylabel('\bfx_3','rotation',0);
    set (gca,'FontSize',fontsize,'FontWeight','bold','LineWidth',2)
    axis image;
    axis off;

    if outline
        showOutline(RUC.mats,RUC.h,RUC.l);
    end
    
    set(gcf,'color','w');
    
   % -- Save plot to Word output file or separate .bmp
   if OutInfo.Format == "doc"
       WriteWordFig([pwd,'/', OutInfo.OutFile, '.doc'])
   elseif OutInfo.Format == "txt"
       Plotname = [OutInfo.PlotFile,' - eps ',AddPly,'.bmp'];
       saveas(gcf,Plotname);
   end   
   
    % -- Plot gam23 and gam13 for HFGMC  (Chapters 6 and 7)
    if (props.micro == "HFGMC")

        fontsize = SetupSubplot(OutInfo, 'Local gamma23 and gamma13');
        
        if OutInfo.Format == "txt"
            fontsize = 28;
        end
        
        %------------------- Eps23 -------------------
        subplot(1,2,1)
        colormap(jet);
        pcolor( X, Y,gam23), shading interp;
        % -- Stops constant small value from being treated as zero
        caxis([min(min(gam23))-1e-6 max(max(gam23))+1e-6]);
        %caxis([-93.62 169.33]);  % -- Specify colorbar limits for this plot
        c=colorbar; 
        c.FontSize=fontsize;
        c.FontWeight='bold';
        title('Strain 23','FontSize',fontsize);
        xlabel('\bfx_2'); 
        ylabel('\bfx_3','rotation',0);
        set (gca,'FontSize',fontsize,'FontWeight','bold','LineWidth',2)
        axis image;
        axis off;

        if outline
            showOutline(RUC.mats,RUC.h,RUC.l);
        end
        
        %------------------- Eps13 -------------------
        subplot(1,2,2)
        pcolor( X, Y,gam13), shading interp;
        % -- Stops constant small value from being treated as zero
        caxis([min(min(gam13))-1e-6 max(max(gam13))+1e-6]);
        %caxis([-93.62 169.33]);  % -- Specify colorbar limits for this plot
        c=colorbar; 
        c.FontSize=fontsize;
        c.FontWeight='bold';
        title('Strain 13','FontSize',fontsize);
        xlabel('\bfx_2'); 
        ylabel('\bfx_3','rotation',0);
        set (gca,'FontSize',fontsize,'FontWeight','bold','LineWidth',2)
        axis image;
        axis off;

        if outline
            showOutline(RUC.mats,RUC.h,RUC.l);
        end
        
        set(gcf,'color','w');
        
        % -- Save plot to Word output file or separate .bmp
        if OutInfo.Format == "doc"
            WriteWordFig([pwd,'/', OutInfo.OutFile, '.doc'])
        elseif OutInfo.Format == "txt"
            Plotname = [OutInfo.PlotFile,' - gam23-gam13 ',AddPly,'.bmp'];
            saveas(gcf,Plotname);
        end   
        
    end   
   
    %%%============ Added for a Ch5 Problems, and used in Ch 6 & 7 ==============
    % -- Calculate von Mises Stress and pressure
    %PlotSvm = true;
    PlotSvm = false;
    if PlotSvm
        count = 0;
        for j = 1:2*ng
            g = round(j/2);
            for i = 1:2*nb
                b = round(i/2);
                count = count + 1;
                Svm(i,j) = sqrt( (sig11(i,j) - sig22(i,j))^2 + ...
                                 (sig22(i,j) - sig33(i,j))^2 + ...
                                 (sig11(i,j) - sig33(i,j))^2 + ...
                                  6*sig23(i,j)^2 + ...
                                  6*sig13(i,j)^2 + ...
                                  6*sig12(i,j)^2 ) / sqrt(2);
                press(i,j) = -( sig11(i,j) + sig22(i,j) + sig33(i,j) )/3; 
            end 
        end 

        % -- Plot von Mises Stress and pressure
        fontsize = SetupSubplot(OutInfo, 'Local Svm and Pressure');

        if OutInfo.Format == "txt"
            fontsize = 28;
        end

        %------------------- von Mises Stress -------------------
        subplot(1,2,1)
        colormap(jet);
        pcolor( X, Y, Svm), shading interp;
%         caxis([36.46 373.68]);  % -- Specify colorbar limits for this plot
%         caxis([0 200]);  % -- Specify colorbar limits for this plot
%        caxis([0 280]);  % -- Specify colorbar limits for this plot
%        caxis([0 45]);  % -- Specify colorbar limits for this plot
        c=colorbar; 
        c.FontSize=fontsize;
        c.FontWeight='bold';
        title('von Mises stress','FontSize',fontsize);
        xlabel('\bfx_2'); 
        ylabel('\bfx_3','rotation',0);
        set (gca,'FontSize',fontsize,'FontWeight','bold','LineWidth',2)
        axis image;
        axis off;

        if outline
            showOutline(RUC.mats,RUC.h,RUC.l);
        end

        %------------------- Pressure -------------------
        subplot(1,2,2)
        pcolor( X, Y, press), shading interp;
%         caxis([-233.15 44.39]); % -- Specify colorbar limits for this plot
%         caxis([-100 10]); % -- Specify colorbar limits for this plot
%        caxis([-170 20]); % -- Specify colorbar limits for this plot
%        caxis([-25 25]); % -- Specify colorbar limits for this plot
        c=colorbar; 
        c.FontSize=fontsize;
        c.FontWeight='bold';
        title('Pressure','FontSize',fontsize);
        xlabel('\bfx_2'); 
        ylabel('\bfx_3','rotation',0);
        set (gca,'FontSize',fontsize,'FontWeight','bold','LineWidth',2)
        axis image;
        axis off;

        if outline
             showOutline(RUC.mats,RUC.h,RUC.l);
        end

        set(gcf,'color','w');

        % -- Save plot to Word output file or separate .bmp
        if OutInfo.Format == "doc"
            WriteWordFig([pwd,'/', OutInfo.OutFile, '.doc'])
        elseif OutInfo.Format == "txt"
            Plotname = [OutInfo.PlotFile,' - Svm-press ',AddPly,'.bmp'];
            saveas(gcf,Plotname);
        end   

        % -- Calculate the avg von Mises stress & pressure (Ch5 Exercise 5.9)
        % User to add this code for Exercise 5.9


    end
    %============  End added for a Ch5 Problems ==============
    
    %============ Added for Ch7 ==============
    if SubcellFailure
       % -- Plot failed subcells (material num > 100)

        figure('units','normalized','outerposition',[0 0 1 1])

        %------------------- Failure -------------------

        colormap(jet);
        pcolor( X, Y, mats), shading faceted;
        if min(min(mats)) >= 15
            matsmin = 14;
        else
            matsmin = min(min(mats));
        end
        caxis([matsmin 16]);
        %c=colorbar; 
        %c.FontSize=fontsize;
        %c.FontWeight='bold';
        title('Subcell Failure','FontSize',fontsize);
        xlabel('\bfx_2'); 
        ylabel('\bfx_3','rotation',0);
        set (gca,'FontSize',fontsize,'FontWeight','bold','LineWidth',2)
        axis image;
        axis off;    

        % -- Save plot to Word output file or separate .bmp
        if OutInfo.Format == "doc"
            WriteWordFig([pwd,'/', OutInfo.OutFile, '.doc'])
        elseif OutInfo.Format == "txt"
            Plotname = [OutInfo.PlotFile,' - Failure ','.bmp'];
            saveas(gcf,Plotname);
        end   
    end
    %============ End added for Ch7 ==============

end 

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function [fontsize] = SetupSubplot(OutInfo, title)

if OutInfo.Format == "doc"
    fontsize = 12;
    figure
    dim = [.4 .7 .3 .3];
else
    figure('units','normalized','outerposition',[0 0 1 1])
    fontsize = 18;
    dim = [.45 .7 .3 .3];
end
a = annotation('textbox',dim,'String',title,'FitBoxToText','on');    
a.FontSize = fontsize - 2;
a.FontWeight = 'bold';
a.LineStyle = 'none';
    
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
    
function []= showOutline(RUC,H,L )

% -- Outlining interface and fibers inside RUC

[m,n]=size(RUC);
px=0;py=0;
hold on;
xlim([0 sum(L)]);ylim([0 sum(H)]);
daspect([1,1,1]);
LW =1.;   % Inside border line width dimension
BW =1.5; % Outside RUC border line width dimension
% Initialize arrays
X=zeros(10000,1);Y=zeros(10000,1);
% First horizontal Lines
Jcntrl=0;
for i = 2: m
    for j=2:n+1
        Jcntrl=Jcntrl+1;
        if (RUC(i-1,j-1) == RUC(i,j-1))
            X(Jcntrl)=NaN;Y(Jcntrl)=NaN;            
        else
            X(Jcntrl)= px ; X(Jcntrl+1) = px+L(j-1);
            Y(Jcntrl)= py+H(i-1);Y(Jcntrl+1) = py+H(i-1);
            Jcntrl=Jcntrl+1;
        end
        px=px+L(j-1);      
    end
    px=0;
    Jcntrl=Jcntrl+1;
    X(Jcntrl)=NaN;Y(Jcntrl)=NaN;
    py=py+(H(i-1));
end
plot(X,Y,'k','LineWidth',LW)
% Second Vertical Lines
px=0;py=0;
for j = 2: n
    for i = 1:m
        Jcntrl = Jcntrl+1;
        if (RUC(i,j-1) == RUC(i,j))
            X(Jcntrl)=NaN;Y(Jcntrl)=NaN;           
        else
            X(Jcntrl)= px+L(j-1) ; X(Jcntrl+1) = px +L(j-1);
            Y(Jcntrl)= py ;Y(Jcntrl+1) = py+H(i);
            Jcntrl=Jcntrl+1;
        end
        py=py+H(i);        
    end    
    py=0;
    Jcntrl = Jcntrl+1;
    X(Jcntrl)=NaN;Y(Jcntrl)=NaN;
    px=px+ L(j-1);
end
X=X(1:Jcntrl);Y=Y(1:Jcntrl);
Idx = ~isnan(X); % keep all non-nans
Idx(logical([1;diff(~Idx)])) = true; % keep first of consecutive nans
X = X(Idx);                          % remove unwanted values from x
Y = Y(Idx);                          % remove unwanted values from y
plot(X,Y,'LineWidth',LW,'Color','k')
axis off;
% Draw RUC Border
% rectangle('Position',[0 0 sum(L) sum(H)],'EdgeColor','k',...
%     'LineWidth',BW);
end