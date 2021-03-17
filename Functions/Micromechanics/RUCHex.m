function[fvr,ivr,RUC]= RUCHex(Fvr,Ivr,Ni,Nw,Plot)
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

% Program that generates a Hexagonal packing CMC Ruc with 
% interface.
% Fiber radius and Interface thickness are user difinable.
% Developed by Pappu L.N. Murthy, NASA GRC
% Ref: Murthy, P.L.N., Bednarcyk, B.A., and Mital, S.K. (2017) “A Compilation of MATLAB 
%      Scripts and Functions for MAC/GMC Analysis” NASA/TM-2017-219500.

% format long;clear all;
if (nargin==3)
    Ni=1;
end
if (Ni == 0)
    Ivr = 0;
end
if (Ivr == 0)
    Ni = 0;
end
M=2;F=1;I=3;
W=1;H=W/sqrt(3);
TA = W*H;
FA = Fvr*TA;
Radf = sqrt(FA/(2*pi));
IA = Ivr*TA;
Radi = sqrt((FA+IA)/(2*pi));
flag=0;
% Determine if fiber radius and interface radius are same
if ( abs(Radf-Radi)>=eps )
    flag=1;
end
if(flag)
    ti=Radi-Radf;
    Radi=round(Ni*Radi/ti);
    Radf=round(Ni*Radf/ti);
    ny=round(Ni*W/ti);
else
    ny=Nw;
    Radf= round( sqrt ( ((ny*ny/sqrt(3))*Fvr)/(2*pi)));
    Radi=Radf;
end
if (mod(ny,2))
    ny=ny-1;
end
nx=fix(1.732\ny);
if(mod(nx,2))
    nx=nx-1;
end
Lbase=ones(1,ny);Hbase=ones(1,nx);
base = M*ones(nx/2,ny/2);
% Rectangular RUC with two of the opposite corners made of 
% fibers with interface.
% Fiber volume ratio Optimization
Fact=(.95:0.001:1.2);
fv=zeros(1,length(Fact));
iv=fv;Raf=fv;Rai=fv;
for ii = 1:length(Fact)
    Rf=Fact(ii)*Radf;Ri=Fact(ii)*Radi;
    for i=1:nx/2
        for j=1:ny/2
            dist=sqrt(i^2+j^2);
            dist2=sqrt((nx/2-i+1)^2 + (ny/2-j+1)^2);
            if (dist<=Rf || dist2<=Rf)
                base(i,j)=F;
            elseif (dist>Rf && dist<=Ri || dist2>Rf && dist2<=Ri)
                base(i,j)=I;
            end
        end
    end
    Ma= base==M;Marea=sum(Ma(:));
    Fa= base==F;Farea=sum(Fa(:));
    Ia= base==I;Iarea=sum(Ia(:));
    Ta=Marea+Farea+Iarea;
    fvr= Farea/Ta;ivr=Iarea/Ta;
    fv(ii)=fvr;
    iv(ii)=ivr;
    Raf(ii)=Rf;Rai(ii)=Ri;
    if (fv(ii)>=Fvr)
        break
    end
end
Rf=Raf(ii);
if (ii>1)
    if ( abs(Fvr-fv(ii)) > abs(Fvr-fv(ii-1)) )
        Rf=Raf(ii-1);
    end
end

if (flag)
    % Interphase Optimzation
    base(base(:)==I)=M;
    Fact=(Rf/Radi:0.001:1.2);
    for ii = 1:length(Fact)
        Ri=Fact(ii)*Radi;
        for i=1:nx/2
            for j=1:ny/2
                dist=sqrt(i^2+j^2);
                dist2=sqrt((nx/2-i+1)^2 + (ny/2-j+1)^2);
                if (dist<=Rf || dist2<=Rf)
                    base(i,j)=F;
                elseif (dist>Rf && dist<=Ri || dist2>Rf && dist2<=Ri)
                    base(i,j)=I;
                end
            end
        end
        base2=flipud(base);base12=[base;base2];
        base13=fliplr(base12);
        base14 = [base12,base13];
        Ma= base14==M;Marea=sum(Ma(:));
        Fa= base14==F;Farea=sum(Fa(:));
        Ia= base14==I;Iarea=sum(Ia(:));
        Ta=Marea+Farea+Iarea;
        fvr= Farea/Ta; ivr=Iarea/Ta;
        Rai(ii)=Ri;
        fv(ii)=fvr;
        iv(ii)=ivr;
        if (iv(ii)>=Ivr)
            break
        end
    end
    Ri=Rai(ii);
    if (ii>1)
        if ( abs(Ivr-iv(ii)) > abs(Ivr-iv(ii-1)) )
            Ri=Rai(ii-1);
        end
    end
    
    %----------------------------------------------------------
    % Generate current Ruc;
    base = M*ones(nx/2,ny/2); % Rectangular RUC with two of the opposite
    for i=1:nx/2
        for j=1:ny/2
            dist=sqrt(i^2+j^2);
            dist2=sqrt((nx/2-i+1)^2 + (ny/2-j+1)^2);
            if (dist<=Rf || dist2<=Rf)
                base(i,j)=F;
            elseif (dist>Rf && dist<=Ri || dist2>Rf && dist2<=Ri)
                base(i,j)=I;
            end
        end
    end
    base2=flipud(base);base12=[base;base2];
    base13=fliplr(base12);
    base14 = [base12,base13];
    Ma= base14==M;Marea=sum(Ma(:));
    Fa= base14==F;Farea=sum(Fa(:));
    Ia= base14==I;Iarea=sum(Ia(:));
    Ta=Marea+Farea+Iarea;
    fvr= Farea/Ta; mvr=Marea/Ta; ivr=Iarea/Ta;
%----------------------------------------------------------
else
    % PMC. No interface.
    base = M*ones(nx/2,ny/2); % Rectangular RUC with two of the opposite
    for i=1:nx/2
        for j=1:ny/2
            dist=sqrt(i^2+j^2);
            dist2=sqrt((nx/2-i+1)^2 + (ny/2-j+1)^2);
            if (dist<=Rf || dist2<=Rf)
                base(i,j)=F;
            end
        end
    end
    base2=flipud(base);base12=[base;base2];
    base13=fliplr(base12);
    base14 = [base12,base13];
    Ma= base14==M;Marea=sum(Ma(:));
    Fa= base14==F;Farea=sum(Fa(:));
    Ta=Marea+Farea;
    fvr= Farea/Ta; mvr=Marea/Ta;
end
RUC = base14;

H = ones(1,nx);
L = ones(1,ny);

% -- plot
if (Plot)
    [m,n]=size(RUC); Nmat = length(unique(RUC));
    px=0;py=0;
    switch Nmat
        case 1
            cmap=[0 0.5 1];
        case 2
            cmap=[0 0.5 1;1 1 0.5];
        case 3
            cmap= [0 0.5 1;1 1 0.5;0 0 0];
        case 4
            cmap= [0 0.5 1;1 1 0.5;0 0 0;0 1 1];
    end
    for i = 1:m
        for j=1:n
            rectangle('Position',[px py L(j) H(i)],'FaceColor',...
            cmap(RUC( i ,j),:));
            px=px+L (j);
        end
        px=0;
        py=py+H (i);
    end
    xlim([0 sum(L )]);ylim([0 sum(H )]);
    daspect([1,1,1]);
    txt = [num2str(m), 'x' , num2str(n),' Hex RUC, Vf = ', num2str(fvr), ', Vi = ', num2str(ivr)];
    title(txt);
    axis off;
    drawnow
%     dbox = warndlg('Click OK when finished inspecting RUC','RUC Inspection');
%     uiwait(dbox);
end



