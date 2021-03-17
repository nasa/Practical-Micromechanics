function [Cstar, CTEs, As] = HFGMC(Constits, RUC)
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
% Purpose: Calculate the effective stiffness, effective CTES, and the mechanical and 
%          thermal strain concentration tensors for the HFGMC micromechanics theory 
% Input:
% - Constits: Struct containing constituent material information
% - RUC: Struct containing repeating unit cell (RUC) information
% Output:
% - Cstar: Effective stiffness matrix
% - CTEs: Effective CTEs
% - As: Struct containing subcell concentration tensors
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    NB = RUC.NB;
    NG = RUC.NG;

    % -- Output run time for large RUCs
    if NB*NG > 25*25
        disp(['  ** Getting Effective Properties using HFGMC - RUC = ', num2str(NB), ' by ', num2str(NG)])
        tic
    end
    
    % -- Preallocate props
    props(NB, NG).C = zeros(6);
    
    % -- Extract and store needed consituent properties 
    for IB = 1:NB
        for IG = 1:NG
            switch(RUC.matsCh(IB, IG))
                case 'F'
                    constit = Constits.Fiber;
                case '1'
                    constit = Constits.FiberDam;
                case 'M'
                    constit = Constits.Matrix;
                case '2'
                    constit = Constits.MatrixDam;
                case 'I'
                    constit = Constits.Interface;
                case '3'
                    constit = Constits.InterfaceDam;
                otherwise
                    error(['Invalid RUC.matsch in HFGMC - RUC.matsch = ', RUC.matsch])
            end
            props(IB, IG).C = GetCFromProps(constit);
            props(IB, IG).alpha = [constit.aL; constit.aT; constit.aT; 0; 0; 0];
        end
    end
    
    % -- Subcell and RUC dimensions
    h = RUC.h;
    l = RUC.l;
    sumH = sum(h);
    sumL = sum(l);
    
    % -- Determine the local stiffness matrix K (Eq. 6.47) and the Dij
    %    terms (Eqs. 6.41 - 6.46) per subcell
    [ak, dd] = local(RUC, props);
    
    % -- Determine the strain concentration matrics [As] of the subcells
    
    % -- Pre-allocate
    As(NB, NG).A = zeros(6);
    Omega = 0;
    
    % -- Successively apply 6 global strain components (epsb)
    for i = 1:6
        
        % -- Apply epsb(i) = 1, all other epsb = 0
        dtemp = 0;
        epsb = zeros(6,1);
        if i < 4
            epsb(i) = 1;
        else
            % -- epsb shear strain is tensorial - Apply gamma = 1 to get A_Eng
            epsb(i) = 0.5;           
        end

        % -- Flag for reassembling the global structural stiffness matrix, Omega
        %    This only needs to be done once as changing epsb only affects the
        %    force vector f on the right hand side of Eq. 6.61
        assemble = true;
        if i > 1
            assemble = false;
        end
                
        % -- Assemble the global structural stiffness matrix, Omega, the
        %    force vector, f, and solve for the unknown surf average displacement 
        %    vector, U (Eq. 6.61) 
        [props, Omega, U] = global_Om(assemble, props, RUC, ak, Omega, epsb, dtemp);

        for j = 1:NB
            for k = 1:NG

              % -- Determine the microvariables (from Eqs. 6.33 - 6.36 and 6.40)
              w = micro_variables(j, k, NB, RUC, dd, U);

              % -- Determine column i of current subcell's strain concentration tensor
              %    Does not yet include the global strain contributions, see Eqs. 6.5 - 6.10
              As(j, k).A(1,i) = 0;
              As(j, k).A(2,i) = w{2}.w10;
              As(j, k).A(3,i) = w{3}.w01;
              As(j, k).A(4,i) = w{2}.w01 + w{3}.w10;
              As(j, k).A(5,i) = w{1}.w01;
              As(j, k).A(6,i) = w{1}.w10;

            end 
        end

    end
    
    % -- The concentration matrices of the entire field = [I + As] (see Eqs. 6.5 - 6.10)
    %    Adding in the identity matrix now includes the global strains that appear on the
    %    right hand sides of these equations
    for j = 1:NB
       for k = 1:NG
           As(j, k).A = As(j, k).A + eye(6);
       end
    end
  
    % -- Calculate effective stiffness (Eq. 6.67)
    Cstar = zeros(6);
    for IG = 1:NG
        for IB = 1:NB
            Cstar = Cstar + RUC.h(IB)*RUC.l(IG)*props(IB,IG).C*As(IB, IG).A/(sumH*sumL);
        end
    end
    
    % ------------ THERMAL EFFECTS --------------
    
    %  -- Calculate effective compliance matrix SS
    SS = inv(Cstar);

    % -- Determine the stress concentration matrix B and the transpose BT 
    for IG = 1:NG
        for IB = 1:NB
            B(IB, IG).B = props(IB, IG).C*As(IB, IG).A*SS; % -- Eq. 6.70 
            B(IB, IG).BT = transpose(B(IB, IG).B);
        end
    end

    % -- Calculate the effective CTEs (Eq. 6.69)
    CTEs = zeros(6,1);
    for IG = 1:NG
        for IB = 1:NB
            CTEs = CTEs + h(IB)*l(IG)*B(IB, IG).BT*props(IB, IG).alpha/(sumH*sumL);
        end
    end

    % -- Determine the thermal strain concentration tensor AT
    % -- Apply epsb(i) = 0, dtemp = 1
    dtemp = 1;
    epsb = zeros(6,1);

    % -- Omega does not change; epsb and dtemp only affect f, see Eq. 6.72
    assemble = false;
    
    % -- Assemble f and solve for U
    [props, ~, U] = global_Om(assemble, props, RUC, ak, Omega, epsb, dtemp);

    CTE_check = 0;
    for j = 1:NB
        for k = 1:NG

          w = micro_variables(j, k, NB, RUC, dd, U);

          As(j, k).AT(1,1) = 0;
          As(j, k).AT(2,1) = w{2}.w10;
          As(j, k).AT(3,1) = w{3}.w01;
          As(j, k).AT(4,1) = w{2}.w01 + w{3}.w10;
          As(j, k).AT(5,1) = w{1}.w01;
          As(j, k).AT(6,1) = w{1}.w10;

          % -- CTE_Check (Eq. 6.79) should match CTEs
          CTE_check = CTE_check + -(Cstar^(-1))*props(j, k).C*(As(j, k).AT ...
                      - props(j, k).alpha)*h(j)*l(k)/(sumH*sumL);
          
          % -- Store subcell stiffness and CTE per subcell
          As(j, k).C = props(j,k).C;
          As(j, k).alpha = props(j,k).alpha;   
          
        end 
    end
    
    % -- Output run time for large RUCs
    if NB*NG > 25*25
        ElapsedTime = toc;
        disp(['     HFGMC homogenization time = ', char(num2str(ElapsedTime)), ' seconds']);
    end
    
end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [props, Omega, U] = global_Om(assemble, props, RUC, ak, Omega, epsb, dtemp)

% -- Assemble the global structural stiffness matrix, Omega, the force vector, f, 
%    and solve for the unknown surface average displacement vector, U (Eq. 6.61) 

global NB NG
pinned = false; % -- Flag for pinned RUC corners.  HFGMC originally used this.

NB  =  RUC.NB;
NG  =  RUC.NG;
NEQ = 6*NB*NG + 3*(NB + NG);
NNZ = 8*NB*NG + 20; % -- Estimate of upper bound of non-zeros in A

% -- Calculate thermal stress terms per subcell
for j = 1:NB
    for k = 1:NG
        props(j,k).gama = props(j,k).C*props(j,k).alpha;
    end
end


% -- The global matrix A is stored in sparse format, storing only non-zeros
%    NR and NC are the row and column of each non-zero element in A

% -- Pre-define A, NR, and NC for speed of assembly
A = zeros(NNZ, 1);
NR = zeros(NNZ, 1);
NC = zeros(NNZ, 1);

% -- LHS - only assemble the first time (assemble  =  true)
if assemble

    ic = 0; % -- Counter of actual number of non-zeros in A
    n  = 0;

    for  j = 1:NB - 1
        for  k = 1:NG

            nsub  = j   + NB*(k - 1);
            nsub1 = j + 1 + NB*(k - 1);

            % -- s2m( + )(beta,gama) = s2m( - )(beta + 1,gama), Eq. 6.49

            for  m = 1:3
                
                m1 = 2*m - 1;
                m2 = 2*m;

                n = n + 1;

                nr1 = n;
                nc1 = nu(j + 1,k) + 1;
                a1 = ak(m1,1,nsub) - ak(m2,2,nsub1);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j + 2,k) + 1;
                a1 =  -ak(m2,1,nsub1);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j,k) + 1;
                a1 = ak(m1,2,nsub);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j + 1,k) + 2;
                a1 = ak(m1,3,nsub) - ak(m2,4,nsub1);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j + 2,k) + 2;
                a1 =  -ak(m2,3,nsub1);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j,k) + 2;
                a1 = ak(m1,4,nsub);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j + 1,k) + 3;
                a1 = ak(m1,5,nsub) - ak(m2,6,nsub1);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j + 2,k) + 3;
                a1 =  -ak(m2,5,nsub1);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j,k) + 3;
                a1 = ak(m1,6,nsub);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j,k + 1) + 4;
                a1 = ak(m1,7,nsub);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j + 1,k + 1) + 4;
                a1 =  -ak(m2,7,nsub1);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j,k) + 4;
                a1 = ak(m1,8,nsub);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j + 1,k) + 4;
                a1 =  -ak(m2,8,nsub1);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j,k + 1) + 5;
                a1 = ak(m1,9,nsub);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j + 1,k + 1) + 5;
                a1 =  -ak(m2,9,nsub1);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j,k) + 5;
                a1 = ak(m1,10,nsub);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j + 1,k) + 5;
                a1 =  -ak(m2,10,nsub1);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j,k + 1) + 6;
                a1 = ak(m1,11,nsub);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j + 1,k + 1) + 6;
                a1 =  -ak(m2,11,nsub1);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j,k) + 6;
                a1 = ak(m1,12,nsub);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j + 1,k) + 6;
                a1 =  -ak(m2,12,nsub1);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);
               
            end
            
        end
    end

    %**********************************************************************

    for  k = 1:NG - 1
        for  j = 1:NB

            nsub  = j + NB*(k - 1);
            nsub1 = j + NB* k;

            % -- s3j( + )(beta,gama)  =  s3j( - )(beta,gama + 1), Eq. 6.51

            for  m = 4:6
                
                m1 = 2*m - 1;
                m2 = 2*m;

                n = n + 1;

                nr1 = n;
                nc1 = nu(j + 1,k) + 1;
                a1 = ak(m1,1,nsub);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j + 1,k + 1) + 1;
                a1 =  -ak(m2,1,nsub1);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j,k) + 1;
                a1 = ak(m1,2,nsub);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j,k + 1) + 1;
                a1 =  -ak(m2,2,nsub1);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j + 1,k) + 2;
                a1 = ak(m1,3,nsub);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j + 1,k + 1) + 2;
                a1 =  -ak(m2,3,nsub1);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j,k) + 2;
                a1 = ak(m1,4,nsub);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j,k + 1) + 2;
                a1 =  -ak(m2,4,nsub1);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j + 1,k) + 3;
                a1 = ak(m1,5,nsub);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j + 1,k + 1) + 3;
                a1 =  -ak(m2,5,nsub1);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j,k) + 3;
                a1 = ak(m1,6,nsub);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j,k + 1) + 3;
                a1 =  -ak(m2,6,nsub1);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j,k + 1) + 4;
                a1 = ak(m1,7,nsub) - ak(m2,8,nsub1);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j,k + 2) + 4;
                a1 =  -ak(m2,7,nsub1);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j,k) + 4;
                a1 = ak(m1,8,nsub);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j,k + 1) + 5;
                a1 = ak(m1,9,nsub) - ak(m2,10,nsub1);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j,k + 2) + 5;
                a1 =  -ak(m2,9,nsub1);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j,k) + 5;
                a1 = ak(m1,10,nsub);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j,k + 1) + 6;
                a1 = ak(m1,11,nsub) - ak(m2,12,nsub1);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j,k + 2) + 6;
                a1 =  -ak(m2,11,nsub1);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

                nr1 = n;
                nc1 = nu(j,k) + 6;
                a1 = ak(m1,12,nsub);
                [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            end
        
        end
    end
    
    %**********************************************************************

    % -- apply periodicity conditions of ui at x2 = 0 & x2 = h, Eq. 6.52
    j1 = 1;
    j2 = NB;

    for k = 1:NG % -- Start long k loop

        % -- Pinned corners
        if ( (pinned ~= 0 && k == 1) || (pinned ~= 0 && k == NG) )

            % -- u1( - )(1,gama) = 0;
            n = n + 1;
            nr1 = n;
            nc1 = nu(j1,k) + 1;
            a1 = 1;
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            % -- u2( - )(1,gama)  =  0;
            n = n + 1;
            nr1 = n;
            nc1 = nu(j1,k) + 2;
            a1 = 1;
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            % -- u3( - )(1,gama)  =  0;
            n = n + 1;
            nr1 = n;
            nc1 = nu(j1,k) + 3;
            a1 = 1;
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            % -- u1( + )(n_beta,gama)  =  0;
            n = n + 1;
            nr1 = n;
            nc1 = nu(j2 + 1,k) + 1;
            a1 = 1;
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            % -- u2( + )(n_beta,gama)  =  0;
            n = n + 1;
            nr1 = n;
            nc1 = nu(j2 + 1,k) + 2;
            a1 = 1;
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            % -- u3( + )(n_beta,gama)  =  0;
            n = n + 1;
            nr1 = n;
            nc1 = nu(j2 + 1,k) + 3;
            a1 = 1;
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

        % -- Non-pinned corners
        else 
            
            % -- u1( - )(1,gama)  =  u1( + )(n_beta,gama);
            n = n + 1;
            nr1 = n;
            nc1 = nu(j1,k) + 1;
            a1 = 1;
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2 + 1,k) + 1;
            a1 =  - 1;
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            % -- u2( - )(1,gama)  =  u2( + )(n_beta,gama);
            n = n + 1;
            nr1 = n;
            nc1 = nu(j1,k) + 2;
            a1 = 1;
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2 + 1,k) + 2;
            a1 =  - 1;
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            % -- u3( - )(1,gama)  =  u3( + )(n_beta,gama);
            n = n + 1;
            nr1 = n;
            nc1 = nu(j1,k) + 3;
            a1 = 1;
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2 + 1,k) + 3;
            a1 =  - 1;
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            % -- Apply periodicity conditions of s2j at x2 = 0 & x2 = H, Eq. 6.53
            nsub1 = j1 + NB*(k - 1);
            nsub2 = j2 + NB*(k - 1);

            % -- s21( - )(1,gama)  =  s21( + )(n_beta,gama);
            n = n + 1;
            nr1 = n;
            nc1 = nu(j1 + 1,k) + 1;
            a1 = ak(2,1,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2 + 1,k) + 1;
            a1 =  -ak(1,1,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1,k) + 1;
            a1 = ak(2,2,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2,k) + 1;
            a1 =  -ak(1,2,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1 + 1,k) + 2;
            a1 = ak(2,3,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2 + 1,k) + 2;
            a1 =  -ak(1,3,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1,k) + 2;
            a1 = ak(2,4,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2,k) + 2;
            a1 =  -ak(1,4,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1 + 1,k) + 3;
            a1 = ak(2,5,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2 + 1,k) + 3;
            a1 =  -ak(1,5,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1,k) + 3;
            a1 = ak(2,6,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2,k) + 3;
            a1 =  -ak(1,6,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1,k + 1) + 4;
            a1 = ak(2,7,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2,k + 1) + 4;
            a1 =  -ak(1,7,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1,k) + 4;
            a1 = ak(2,8,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2,k) + 4;
            a1 =  -ak(1,8,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1,k + 1) + 5;
            a1 = ak(2,9,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2,k + 1) + 5;
            a1 =  -ak(1,9,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1,k) + 5;
            a1 = ak(2,10,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2,k) + 5;
            a1 =  -ak(1,10,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1,k + 1) + 6;
            a1 = ak(2,11,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2,k + 1) + 6;
            a1 =  -ak(1,11,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1,k) + 6;
            a1 = ak(2,12,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2,k) + 6;
            a1 =  -ak(1,12,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            % -- s22( - )(1,gama)  =  s22( + )(n_beta,gama);

            n = n + 1;
            nr1 = n;
            nc1 = nu(j1 + 1,k) + 1;
            a1 = ak(4,1,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2 + 1,k) + 1;
            a1 =  -ak(3,1,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1,k) + 1;
            a1 = ak(4,2,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2,k) + 1;
            a1 =  -ak(3,2,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1 + 1,k) + 2;
            a1 = ak(4,3,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2 + 1,k) + 2;
            a1 =  -ak(3,3,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1,k) + 2;
            a1 = ak(4,4,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2,k) + 2;
            a1 =  -ak(3,4,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1 + 1,k) + 3;
            a1 = ak(4,5,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2 + 1,k) + 3;
            a1 =  -ak(3,5,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1,k) + 3;
            a1 = ak(4,6,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2,k) + 3;
            a1 =  -ak(3,6,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1,k + 1) + 4;
            a1 = ak(4,7,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2,k + 1) + 4;
            a1 =  -ak(3,7,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1,k) + 4;
            a1 = ak(4,8,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2,k) + 4;
            a1 =  -ak(3,8,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1,k + 1) + 5;
            a1 = ak(4,9,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2,k + 1) + 5;
            a1 =  -ak(3,9,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1,k) + 5;
            a1 = ak(4,10,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2,k) + 5;
            a1 =  -ak(3,10,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1,k + 1) + 6;
            a1 = ak(4,11,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2,k + 1) + 6;
            a1 =  -ak(3,11,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1,k) + 6;
            a1 = ak(4,12,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2,k) + 6;
            a1 =  -ak(3,12,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            % -- s23( - )(1,gama)  =  s23( + )(n_beta,gama);
            n = n + 1;
            nr1 = n;
            nc1 = nu(j1 + 1,k) + 1;
            a1 = ak(6,1,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2 + 1,k) + 1;
            a1 =  -ak(5,1,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1,k) + 1;
            a1 = ak(6,2,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2,k) + 1;
            a1 =  -ak(5,2,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1 + 1,k) + 2;
            a1 = ak(6,3,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2 + 1,k) + 2;
            a1 =  -ak(5,3,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1,k) + 2;
            a1 = ak(6,4,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2,k) + 2;
            a1 =  -ak(5,4,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1 + 1,k) + 3;
            a1 = ak(6,5,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2 + 1,k) + 3;
            a1 =  -ak(5,5,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1,k) + 3;
            a1 = ak(6,6,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2,k) + 3;
            a1 =  -ak(5,6,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1,k + 1) + 4;
            a1 = ak(6,7,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2,k + 1) + 4;
            a1 =  -ak(5,7,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1,k) + 4;
            a1 = ak(6,8,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2,k) + 4;
            a1 =  -ak(5,8,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1,k + 1) + 5;
            a1 = ak(6,9,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2,k + 1) + 5;
            a1 =  -ak(5,9,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1,k) + 5;
            a1 = ak(6,10,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2,k) + 5;
            a1 =  -ak(5,10,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1,k + 1) + 6;
            a1 = ak(6,11,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2,k + 1) + 6;
            a1 =  -ak(5,11,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j1,k) + 6;
            a1 = ak(6,12,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j2,k) + 6;
            a1 =  -ak(5,12,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);
            
        end

    end % -- End long k loop

    %--------------------------------------------------------------------

    % -- apply periodicity conditions of ui at x3  =  0 & x3 = L, Eq. 6.54
    k1 = 1;
    k2 = NG;

    for  j = 1:NB % -- Start long j loop

        % -- Pinned corners
        if ((pinned ~= 0 && j == 1) || (pinned ~= 0 && j == NB) )

            % -- u1( - )(beta,1)  =  0;
            n = n + 1;
            nr1 = n;
            nc1 = nu(j,k1) + 4;
            a1 = 1;
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            % -- u2( - )(beta,1)  =  0;
            n = n + 1;
            nr1 = n;
            nc1 = nu(j,k1) + 5;
            a1 = 1;
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            % -- u3( - )(beta,1)  =  0;
            n = n + 1;
            nr1 = n;
            nc1 = nu(j,k1) + 6;
            a1 = 1;
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            % -- u1( + )(beta,n_gama)  =  0;
            n = n + 1;
            nr1 = n;
            nc1 = nu(j,k2 + 1) + 4;
            a1 = 1;
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            % -- u2( + )(beta,n_gama)  =  0;
            n = n + 1;
            nr1 = n;
            nc1 = nu(j,k2 + 1) + 5;
            a1 = 1;
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            % -- u3( + )(beta,n_gama)  =  0;
            n = n + 1;
            nr1 = n;
            nc1 = nu(j,k2 + 1) + 6;
            a1 = 1;
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

        % -- Non-pinned corners
        else

            % -- u1( - )(beta,1)  =  u1( + )(beta,n_gama);
            n = n + 1;
            nr1 = n;
            nc1 = nu(j,k1) + 4;
            a1 = 1;
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2 + 1) + 4;
            a1 =  - 1;
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            % -- u2( - )(beta,1)  =  u2( + )(beta,n_gama);
            n = n + 1;
            nr1 = n;
            nc1 = nu(j,k1) + 5;
            a1 = 1;
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2 + 1) + 5;
            a1 =  -1;
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            % -- u3( - )(beta,1)  =  u3( + )(beta,n_gama);
            n = n + 1;
            nr1 = n;
            nc1 = nu(j,k1) + 6;
            a1 = 1;
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2 + 1) + 6;
            a1 =  -1;
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            % -- Apply periodicity conditions of s3j at x3 = 0 & x3  = L, Eq. 6.55
            nsub1 = j + NB*(k1 - 1);
            nsub2 = j + NB*(k2 - 1);

            % -- s31( - )(beta,1)  =  s31( + )(beta,n_gama);
            n = n + 1;
            nr1 = n;
            nc1 = nu(j + 1,k1) + 1;
            a1 = ak(8,1,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j + 1,k2) + 1;
            a1 =  -ak(7,1,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k1) + 1;
            a1 = ak(8,2,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2) + 1;
            a1 =  -ak(7,2,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j + 1,k1) + 2;
            a1 = ak(8,3,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j + 1,k2) + 2;
            a1 =  -ak(7,3,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k1) + 2;
            a1 = ak(8,4,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2) + 2;
            a1 =  -ak(7,4,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j + 1,k1) + 3;
            a1 = ak(8,5,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j + 1,k2) + 3;
            a1 =  -ak(7,5,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k1) + 3;
            a1 = ak(8,6,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2) + 3;
            a1 =  -ak(7,6,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k1 + 1) + 4;
            a1 = ak(8,7,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2 + 1) + 4;
            a1 =  -ak(7,7,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k1) + 4;
            a1 = ak(8,8,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2) + 4;
            a1 =  -ak(7,8,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k1 + 1) + 5;
            a1 = ak(8,9,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2 + 1) + 5;
            a1 =  -ak(7,9,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k1) + 5;
            a1 = ak(8,10,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2) + 5;
            a1 =  -ak(7,10,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k1 + 1) + 6;
            a1 = ak(8,11,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2 + 1) + 6;
            a1 =  -ak(7,11,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k1) + 6;
            a1 = ak(8,12,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2) + 6;
            a1 =  -ak(7,12,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            % -- s32( - )(beta,1)  =  s32( + ) (beta,n_gama);
            n = n + 1;
            nr1 = n;
            nc1 = nu(j + 1,k1) + 1;
            a1 = ak(10,1,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j + 1,k2) + 1;
            a1 =  -ak(9,1,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k1) + 1;
            a1 = ak(10,2,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2) + 1;
            a1 =  -ak(9,2,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j + 1,k1) + 2;
            a1 = ak(10,3,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j + 1,k2) + 2;
            a1 =  -ak(9,3,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k1) + 2;
            a1 = ak(10,4,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2) + 2;
            a1 =  -ak(9,4,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j + 1,k1) + 3;
            a1 = ak(10,5,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j + 1,k2) + 3;
            a1 =  -ak(9,5,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k1) + 3;
            a1 = ak(10,6,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);
 
            nr1 = n;
            nc1 = nu(j,k2) + 3;
            a1 =  -ak(9,6,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k1 + 1) + 4;
            a1 = ak(10,7,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2 + 1) + 4;
            a1 =  -ak(9,7,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k1) + 4;
            a1 = ak(10,8,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2) + 4;
            a1 =  -ak(9,8,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k1 + 1) + 5;
            a1 = ak(10,9,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2 + 1) + 5;
            a1 =  -ak(9,9,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k1) + 5;
            a1 = ak(10,10,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2) + 5;
            a1 =  -ak(9,10,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k1 + 1) + 6;
            a1 = ak(10,11,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2 + 1) + 6;
            a1 =  -ak(9,11,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k1) + 6;
            a1 = ak(10,12,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2) + 6;
            a1 =  -ak(9,12,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            % -- s33( - )(beta,1)  =  s33( + )(beta,n_gama);
            n = n + 1;
            nr1 = n;
            nc1 = nu(j + 1,k1) + 1;
            a1 = ak(12,1,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j + 1,k2) + 1;
            a1 =  -ak(11,1,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k1) + 1;
            a1 = ak(12,2,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2) + 1;
            a1 =  -ak(11,2,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j + 1,k1) + 2;
            a1 = ak(12,3,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j + 1,k2) + 2;
            a1 =  -ak(11,3,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k1) + 2;
            a1 = ak(12,4,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2) + 2;
            a1 =  -ak(11,4,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j + 1,k1) + 3;
            a1 = ak(12,5,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j + 1,k2) + 3;
            a1 =  -ak(11,5,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k1) + 3;
            a1 = ak(12,6,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2) + 3;
            a1 =  -ak(11,6,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k1 + 1) + 4;
            a1 = ak(12,7,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2 + 1) + 4;
            a1 =  -ak(11,7,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k1) + 4;
            a1 = ak(12,8,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2) + 4;
            a1 =  -ak(11,8,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k1 + 1) + 5;
            a1 = ak(12,9,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2 + 1) + 5;
            a1 =  -ak(11,9,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k1) + 5;
            a1 = ak(12,10,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2) + 5;
            a1 =  -ak(11,10,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k1 + 1) + 6;
            a1 = ak(12,11,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2 + 1) + 6;
            a1 =  -ak(11,11,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);
 
            nr1 = n;
            nc1 = nu(j,k1) + 6;
            a1 = ak(12,12,nsub1);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

            nr1 = n;
            nc1 = nu(j,k2) + 6;
            a1 =  -ak(11,12,nsub2);
            [NR, NC, A, ic] =  fn(nr1,nc1,a1,ic,NR,NC,A);

        end

    end % -- End long j loop
    
    % -- Define Omega in matlab sparse format
    Omega = sparse( NR(1:ic), NC(1:ic), A(1:ic) );

end
%^^^^^^^^^^^^^^^^^^^^^^^^ end of lhs ^^^^^^^^^^^^^^^^^^^^^^^

% -- Write out sparse matrix
% filename = 'Sparse_Omega.txt'; 
% fid = fopen(filename, 'wt');
% for I = 1: ic
%     text = [num2str(NR(I)), ', ', num2str(NC(I)),  ', ', num2str(A(I))];
%     fprintf(fid, '%s \n', text);
% end
% fclose(fid);

%^^^^^^^^^^^^^^^^^^^^^^^^ start rhs ^^^^^^^^^^^^^^^^^^^^^^^^
    
f = zeros(NEQ,1);
n = 0;

% -- s2j( + )(beta,gama)  =  s2j( - )(beta + 1,gama);
for  j = 1:NB - 1
    for  k = 1:NG
                   
        % -- C and gama for this subcell and next subcell
        C = props(j,k).C;
        gama = props(j,k).gama;
        CN = props(j + 1,k).C;
        gamaN = props(j + 1,k).gama;
        
        n = n + 1;
        f(n) = -2*(C(6,6) - CN(6,6))*epsb(6);

        n = n + 1;
        f(n) = (gama(2) - gamaN(2))*dtemp - (C(1,2) - CN(1,2))*epsb(1) ...
             - (C(2,2) - CN(2,2))*epsb(2) - (C(2,3) - CN(2,3))*epsb(3);

        n = n + 1;
        f(n) =  -2*(C(4,4) - CN(4,4))*epsb(4);
        
    end
end


% -- s3j( + )(beta,gama)  =  s3j( - )(beta + 1,gama);
for  k = 1:NG - 1
    for  j = 1:NB

        % -- C and gama for this subcell and next subcell
        C = props(j,k).C;
        gama = props(j,k).gama;
        CN = props(j,k + 1).C;
        gamaN = props(j,k + 1).gama;

        n = n + 1;
        f(n) = -2*(C(5,5) - CN(5,5))*epsb(5);

        n = n + 1;
        f(n) = -2*(C(4,4) - CN(4,4))*epsb(4);

        n = n + 1;

        f(n) = (gama(3) - gamaN(3))*dtemp - (C(1,3) - CN(1,3))*epsb(1) ...
             - (C(2,3) - CN(2,3))*epsb(2) - (C(3,3) - CN(3,3))*epsb(3);
         
    end
end

%********************************************!************************

%  -- Apply periodicity conditions at x2 = 0 & x2 = H
for  k = 1:NG

    % -- Pinned corners
    if ((pinned ~=0 && k == 1) || (pinned ~= 0) && (k == NG) )
        n = n + 6;
        
    else
        n = n + 3;

        % -- C and gama for first and last subcell
        C = props(1,k).C;
        gama = props(1,k).gama;
        CN = props(NB,k).C;
        gamaN = props(NB,k).gama;

        n = n + 1;
        f(n) =  -2*(C(6,6) - CN(6,6))*epsb(6);

        n = n + 1;
        f(n) = (gama(2) - gamaN(2))*dtemp - (C(1,2) - CN(1,2))*epsb(1) ...
             - (C(2,2) - CN(2,2))*epsb(2) - (C(2,3) - CN(2,3))*epsb(3);

        n = n + 1;
        f(n) =  -2*(C(4,4) - CN(4,4))*epsb(4);

    end
end

% -- Apply periodicity conditions at x3 = 0 & x3 = L
for  j = 1:NB

    % -- Pinned corners
    if ( (pinned ~= 0) && (j == 1) ) || (pinned ~= 0 && j == NB)
        n = n + 6;
        
    else
        n = n + 3;

        % -- s3j( - )(beta,1)  =  s3j( + )(beta,n_gama);
                   
        % -- C and gama for 1st and last subcell
        C = props(j,1).C;
        gama = props(j,1).gama;
        CN = props(j,NG).C;
        gamaN = props(j,NG).gama;                   

        n = n + 1;
        f(n) =  -2*(C(5,5) - CN(5,5))*epsb(5);

        n = n + 1;
        f(n) =  -2*(C(4,4) - CN(4,4))*epsb(4);

        n = n + 1;
        f(n) = (gama(3) - gamaN(3))*dtemp - (C(1,3) - CN(1,3))*epsb(1) ...
             - (C(2,3) - CN(2,3))*epsb(2) - (C(3,3) - CN(3,3))*epsb(3);

    end
end

% % -- Write out rhs vector
% filename = 'RHS.txt'; 
% fid = fopen(filename, 'wt');
% for I = 1: n
%     text = num2str(f(I));
%     fprintf(fid, '%s \n', text);
% end
% fclose(fid);
%^^^^^^^^^^^^^^^^^^^^^^^^ end of rhs ^^^^^^^^^^^^^^^^^^^^^^^

% -- Turn off MATLAB warning about Omega (clutters output)
warning('off','MATLAB:nearlySingularMatrix')

% -- Solve the system Omega*U = f (Eq. 6.61)
U = Omega\f;

end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ak, dd] = local(RUC, props)

% -- Determine the local stiffness matrix, K, and the terms Dij (Eqs. 6.41 - 6.47)

    NB = RUC.NB;
    NG = RUC.NG;
    nsubs = NB*NG;
    
    xh = RUC.h;
    xl = RUC.l;

    % -- Preallocate
    ak = zeros(12,12,nsubs);
    dd = zeros(3,3,nsubs);

    for j = 1:NB
        for k = 1:NG
            
            nsub = j + NB*(k - 1);
            
            C = props(j, k).C;

            coef1 = 4*(C(6,6)/xh(j)^2 + C(5,5)/xl(k)^2);
            d(1,2) = 2*C(6,6)/(xh(j)^2*coef1); % -- Eq. 6.41
            d(1,3) = 2*C(5,5)/(xl(k)^2*coef1); % -- Eq. 6.42

            coef2 = 4*(C(2,2)/xh(j)^2 + C(4,4)/xl(k)^2);
            d(2,2) = 2*C(2,2)/(xh(j)^2*coef2); % -- Eq. 6.43
            d(2,3) = 2*C(4,4)/(xl(k)^2*coef2); % -- Eq. 6.44

            coef3 = 4*(C(4,4)/xh(j)^2 + C(3,3)/xl(k)^2);
            d(3,2) = 2*C(4,4)/(xh(j)^2*coef3); % -- Eq. 6.45
            d(3,3) = 2*C(3,3)/(xl(k)^2*coef3); % -- Eq. 6.46

            for l = 1:3
                f(l,2) = 2/xh(j)^2 - 4*d(l,2)/xh(j)^2;
                f(l,3) =           - 4*d(l,3)/xh(j)^2;

                g(l,2) =           - 4*d(l,2)/xl(k)^2;
                g(l,3) = 2/xl(k)^2 - 4*d(l,3)/xl(k)^2;
            end

            %----------------------------------------------------------
            % -- Elements of the K matrix
            %----------------------------------------------------------

            ak(1,1,nsub) =  C(6,6)/xh(j) + C(6,6)*3*xh(j)*f(1,2)/2;
            ak(1,2,nsub) =  -C(6,6)/xh(j) + C(6,6)*3*xh(j)*f(1,2)/2;
            ak(1,7,nsub) =  C(6,6)*3*xh(j)*f(1,3)/2;
            ak(1,8,nsub) =  C(6,6)*3*xh(j)*f(1,3)/2;

            ak(2,1,nsub) =  C(6,6)/xh(j) - C(6,6)*3*xh(j)*f(1,2)/2;
            ak(2,2,nsub) =  -C(6,6)/xh(j) - C(6,6)*3*xh(j)*f(1,2)/2;
            ak(2,7,nsub) =  -C(6,6)*3*xh(j)*f(1,3)/2;
            ak(2,8,nsub) =  -C(6,6)*3*xh(j)*f(1,3)/2;

            ak(3,3 ,nsub) =  C(2,2)/xh(j) + C(2,2)*3*xh(j)*f(2,2)/2;
            ak(3,4 ,nsub) =  -C(2,2)/xh(j) + C(2,2)*3*xh(j)*f(2,2)/2;
            ak(3,9 ,nsub) =  C(2,2)*3*xh(j)*f(2,3)/2;
            ak(3,10,nsub) =  C(2,2)*3*xh(j)*f(2,3)/2;
            ak(3,11,nsub) =  C(2,3)/xl(k);
            ak(3,12,nsub) =  -C(2,3)/xl(k);

            ak(4,3 ,nsub) =  C(2,2)/xh(j) - C(2,2)*3*xh(j)*f(2,2)/2;
            ak(4,4 ,nsub) =  -C(2,2)/xh(j) - C(2,2)*3*xh(j)*f(2,2)/2;
            ak(4,9 ,nsub) =  -C(2,2)*3*xh(j)*f(2,3)/2;
            ak(4,10,nsub) =  -C(2,2)*3*xh(j)*f(2,3)/2;
            ak(4,11,nsub) =  C(2,3)/xl(k);
            ak(4,12,nsub) =  -C(2,3)/xl(k);

            ak(5,5 ,nsub) =  C(4,4)/xh(j) + C(4,4)*3*xh(j)*f(3,2)/2;
            ak(5,6 ,nsub) =  -C(4,4)/xh(j) + C(4,4)*3*xh(j)*f(3,2)/2;
            ak(5,9 ,nsub) =  C(4,4)/xl(k);
            ak(5,10,nsub) =  -C(4,4)/xl(k);
            ak(5,11,nsub) =  C(4,4)*3*xh(j)*f(3,3)/2;
            ak(5,12,nsub) =  C(4,4)*3*xh(j)*f(3,3)/2;

            ak(6,5 ,nsub) =  C(4,4)/xh(j) - C(4,4)*3*xh(j)*f(3,2)/2;
            ak(6,6 ,nsub) =  -C(4,4)/xh(j) - C(4,4)*3*xh(j)*f(3,2)/2;
            ak(6,9 ,nsub) =  C(4,4)/xl(k);
            ak(6,10,nsub) =  -C(4,4)/xl(k);
            ak(6,11,nsub) =  -C(4,4)*3*xh(j)*f(3,3)/2;
            ak(6,12,nsub) =  -C(4,4)*3*xh(j)*f(3,3)/2;

            ak(7,1 ,nsub) =  C(5,5)*3*xl(k)*g(1,2)/2;
            ak(7,2 ,nsub) =  C(5,5)*3*xl(k)*g(1,2)/2;
            ak(7,7 ,nsub) =  C(5,5)/xl(k) + C(5,5)*3*xl(k)*g(1,3)/2;
            ak(7,8 ,nsub) =  -C(5,5)/xl(k) + C(5,5)*3*xl(k)*g(1,3)/2;

            ak(8,1 ,nsub) =  - C(5,5)*3*xl(k)*g(1,2)/2;
            ak(8,2 ,nsub) =  - C(5,5)*3*xl(k)*g(1,2)/2;
            ak(8,7 ,nsub) =  C(5,5)/xl(k) - C(5,5)*3*xl(k)*g(1,3)/2;
            ak(8,8 ,nsub) =  -C(5,5)/xl(k) - C(5,5)*3*xl(k)*g(1,3)/2;

            ak(9,3 ,nsub) =  C(4,4)*3*xl(k)*g(2,2)/2;
            ak(9,4 ,nsub) =  C(4,4)*3*xl(k)*g(2,2)/2;
            ak(9,5 ,nsub) =  C(4,4)/xh(j);
            ak(9,6 ,nsub) =  -C(4,4)/xh(j);
            ak(9,9 ,nsub) =  C(4,4)/xl(k) + C(4,4)*3*xl(k)*g(2,3)/2;
            ak(9,10,nsub) =  -C(4,4)/xl(k) + C(4,4)*3*xl(k)*g(2,3)/2;

            ak(10,3 ,nsub) =  -C(4,4)*3*xl(k)*g(2,2)/2;
            ak(10,4 ,nsub) =  -C(4,4)*3*xl(k)*g(2,2)/2;
            ak(10,5 ,nsub) =  C(4,4)/xh(j);
            ak(10,6 ,nsub) =  -C(4,4)/xh(j);
            ak(10,9 ,nsub) =  C(4,4)/xl(k) - C(4,4)*3*xl(k)*g(2,3)/2;
            ak(10,10,nsub) =  - C(4,4)/xl(k) - C(4,4)*3*xl(k)*g(2,3)/2;

            ak(11,3 ,nsub) =  C(2,3)/xh(j);
            ak(11,4 ,nsub) =  -C(2,3)/xh(j);
            ak(11,5 ,nsub) =  C(3,3)*3*xl(k)*g(3,2)/2;
            ak(11,6 ,nsub) =  C(3,3)*3*xl(k)*g(3,2)/2;
            ak(11,11,nsub) =  C(3,3)/xl(k) + C(3,3)*3*xl(k)*g(3,3)/2;
            ak(11,12,nsub) =  -C(3,3)/xl(k) + C(3,3)*3*xl(k)*g(3,3)/2;

            ak(12,3 ,nsub) =  C(2,3)/xh(j);
            ak(12,4 ,nsub) =  -C(2,3)/xh(j);
            ak(12,5 ,nsub) =  -C(3,3)*3*xl(k)*g(3,2)/2;
            ak(12,6 ,nsub) =  -C(3,3)*3*xl(k)*g(3,2)/2;
            ak(12,11,nsub) =  C(3,3)/xl(k) - C(3,3)*3*xl(k)*g(3,3)/2;
            ak(12,12,nsub) =  -C(3,3)/xl(k) - C(3,3)*3*xl(k)*g(3,3)/2;

            % -- Store d per subcell for return
            for  l = 1:3
                for  m = 1:3
                    dd(l, m, nsub) = d(l, m);
                end
            end
        end
        
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [w] = micro_variables(j,k,NB,RUC,dd,x)

% -- Determine the microvariables, w, for a given subcell (beta, gama)

    nsub = j + NB*(k - 1);

    n  = nu(j    ,k);
    n2 = nu(j + 1,k);
    n3 = nu(j    ,k + 1);
    
    % -- Eq. 6.40
    w{1}.w00 = dd(1,2,nsub)*(x(n2 + 1) + x(n + 1)) + dd(1,3,nsub)*(x(n3 + 4) + x(n + 4));
    w{2}.w00 = dd(2,2,nsub)*(x(n2 + 2) + x(n + 2)) + dd(2,3,nsub)*(x(n3 + 5) + x(n + 5));
    w{3}.w00 = dd(3,2,nsub)*(x(n2 + 3) + x(n + 3)) + dd(3,3,nsub)*(x(n3 + 6) + x(n + 6));

    % -- Eq. 6.33
    w{1}.w10 = (x(n2 + 1) - x(n + 1))/RUC.h(j);
    w{2}.w10 = (x(n2 + 2) - x(n + 2))/RUC.h(j);
    w{3}.w10 = (x(n2 + 3) - x(n + 3))/RUC.h(j);

    % -- Eq. 6.34
    w{1}.w01 = (x(n3 + 4) - x(n + 4))/RUC.l(k);
    w{2}.w01 = (x(n3 + 5) - x(n + 5))/RUC.l(k);
    w{3}.w01 = (x(n3 + 6) - x(n + 6))/RUC.l(k);

    % -- Eq. 6.35
    w{1}.w20 = (2*(x(n2 + 1) + x(n + 1)) - 4*w{1}.w00)/RUC.h(j)^2;
    w{2}.w20 = (2*(x(n2 + 2) + x(n + 2)) - 4*w{2}.w00)/RUC.h(j)^2;
    w{3}.w20 = (2*(x(n2 + 3) + x(n + 3)) - 4*w{3}.w00)/RUC.h(j)^2;

    % -- Eq. 6.36
    w{1}.w02 = (2*(x(n3 + 4) + x(n + 4)) - 4*w{1}.w00)/RUC.l(k)^2;
    w{2}.w02 = (2*(x(n3 + 5) + x(n + 5)) - 4*w{2}.w00)/RUC.l(k)^2;
    w{3}.w02 = (2*(x(n3 + 6) + x(n + 6)) - 4*w{3}.w00)/RUC.l(k)^2;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [NR, NC, A, ic] = fn(nr1, nc1, a1, ic, NR, NC, A)

% -- Place value, row, and column in sparsely-stored matrix A

    tiny = 1.d-20;

    % -- Ignore very small numbers
    if abs(a1) < tiny
        return
    end

    ic = ic + 1;
    NR(ic) = nr1;
    NC(ic) = nc1;
    A(ic) = a1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nuresult = nu(j, k)

    % -- Utility function
    % -- For a given(j,k), the unknowns are at nuresult + 1,2,...;

    global NB NG

    if (j <= NB && k <= NG)
        nsub = j + NB*(k - 1);
        nuresult = 6*(nsub - 1);

    elseif (j == NB + 1)
        nuresult = 6*NB*NG + 3*(k-1);

    elseif(k == NG + 1)
        nuresult = 6*NB*NG + 3*NG + 3*(j - 1) - 3;

    end

end



