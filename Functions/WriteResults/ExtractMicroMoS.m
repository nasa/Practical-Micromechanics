function [DataCell] = ExtractMicroMoS(CN, Type, MoS, props, DataCell, LayerZ)
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
% Purpose: Setup cell for writing micromechanics-based MoS to output file
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if nargin == 5
    LayerZ = 0;
end

row = 1;

if Type == "FM"

    row = row + 1;
    Name = 'F';
    DataCell{row, 1} = num2str(LayerZ);
    DataCell{row, 2} = Name;

    if CN < 3
      for i = 1:6
        DataCell{row, i+2} = num2str(MoS.MicroF.Crit{CN}.MoS(i));
      end
    else
      DataCell{row, 3} = num2str(MoS.MicroF.Crit{CN}.MoS(1));
    end

    row = row + 1;
    Name = 'M';
    DataCell{row, 1} = ' ';
    DataCell{row, 2} = Name;

    if CN < 3
       for i = 1:6
         DataCell{row, i+2} = num2str(MoS.MicroM.Crit{CN}.MoS(i));
       end
    else
       DataCell{row, 3} = num2str(MoS.MicroM.Crit{CN}.MoS(1));
    end

else
  for b = 1:props.RUC.NB
     for g = 1:props.RUC.NG
       row = row + 1;
       Subcell = [num2str(b),',',num2str(g),' - ',props.RUC.matsCh(b,g)];
       if (b == 1) && (g == 1)
           DataCell{row, 1} = num2str(LayerZ);
       else
           DataCell{row, 1} = ' ';
       end
       DataCell{row, 2} = Subcell;
       if CN < 3
          for i = 1:6
            DataCell{row, i+2} = num2str(MoS.Micro(b,g).Crit{CN}.MoS(i));
          end
       else
          DataCell{row, 3} = num2str(MoS.Micro(b,g).Crit{CN}.MoS(1));
       end
     end
   end
end

end