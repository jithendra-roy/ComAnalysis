%Nominal moment of section calculator
clc;
clear all;
filename = 'SectionInfo.xlsx';
sheetname = 'Data';
cell = 'A121:G121';         % change the number accordingly---(S.No+1)
O = xlsread(filename,sheetname,cell);
b        = O(1,1);      % b = width of section
d        = O(1,2);      % d = effective depth
dc       = O(1,3);      % dc = cover to compressive steel
Ast      = O(1,4);      % Ast = Area of tension reinforcement
Asc      = O(1,5);      % Asc = Area of compression reinforcement
fy       = O(1,6);      % fy = characteristic strength of steel
fck      = O(1,7);      % sigmancbc = characteristic compressive stress of concrete
if fy == 250
    fsc = 0.0035*(0.53*d-dc)/(0.53*d);
    xu = (0.87*fy*Ast-fsc*Asc)/(0.36*fck*b);
    xumax = 0.53*d;
    if xu<xumax        %Under reinforced
        fprintf("Under reinforced\n")
        M = 0.36*(xu/d)*(1-0.42*(xu/d))*b*d*d*fck+fsc*Asc*(d-dc);
    else
        fprintf("over reinforced\n")
        M = 0.36*0.53*(1-0.42*0.53)*b*d*d*fck+fsc*Asc*(d-dc);
    end
elseif fy == 415
    fsc = 0.0035*(0.48*d-dc)/(0.48*d);
    xu = (0.87*fy*Ast-fsc*Asc)/(0.36*fck*b);
    xumax = 0.48*d;
    if xu<xumax        %Under reinforced
        fprintf("Under reinforced\n")
        M = 0.36*(xu/d)*(1-0.42*(xu/d))*b*d*d*fck+fsc*Asc*(d-dc);
    else
        fprintf("over reinforced\n")
        M = 0.36*0.48*(1-0.42*0.48)*b*d*d*fck+fsc*Asc*(d-dc);
    end
elseif fy == 500
    %fsc = 0.0035*(0.46*d-dc)/(0.46*d);
    fsc = 0.85*fy;
    xu = (0.87*fy*Ast-fsc*Asc)/(0.36*fck*b);
    xumax = 0.46*d;
    if xu<xumax        %Under reinforced
        fprintf("Under reinforced\n")
        M = 0.36*(xu/d)*(1-0.42*(xu/d))*b*d*d*fck+fsc*Asc*(d-dc);
    else
        fprintf("over reinforced\n")
        M = 0.36*0.46*(1-0.42*0.46)*b*d*d*fck+fsc*Asc*(d-dc);
    end
end
fprintf("b = %f\n",b)
fprintf("d = %f\n",d)
fprintf("fsc = %f\n",fsc)
fprintf("xu = %f\n",xu)
fprintf("xumax = %f\n",xumax)
fprintf("M = %f\n",M/10^7)