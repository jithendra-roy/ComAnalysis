clc;
N=input('number of nodes in structure= ');
% number of nodes=17*5
%for portal no of nodes=4
M=input('number of members in frame= ');
% total members=335
%for portal members=3
a=input('total number of restrained degree of freedoms= ');
% foundation level=17*6
%for portal a=12
b=(6*N)-a;
% Input data file test in sheet x will contain joint x coordinates
X=xlsread('test.xlsx','x','A1:A102'); % X=xlsread('F Port.xlsx','coordinates','B2;B5')
% Input data file test in sheet y will contain joint y coordinates
Y=xlsread('test.xlsx','y','A1:A102'); % Y=xlsread('F Port.xlsx','coordinates','C2:C5')
% Input data file test in sheet z will contain joint Z coordinates
Z=xlsread('test.xlsx','z','A1:A102'); % Z=xlsread('F Port.xlsx','coordinates','D2:D5')
% Input data file test in sheet sectionproperties will contain Member properties
O=xlsread('test.xlsx','sectionproperties','A2:P336'); % O=xlsread('F Port.xlsx','sectionproperties','A2:P4')
Smi=zeros(12,12);
%Smi=Element stiffness matrix for member i IN MEMBER AXIS
Sms=zeros(12,12);
% Member stiffnes matrix for structure axis
Ss=zeros(6*N,6*N);
%Ss=system stiffness matrix
xj=zeros(12,1);
%xj is the displacement and rotations of each member for calculating member force
R=xlsread('test.xlsx','resDofDis','A1:B102'); % R=xlsread('F Port.xlsx','resDofDis','A1:B12');
%input data file test in sheet resDofDis will contain restrained dof`s
R1=xlsread('test.xlsx','freeDofDis','A1:B408'); %R1=xlsread('F Port.xlsx','freeDofDis','A1:A12');
%input data file test in sheet freeDofDis will contain free dof`s
Eq=xlsread('test.xlsx','freeDofLoad','A1:B408'); %Eq=xlsread('F Port.xlsx','freeDofload','A1:B12');
%input data file test in sheet freeDofLoad will contain equivalent load on
%dof`s
Rn=xlsread('test.xlsx','resDofLoad','A1:B102'); %Rn=xlsread('F Port','resDofload','A1:B12');
%input data file test in sheet resDofLoad will contain restained dof and
%corresponding eq loads

y=a+1;
Dk=R(1:a,2);
%known deformation of joints corrosponding to external reaction terms%
Ak=Eq(1:b,2);
% system loads corrosponding to joint degrees of freedom
B=R1(1:b,1);
% Position of joint degrees of freedoms
%formation of system stiffness matrix
for i=1:M %member number
    j=O(i,2); %node 1
    k=O(i,3); %node 2
    w1=O(i,10); % theeta 1
    w2=O(i,11); %theeta 2
    L=O(i,12); % length of member
    AE=O(i,13); % product of A and E
    EIz=O(i,14); % product of E and Iz
    EIy=O(i,15); % product of E and Iy
    GIx=O(i,16); % product of G and J
    C1=(O(i,7)-O(i,4))/O(i,12); % lambda X
    C2=(O(i,8)-O(i,5))/O(i,12); % lambda Y
    C3=(O(i,9)-O(i,6))/O(i,12); % lambda Z
    if (O(i,5) == O(i,8) && O(i,6) == O(i,9) && (O(i,7)-O(i,4))/O(i,12)==1) % member with lambda X positive
        Di=[1 0 0;0 0 1;0 -1 0];
    elseif(O(i,5) == O(i,8) && O(i,6) == O(i,9) && (O(i,7)-O(i,4))/O(i,12)==-1) %member with lambda X negative
        Di=[-1 0 0;0 0 1;0 1 0];
    elseif (O(i,4)== O(i,7) && O(i,6)== O(i,9) && (O(i,8)-O(i,5))/O(i,12)==1)% member with lambda Y positive
       Di = [0 1 0;0 0 1;1 0 0];
    elseif (O(i,4)== O(i,7) && O(i,6)== O(i,9) && (O(i,8)-O(i,5))/O(i,12)==-1)% member with lambda Y negative
       Di = [0 -1 0;0 0 1;-1 0 0];
    elseif (O(i,4)== O(i,7) && O(i,5)== O(i,8) && (O(i,9)-O(i,6))/O(i,12)==1) % member with lambda Z positive
         Di = [0 0 1;-1 0 0;0 -1 0];
    elseif (O(i,4)== O(i,7) && O(i,5)== O(i,8) && (O(i,9)-O(i,6))/O(i,12)==-1) % member with lambda Z negative
       Di = [0 0 -1;-1 0 0;0 1 0];
    elseif((O(i,7)-O(i,4))/O(i,12)==0 && (O(i,8)-O(i,5))/O(i,12)>0 ) % for yz plane and lambda y >0
        Di =[C1 C2 C3;-1 0 0;0 -C3 C2];
    elseif((O(i,7)-O(i,4))/O(i,12)==0 && (O(i,8)-O(i,5))/O(i,12)<0 ) % for yz plane and lambda y<0
        Di=[C1 C2 C3;1 0 0;0 C3 -C2];
    elseif((O(i,7)-O(i,4))/O(i,12)>0 && (O(i,8)-O(i,5))/O(i,12)==0) % for xz plane and lambda x>0
        Di=[C1 C2 C3;0 1 0;-C3 0 C1];
    elseif((O(i,7)-O(i,4))/O(i,12)<0 && (O(i,8)-O(i,5))/O(i,12)==0) % for xz plane and lambda x<0
        Di=[C1 C2 C3;0 -1 0;C3 0 -C1];
    end
    Dt=[Di Di-Di Di-Di Di-Di; Di-Di Di Di-Di Di-Di; Di-Di Di-Di Di Di-Di;Di-Di Di-Di Di-Di Di]; %Transformation matrix
    
    %Smi is member stiffness matrix with local coordinates
    Smi=[AE/L 0 0 0 0 0 -AE/L 0 0 0 0 0
        0 12*EIz/(L^3) 0 0 0 6*EIz/(L^2) 0 -12*EIz/(L^3) 0 0 0 6*EIz/(L^2)
        0 0 12*EIy/(L^3) 0 -6*EIy/(L^2) 0 0 0 -12*EIy/(L^3) 0 -6*EIy/(L^2) 0
        0 0 0 GIx/L 0 0 0 0 0 -GIx/L 0 0
        0 0 -6*EIy/(L^2) 0 4*EIy/L 0 0 0 6*EIy/(L^2) 0 2*EIy/L 0
        0 6*EIz/(L^2) 0 0 0 4*EIz/L 0 -6*EIz/(L^2) 0 0 0 2*EIz/L
        -AE/L 0 0 0 0 0 AE/L 0 0 0 0 0
        0 -12*EIz/(L^3) 0 0 0 -6*EIz/(L^2) 0 12*EIz/(L^3) 0 0 0 -6*EIz/(L^2)
        0 0 -12*EIy/(L^3) 0 6*EIy/(L^2) 0 0 0 12*EIy/(L^3) 0 6*EIy/(L^2) 0
        0 0 0 -GIx/L 0 0 0 0 0 GIx/L 0 0
        0 0 -6*EIy/(L^2) 0 2*EIy/L 0 0 0 6*EIy/(L^2) 0 4*EIy/L 0
        0 6*EIz/(L^2) 0 0 0 2*EIz/L 0 -6*EIz/(L^2) 0 0 0 4*EIz/L ];    %formation of element stiffness matrix IN MEMBER AXIS
    
    Smsi= Dt.'*Smi*Dt; 
    Smsi1=Smsi(1:6,1:6);
    Smsi2=Smsi(1:6,7:12);
    Smsi3=Smsi(7:12,1:6);
    Smsi4=Smsi(7:12,7:12);
    Ss(((6*j)-5):(6*j),((6*j)-5):(6*j))=(Smsi1)+Ss(((6*j)-5):(6*j),((6*j)-5):(6*j));
    Ss(((6*j)-5):(6*j),((6*k)-5):(6*k))=(Smsi2)+Ss(((6*j)-5):(6*j),((6*k)-5):(6*k));
    Ss(((6*k)-5):(6*k),((6*j)-5):(6*j))=(Smsi3)+Ss(((6*k)-5):(6*k),((6*j)-5):(6*j));
    Ss(((6*k)-5):(6*k),((6*k)-5):(6*k))=(Smsi4)+Ss(((6*k)-5):(6*k),((6*k)-5):(6*k)); 
end
%Rearrangement of stiffness matrix according to restrained degrees of
%freedom positon
for c=1:a
    d=R(c,1);   %dof number for restrained terms
    for C=1:a
        D=R(C,1);
        Srr(c,C)=Ss(d,D);
    end
    for C1=1:b
        D1=R1(C1,1);
        Srk(c,C1)=Ss(d,D1);
    end
end
for e=1:b
    f=R1(e,1);   %dof number for free joint terms
    for E=1:b
        F=R1(E,1);
        Skk(e,E)=Ss(f,F);
    end
    for E1=1:a
        F1=R(E1,1);
        Skr(e,E1)=Ss(f,F1);
    end
end
Dj=(inv(Skk))*(Ak-(Skr*Dk));  % joint deformations corrosponding to known system joint loads
Ar=-(Rn(1:a,2))+(Srk*Dj)+(Srr*Dk);         % system reaction terms
As=zeros((6*N),1);            % full system joint load matrix
Ds=zeros((6*N),1);            % full system joint deformation matrix
for g=1:a
    As(R(g,1),1)=Ar(g,1);
    Ds(R(g,1),1)=Dk(g,1);
end
for h=1:(6*N)-a
    As(R1(h,1),1)=Ak(h,1);
    Ds(R1(h,1),1)=Dj(h,1);
end
Ds;           %Displacement corrosponding each joint dof
Ar;           %Forces at all reaction members
fprintf(fopen('Ds.txt','w'), '%f\n', Ds);
fprintf(fopen('Ar.txt','w'), '%f\n', Ar);
