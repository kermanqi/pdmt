clear
time=1;

for i=1:5001
	node(i,1)=(i-1);
end
a1=size(node,1);
node=node*0.02;
for i=1:5000
    element(i,1)=i;
    element(i,2)=i+1;
end
a22=size(element,1);

ii1=0;
ii2=0;
for i=1:5001
    xx=node(i,1);
    if xx>30 && xx<70
        ii1=ii1+1;
        node_PD(ii1)=i;
    else
        ii2=ii2+1;
        node_FEM(ii2)=i;
    end
end
a1_PD=size(node_PD,2);
a1_FEM=size(node_FEM,2);

tic
kkk=1;
E0=210000;
alf=1.0e-6;
mmm=3;
yuneilizi=6;
detx=0.02;
detla=(mmm+0.01)*detx;
fai=1;
A0=0.0004;
vx=detx*A0;
c0=2/6/vx;
%c0=1/A0/detla;
yunei=zeros(a1,yuneilizi);
for i=1:a1
    xx1=node(i,1);
    j1=0;
    for j=1:a1
            xx2=node(j,1);
            ll=abs(xx2-xx1);
            if ll<detla && ll>0
                j1=j1+1;
                yunei(i,j1)=j;
            end
    end
end


c00=zeros(a1,1);
for i=1:a1
    for j=1:yuneilizi
        if yunei(i,j)>0
            c00(i,1)=c0*yuneilizi/j;
          %  c00(i,1)=c0;
        end
    end
end



kz1_t=zeros(a1);
fz1=zeros(a1,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FEM温度场
for i=1:a22
    i1=element(i,1);
    i2=element(i,2);
    xxi=node(i1,1);
    xxj=node(i2,1);
    le=xxj-xxi;
    ke_t=kkk*A0/le*[1 -1;-1 1];
    line=[i1 i2];
    kz1_t(line,line)=kz1_t(line,line)+ke_t;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PD温度场
kz1_t(node_PD,:)=[0];
for ii=1:a1_PD
	i=node_PD(ii);
    xxi=node(i,1);
    for j=1:6
        jj=yunei(i,j);
        if jj>0
            xxj=node(jj,1);  
            le=abs(xxj-xxi);
            ke_t=kkk*c0/le/le*[1 -1;-1 1]*vx*vx;
            line=[i jj];
            kz1_t(i,line)=kz1_t(i,line)+ke_t(1,:);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kz1_t(1,1)=kz1_t(1,1)*1.0e8;
kz1_t(5001,5001)=kz1_t(5001,5001)*1.0e8;
fz1(5001,1)=kz1_t(5001,5001)*100;
faiz=kz1_t\fz1;

kz1=zeros(a1);
pz1=zeros(a1,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FEM位移场
for i=1:a22
    i1=element(i,1);
    i2=element(i,2);
    xxi=node(i1,1);
    xxj=node(i2,1);
    le=xxj-xxi;
    ke=E0*A0/le*[1 -1;-1 1];
    line=[i1 i2];
    kz1(line,line)=kz1(line,line)+ke;
    fai_aver=(faiz(i1,1)+faiz(i2,1))/2;
    pz1(line,1)=pz1(line,1)+[-1/le;1/le]*E0*alf*fai_aver*A0*le;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PD位移场
kz1(node_PD,:)=[0];
pz1(node_PD,1)=[0];
for ii=1:a1_PD
	i=node_PD(ii);
    xxi=node(i,1);
    for j=1:6
        jj=yunei(i,j);
        if jj>0
            xxj=node(jj,1);  
            le=abs(xxj-xxi);
            ke=E0*c0/le/le*[1 -1;-1 1]*vx*vx;
            line=[i jj];
            kz1(i,line)=kz1(i,line)+ke(1,:);
            fai_aver=(faiz(i,1)+faiz(jj,1))/2;
            pz1(i,1)=pz1(i,1)+c0*[-1/(xxj-xxi)]*E0*alf*fai_aver*vx*vx;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kz1(1,1)=kz1(1,1)*1.0e8;
uz=kz1\pz1;
