clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%坐标
Coords=load('Coords.txt');
size_node=size(Coords,1);
Elems=load('Elems.txt');
size_element=size(Elems,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tt=1;
E0=22000;
mu=1/3;
DD=[1 mu 0;mu 1 0;0 0 (1-mu)/2];
DD=DD*E0/(1-mu*mu);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 边界节点
ii_left=0;
boundary_left=[];
ii_right=0;
boundary_right=[];
ii_bottom=0;
for i=1:size_node
    xx=Coords(i,1);
    if xx<0.5
        ii_left=ii_left+1;
        boundary_left(ii_left,1)=i;
    end
    if xx>199.5
        ii_right=ii_right+1;
        boundary_right(ii_right,1)=i;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PD粒子
PD_lizi=[];
ii=0;
for i=1:size_node
    xx=Coords(i,1);
    yy=Coords(i,2);
    if xx>49.5 && xx<150.5
        if yy>29.5 && yy<70.5
            ii=ii+1;
            PD_lizi(ii,1)=i;
        end
    end
end
size_PD=size(PD_lizi,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 域内粒子
dettx=1.0;
detla=3.015*dettx;
vx=dettx*dettx*tt;
c=6*E0/(pi*tt*detla*detla*detla);
horizon_parrticle=zeros(size_PD,28);
for i=1:size_PD
    ii=PD_lizi(i,1);
    xxi=Coords(ii,1);
    yyi=Coords(ii,2);
    jj=0;
    for j=1:size_node
        if ii~=j
            xxj=Coords(j,1);
            yyj=Coords(j,2); 
            ll=sqrt((xxi-xxj)*(xxi-xxj)+(yyi-yyj)*(yyi-yyj));
            if ll<detla
                jj=jj+1;
                horizon_parrticle(i,jj)=j;
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FEM求解
kz1=zeros(size_node*2,size_node*2);
pz1=zeros(size_node*2,1);
uz1=zeros(size_node*2,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%积分点位置及权函数
pointr(1)=-1.0/sqrt(3);
pointr(2)=1.0/sqrt(3);
weightr(1)=1.0;
weightr(2)=1.0;
points(1)=-1.0/sqrt(3);
points(2)=1.0/sqrt(3);
weights(1)=1.0;
weights(2)=1.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:size_element
    i1=Elems(i,1);
    i2=Elems(i,2);
    i3=Elems(i,3);
    i4=Elems(i,4);
    xxi=Coords(i1,1);
    yyi=Coords(i1,2);
    xxj=Coords(i2,1);
    yyj=Coords(i2,2);
    xxm=Coords(i3,1);
    yym=Coords(i3,2);
    xxn=Coords(i4,1);
    yyn=Coords(i4,2);
    ke=zeros(8);
	enode=[xxi yyi; xxj yyj; xxm yym; xxn yyn];               
    for ii=1:2
        for jj=1:2
            rr=pointr(ii);
            ss=points(jj);
            dN(1,1)=-(1-ss)/4;
            dN(1,2)=(1-ss)/4;
            dN(1,3)=(1+ss)/4;
            dN(1,4)=-(1+ss)/4;
            dN(2,1)=-(1-rr)/4;
            dN(2,2)=-(1+rr)/4;
            dN(2,3)=(1+rr)/4;
            dN(2,4)=(1-rr)/4;
            J=dN*enode;
            detJ=det(J);
            dNxyz=inv(J)*dN;
            dNdx=dNxyz(1,1);
            dNdy=dNxyz(2,1);
            B1=[dNdx 0;0 dNdy;dNdy dNdx];
            dNdx=dNxyz(1,2);
            dNdy=dNxyz(2,2);
            B2=[dNdx 0;0 dNdy;dNdy dNdx];
            dNdx=dNxyz(1,3);
            dNdy=dNxyz(2,3);
            B3=[dNdx 0;0 dNdy;dNdy dNdx];
            dNdx=dNxyz(1,4);
            dNdy=dNxyz(2,4);
            B4=[dNdx 0;0 dNdy;dNdy dNdx];
            BB=[B1 B2 B3 B4];
            ke=ke+BB'*DD*BB*tt*detJ*weightr(ii)*weights(jj);
        end
    end
    for ij=1:4
        kz1(Elems(i,ij)*2-1,Elems(i,1)*2-1)=kz1(Elems(i,ij)*2-1,Elems(i,1)*2-1)+ke(ij*2-1,1);
        kz1(Elems(i,ij)*2-1,Elems(i,1)*2-0)=kz1(Elems(i,ij)*2-1,Elems(i,1)*2-0)+ke(ij*2-1,2);
        kz1(Elems(i,ij)*2-0,Elems(i,1)*2-1)=kz1(Elems(i,ij)*2-0,Elems(i,1)*2-1)+ke(ij*2-0,1);
        kz1(Elems(i,ij)*2-0,Elems(i,1)*2-0)=kz1(Elems(i,ij)*2-0,Elems(i,1)*2-0)+ke(ij*2-0,2);                                     
        kz1(Elems(i,ij)*2-1,Elems(i,2)*2-1)=kz1(Elems(i,ij)*2-1,Elems(i,2)*2-1)+ke(ij*2-1,3);
        kz1(Elems(i,ij)*2-1,Elems(i,2)*2-0)=kz1(Elems(i,ij)*2-1,Elems(i,2)*2-0)+ke(ij*2-1,4);
        kz1(Elems(i,ij)*2-0,Elems(i,2)*2-1)=kz1(Elems(i,ij)*2-0,Elems(i,2)*2-1)+ke(ij*2-0,3);
        kz1(Elems(i,ij)*2-0,Elems(i,2)*2-0)=kz1(Elems(i,ij)*2-0,Elems(i,2)*2-0)+ke(ij*2-0,4);
        kz1(Elems(i,ij)*2-1,Elems(i,3)*2-1)=kz1(Elems(i,ij)*2-1,Elems(i,3)*2-1)+ke(ij*2-1,5);
        kz1(Elems(i,ij)*2-1,Elems(i,3)*2-0)=kz1(Elems(i,ij)*2-1,Elems(i,3)*2-0)+ke(ij*2-1,6);
        kz1(Elems(i,ij)*2-0,Elems(i,3)*2-1)=kz1(Elems(i,ij)*2-0,Elems(i,3)*2-1)+ke(ij*2-0,5);
        kz1(Elems(i,ij)*2-0,Elems(i,3)*2-0)=kz1(Elems(i,ij)*2-0,Elems(i,3)*2-0)+ke(ij*2-0,6);
        kz1(Elems(i,ij)*2-1,Elems(i,4)*2-1)=kz1(Elems(i,ij)*2-1,Elems(i,4)*2-1)+ke(ij*2-1,7);
        kz1(Elems(i,ij)*2-1,Elems(i,4)*2-0)=kz1(Elems(i,ij)*2-1,Elems(i,4)*2-0)+ke(ij*2-1,8);
        kz1(Elems(i,ij)*2-0,Elems(i,4)*2-1)=kz1(Elems(i,ij)*2-0,Elems(i,4)*2-1)+ke(ij*2-0,7);
        kz1(Elems(i,ij)*2-0,Elems(i,4)*2-0)=kz1(Elems(i,ij)*2-0,Elems(i,4)*2-0)+ke(ij*2-0,8);
    end   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PD总体刚度矩阵
PD_lizi1=[PD_lizi*2-1;PD_lizi*2];
kz1([PD_lizi1],:)=[0];
for i=1:size_PD
    ii=PD_lizi(i,1);
    xxi=Coords(ii,1);
    yyi=Coords(ii,2);
    for j=1:28
        jj=horizon_parrticle(i,j);
        if jj>0
            xxj=Coords(jj,1);
            yyj=Coords(jj,2);
            l=sqrt((xxi-xxj)*(xxi-xxj)+(yyi-yyj)*(yyi-yyj));
            dex=(xxi-xxj);
            dey=(yyi-yyj);
            kk_PD=zeros(2,4);
            kk_PD(1,1)=vx*vx*c*dex*dex/l;
            kk_PD(1,2)=vx*vx*c*dex*dey/l;
            kk_PD(1,3)=-vx*vx*c*dex*dex/l;
            kk_PD(1,4)=-vx*vx*c*dex*dey/l;
            kk_PD(2,1)=vx*vx*c*dex*dey/l;
            kk_PD(2,2)=vx*vx*c*dey*dey/l;
            kk_PD(2,3)=-vx*vx*c*dex*dey/l;
            kk_PD(2,4)=-vx*vx*c*dey*dey/l;
            kz1(ii*2-1,ii*2-1)=kz1(ii*2-1,ii*2-1)+kk_PD(1,1);
            kz1(ii*2-1,ii*2)=kz1(ii*2-1,ii*2)+kk_PD(1,2);
            kz1(ii*2-1,jj*2-1)=kz1(ii*2-1,jj*2-1)+kk_PD(1,3);
            kz1(ii*2-1,jj*2)=kz1(ii*2-1,jj*2)+kk_PD(1,4);
            kz1(ii*2,ii*2-1)=kz1(ii*2,ii*2-1)+kk_PD(2,1);
            kz1(ii*2,ii*2)=kz1(ii*2,ii*2)+kk_PD(2,2);
            kz1(ii*2,jj*2-1)=kz1(ii*2,jj*2-1)+kk_PD(2,3);
            kz1(ii*2,jj*2)=kz1(ii*2,jj*2)+kk_PD(2,4);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 施加边界条件并求解
for i=1:ii_left
    ii=boundary_left(i,1)*2-1;
    kz1(ii,ii)=kz1(ii,ii)*10e8;
    pz1(ii,1)=kz1(ii,ii)*0;
    ii=boundary_left(i,1)*2;
    kz1(ii,ii)=kz1(ii,ii)*10e8;
    pz1(ii,1)=kz1(ii,ii)*0;
end
for i=1:ii_right
    ii=boundary_right(i,1)*2-1;
    kz1(ii,ii)=kz1(ii,ii)*10e8;
    pz1(ii,1)=kz1(ii,ii)*0.1;
    ii=boundary_right(i,1)*2;
    kz1(ii,ii)=kz1(ii,ii)*10e8;
    pz1(ii,1)=kz1(ii,ii)*0;
end
uz1=kz1\pz1;
ii=1;
for i=1:101
    for j=1:201
        ux(i,j)=uz1(ii*2-1);
        ii=ii+1;
    end
end
ii=1;
for i=1:101
    for j=1:201
        uy(i,j)=uz1(ii*2);
        ii=ii+1;
    end
end



