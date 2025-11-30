clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 节点信息
ii=1;
for i=1:101
    for j=1:101
        coords(ii,2)=i-51;
        coords(ii,1)=j-51;
        ii=ii+1;
    end
end
size_node=size(coords,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 单元信息
ii=1;
i1=1;
j1=102;
for i=1:100
    for j=1:100
        elems(ii,1)=i1;
        elems(ii,2)=i1+1;
        elems(ii,3)=j1+1;
        elems(ii,4)=j1;
        ii=ii+1;
        i1=i1+1;
        j1=j1+1;
    end
        i1=i1+1;
        j1=j1+1;
end
size_element=size(elems,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 边界节点
ii_left=0;
boundary_left=[];
ii_right=0;
boundary_right=[];
ii_bottom=0;
boundary_bottom=[];
ii_top=0;
boundary_top=[];
for i=1:size_node
    xx=coords(i,1);
    yy=coords(i,2);
    if xx<-49.5
        ii_left=ii_left+1;
        boundary_left(ii_left,1)=i;
    end
    if xx>49.5
        ii_right=ii_right+1;
        boundary_right(ii_right,1)=i;
    end
    if yy<-49.5 && xx>-49.5 && xx<49.5
        ii_bottom=ii_bottom+1;
        boundary_bottom(ii_bottom,1)=i;
    end
    if yy>49.5 && xx>-49.5 && xx<49.5
        ii_top=ii_top+1;
        boundary_top(ii_top,1)=i;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PD粒子
PD_lizi=[];
ii=0;
for i=1:size_node
    xx=coords(i,1);
    yy=coords(i,2);
    if abs(xx)<30.5 && abs(yy)<30.5
        ii=ii+1;
        PD_lizi(ii,1)=i;
    end
end
size_PD=size(PD_lizi,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 域内粒子
k0=1;   
detla=3.015;
tt=1;
c=6*k0/(pi*tt*detla*detla*detla);
horizon_parrticle=zeros(size_PD,28);
for i=1:size_PD
    ii=PD_lizi(i,1);
    xxi=coords(ii,1);
    yyi=coords(ii,2);
    jj=0;
    for j=1:size_node
        if ii~=j
            xxj=coords(j,1);
            yyj=coords(j,2); 
            ll=sqrt((xxi-xxj)*(xxi-xxj)+(yyi-yyj)*(yyi-yyj));
            if ll<detla
                jj=jj+1;
                horizon_parrticle(i,jj)=j;
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEM热传导矩阵
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

KZ=zeros(size_node,size_node);
PP=zeros(size_node,1);
for i=1:size_element
    i1=elems(i,1);
    i2=elems(i,2);
    i3=elems(i,3);
    i4=elems(i,4);
    xxi=coords(i1,1);
    yyi=coords(i1,2);
    xxj=coords(i2,1);
    yyj=coords(i2,2);
    xxm=coords(i3,1);
    yym=coords(i3,2);
    xxn=coords(i4,1);
    yyn=coords(i4,2);
    ke=zeros(4);
	enode=[xxi yyi; xxj yyj; xxm yym; xxn yyn];               
    for ii11=1:2
        for jj11=1:2
            rr=pointr(ii11);
            ss=points(jj11);
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
            BB=[dNxyz(1,1) dNxyz(1,2) dNxyz(1,3) dNxyz(1,4);dNxyz(2,1) dNxyz(2,2) dNxyz(2,3) dNxyz(2,4)];
            ke=ke+BB'*[k0 0;0 k0]*BB*tt*detJ*weightr(ii11)*weights(jj11);
        end
    end
    for ii=1:4
        KZ(elems(i,ii),elems(i,1))=KZ(elems(i,ii),elems(i,1))+ke(ii,1);
        KZ(elems(i,ii),elems(i,2))=KZ(elems(i,ii),elems(i,2))+ke(ii,2);
        KZ(elems(i,ii),elems(i,3))=KZ(elems(i,ii),elems(i,3))+ke(ii,3);
        KZ(elems(i,ii),elems(i,4))=KZ(elems(i,ii),elems(i,4))+ke(ii,4);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEM热传导矩阵
KZ([PD_lizi],:)=[0];
for i=1:size_PD
    ii=PD_lizi(i,1);
    xxi=coords(ii,1);
    yyi=coords(ii,2);
    for j=1:28
        jj=horizon_parrticle(i,j);
        if jj>0
            xxj=coords(jj,1);
            yyj=coords(jj,2);
            l=sqrt((xxi-xxj)*(xxi-xxj)+(yyi-yyj)*(yyi-yyj));
            kk_PD(1,1)=c/l;
            kk_PD(1,2)=-c/l;
            KZ(ii,ii)=KZ(ii,ii)+kk_PD(1,1);
            KZ(ii,jj)=KZ(ii,jj)+kk_PD(1,2);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 施加边界条件并求解
for i=1:ii_left
    ii=boundary_left(i,1);
    KZ(ii,ii)=KZ(ii,ii)*10e8;
    PP(ii,1)=KZ(ii,ii)*0;
end
for i=1:ii_right
    ii=boundary_right(i,1);
    KZ(ii,ii)=KZ(ii,ii)*10e8;
    PP(ii,1)=KZ(ii,ii)*0;
end
for i=1:ii_bottom
    ii=boundary_bottom(i,1);
    KZ(ii,ii)=KZ(ii,ii)*10e8;
    PP(ii,1)=KZ(ii,ii)*0;
end
for i=1:ii_top
    ii=boundary_top(i,1);
    KZ(ii,ii)=KZ(ii,ii)*10e8;
    PP(ii,1)=KZ(ii,ii)*100;
end

TT=KZ\PP;
TT3=zeros(101);
ii=0;
for i=1:101
    for j=1:101
        ii=ii+1;
        TT3(i,j)=TT(ii);
    end
end




