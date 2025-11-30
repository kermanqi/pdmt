clear
nt=1000;   %计算次数
maxx=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%坐标
Coords=load('node.txt');
Elems=load('element.txt');
Coords=Coords*1.0;

ux0=0.064;
dux=0.0005;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%材料参数

E11=159960;
E22=8960;
MU12=1/3.0;
MU21=1/3*E22/E11;
G12=MU12*E22/(1.0-MU12*MU21);
Q11=E11/(1.0-MU12*MU21);
Q22=E22/(1.0-MU12*MU21);
Q12=MU12*E22/(1.0-MU12*MU21);
Q21=Q12;
Q66=G12;
DD=[Q11 Q12 0;Q21 Q22 0; 0 0 G12];
jiaodu=90;
TTT11=cos(jiaodu*pi/180)*cos(jiaodu*pi/180);
TTT12=sin(jiaodu*pi/180)*sin(jiaodu*pi/180);
TTT13=-2*sin(jiaodu*pi/180)*cos(jiaodu*pi/180);
TTT21=sin(jiaodu*pi/180)*sin(jiaodu*pi/180);
TTT22=cos(jiaodu*pi/180)*cos(jiaodu*pi/180);
TTT23=2*sin(jiaodu*pi/180)*cos(jiaodu*pi/180);
TTT31=sin(jiaodu*pi/180)*cos(jiaodu*pi/180);
TTT32=-sin(jiaodu*pi/180)*cos(jiaodu*pi/180);
TTT33=cos(jiaodu*pi/180)*cos(jiaodu*pi/180)-sin(jiaodu*pi/180)*sin(jiaodu*pi/180);
TTT_inv=[TTT11 TTT12 TTT13; TTT21 TTT22 TTT23; TTT31 TTT32 TTT33];
DD=TTT_inv*DD*(TTT_inv');
clear TTT_inv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a1=size(Coords,1);
a3=size(Elems,1);
yuneilizi=28;
%vnn=pi*tt*det*det/28;
vnn=1.0;
kz1=zeros(a1*2,a1*2);
d=zeros(a1*2,1);
v=zeros(a1*2,1);
a=zeros(a1*2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pd参数
tt=1;
dettx=1.0;
detla=3.015*dettx;
vx=dettx*dettx*tt;
ccf1=2.0*E11*(E11-E22)/((E11-E22/9.0))/vnn;
ccm=8.0*E11*E22/((E11-E22/9.0)*pi*tt*(detla^3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%边界条件
set_zuo=[1:100:9901];
set_zuo1=set_zuo*2-1;
set_zuo2=set_zuo*2;
size_zuo=size(set_zuo,2);

set_you=[100:100:10000];
set_you1=set_you*2-1;
set_you2=set_you*2;
size_you=size(set_you,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PD区域
ii=0;
for i=1:a1
    xx=Coords(i,1);
    if xx>20*dettx && xx<81*dettx
        ii=ii+1;
        node_PD(ii)=i;
    end
end
a2=size(node_PD,2);
node_PD3=[node_PD*2-1,node_PD*2];
node_PD3=sort(node_PD3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%影响域
yunei=zeros(a2,yuneilizi);
for i=1:a2
    ii1=node_PD(i);
    xx1=Coords(ii1,1);
    yy1=Coords(ii1,2);
     ii=1;
    for j=1:a1
       
        if ii1~=j
            xx2=Coords(j,1);
            yy2=Coords(j,2);
            l=sqrt((xx1-xx2)^2+(yy1-yy2)^2);
            if l<=detla
                yunei(i,ii)=j;
                ii=ii+1;
                if ii>28.5
                    break;
                end
            end
        end
    end
end       

s_juzhen=zeros(a2,yuneilizi);
c_juzhen=zeros(a2,yuneilizi);
for i=1:a2
    xx=Coords(node_PD(i),1);
    yy=Coords(node_PD(i),2);
     ii=1;
    for j=1:yuneilizi
       if yunei(i,j)~=0
            ii2=yunei(i,j);
            xx1=Coords(ii2,1);
            yy1=Coords(ii2,2);
            if abs((xx-xx1))<(0.1*1.0) && abs((yy-yy1))>(0.9*1.0)
                 s_juzhen(i,j)=0.01775862;
                 c_juzhen(i,j)=ccm+ccf1/(12*dettx);
            else
                 s_juzhen(i,j)=0.0067077;
                 c_juzhen(i,j)=ccm;
            end
        end
    end
end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%预制裂纹
damage=zeros(a2,yuneilizi);
yyshuju=50.5*dettx;
for i=1:a2
    ii1=node_PD(i);
    yy=Coords(ii1,1);
    xx=Coords(ii1,2);
    for j=1:yuneilizi
        if yunei(i,j)~=0
            ii2=yunei(i,j);
            yy1=Coords(ii2,1);
            xx1=Coords(ii2,2);
        %    l=sqrt((xx-xx1)*(xx-xx1)+(yy-yy1)*(yy-yy1));
       %     if l<detla
                if yy>yyshuju && yy1<yyshuju
                        detx=abs(xx-xx1);
                        if xx>xx1
                            detx2=xx1+detx*abs(yyshuju-yy1)/abs(yy1-yy);
                        end
                        if xx<=xx1
                            detx2=xx+detx*abs(yyshuju-yy)/abs(yy1-yy);
                        end
                        if detx2<(71+2.5) && detx2>(30-2.5)
                            damage(i,j)=1;
                        end
                end
                 if yy1>yyshuju && yy<yyshuju
                        detx=abs(xx-xx1);
                        if xx>xx1
                            detx2=xx1+detx*abs(yyshuju-yy1)/abs(yy1-yy);
                        end
                        if xx<=xx1
                            detx2=xx+detx*abs(yyshuju-yy)/abs(yy1-yy);
                        end
                        if detx2<(71+2.5) && detx2>(30-2.5)
                            damage(i,j)=1;
                        end
                 end  
        %    end
        end
    end
end
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
for i=1:a3  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%刚度矩阵，质量矩阵总装
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
        ij0=Elems(i,ij);
        kz1(ij0*2-1,i1*2-1)=kz1(ij0*2-1,i1*2-1)+ke(ij*2-1,1);
        kz1(ij0*2-1,i1*2-0)=kz1(ij0*2-1,i1*2-0)+ke(ij*2-1,2);
        kz1(ij0*2-0,i1*2-1)=kz1(ij0*2-0,i1*2-1)+ke(ij*2-0,1);
        kz1(ij0*2-0,i1*2-0)=kz1(ij0*2-0,i1*2-0)+ke(ij*2-0,2);
                                             
        kz1(ij0*2-1,i2*2-1)=kz1(ij0*2-1,i2*2-1)+ke(ij*2-1,3);
        kz1(ij0*2-1,i2*2-0)=kz1(ij0*2-1,i2*2-0)+ke(ij*2-1,4);
        kz1(ij0*2-0,i2*2-1)=kz1(ij0*2-0,i2*2-1)+ke(ij*2-0,3);
        kz1(ij0*2-0,i2*2-0)=kz1(ij0*2-0,i2*2-0)+ke(ij*2-0,4);
                                             
        kz1(ij0*2-1,i3*2-1)=kz1(ij0*2-1,i3*2-1)+ke(ij*2-1,5);
        kz1(ij0*2-1,i3*2-0)=kz1(ij0*2-1,i3*2-0)+ke(ij*2-1,6);
        kz1(ij0*2-0,i3*2-1)=kz1(ij0*2-0,i3*2-1)+ke(ij*2-0,5);
        kz1(ij0*2-0,i3*2-0)=kz1(ij0*2-0,i3*2-0)+ke(ij*2-0,6);
                                                     
        kz1(ij0*2-1,i4*2-1)=kz1(ij0*2-1,i4*2-1)+ke(ij*2-1,7);
        kz1(ij0*2-1,i4*2-0)=kz1(ij0*2-1,i4*2-0)+ke(ij*2-1,8);
        kz1(ij0*2-0,i4*2-1)=kz1(ij0*2-0,i4*2-1)+ke(ij*2-0,7);
        kz1(ij0*2-0,i4*2-0)=kz1(ij0*2-0,i4*2-0)+ke(ij*2-0,8);
    end        
end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


maxx=0;
sunshang=zeros(a1,1);


  kz1([node_PD3],:)=0;


for itime=1:nt
    
    itime


    pz1=zeros(a1*2,1);


    kz2=zeros(a1*2,a1*2);
for ii=1:a2
    ii1=node_PD(ii);
    xxi=Coords(ii1,1);
    yyi=Coords(ii1,2);
    for j=1:28
        jj=yunei(ii,j);
        if jj>0
            if damage(ii,j)<0.5
                xxj=Coords(jj,1);
                yyj=Coords(jj,2);
                l=sqrt((xxi-xxj)*(xxi-xxj)+(yyi-yyj)*(yyi-yyj));
                dex=(xxj-xxi);
                dey=(yyj-yyi);
                c=c_juzhen(ii,j);
                kk_PD=zeros(2,4);
                kk_PD(1,1)=vx*vx*c*dex*dex/l;
                kk_PD(1,2)=vx*vx*c*dex*dey/l;
                kk_PD(1,3)=-vx*vx*c*dex*dex/l;
                kk_PD(1,4)=-vx*vx*c*dex*dey/l;
                kk_PD(2,1)=vx*vx*c*dex*dey/l;
                kk_PD(2,2)=vx*vx*c*dey*dey/l;
                kk_PD(2,3)=-vx*vx*c*dex*dey/l;
                kk_PD(2,4)=-vx*vx*c*dey*dey/l;
                
                kz2(ii1*2-1,ii1*2-1)=kz2(ii1*2-1,ii1*2-1)+kk_PD(1,1);
                kz2(ii1*2-1,ii1*2)=kz2(ii1*2-1,ii1*2)+kk_PD(1,2);
                kz2(ii1*2-1,jj*2-1)=kz2(ii1*2-1,jj*2-1)+kk_PD(1,3);
                kz2(ii1*2-1,jj*2)=kz2(ii1*2-1,jj*2)+kk_PD(1,4);
                kz2(ii1*2,ii1*2-1)=kz2(ii1*2,ii1*2-1)+kk_PD(2,1);
                kz2(ii1*2,ii1*2)=kz2(ii1*2,ii1*2)+kk_PD(2,2);
                kz2(ii1*2,jj*2-1)=kz2(ii1*2,jj*2-1)+kk_PD(2,3);
                kz2(ii1*2,jj*2)=kz2(ii1*2,jj*2)+kk_PD(2,4);
            end
            
        end
    end
end

kz2=kz1+kz2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NEWMARK求解
  %  t=(itime-1)*dt;
for i=1:size_zuo
	ii=set_zuo(i);
    kz2(ii*2-1,ii*2-1)=kz2(ii*2-1,ii*2-1)*1.0e8;
    pz1(ii*2-1,1)=kz2(ii*2-1,ii*2-1)*(-ux0-itime*dux);
  
end
for i=1:size_you
	ii=set_you(i);
    kz2(ii*2-1,ii*2-1)=kz2(ii*2-1,ii*2-1)*1.0e8;
    pz1(ii*2-1,1)=kz2(ii*2-1,ii*2-1)*(ux0+itime*dux);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


uz1=kz2\pz1;
    
    
for i=1:a2
    ii1=node_PD(i);
    xx=Coords(ii1,1);
    yy=Coords(ii1,2);
    weiyix=uz1(ii1*2-1,1);
    weiyiy=uz1(ii1*2,1);
    damage1=0;
    damage2=0;
    for j=1:28
        if yunei(i,j)~=0
            ii2=yunei(i,j);

                 if abs(damage(i,j))>0.9
                    damage1=damage1+1;
                  end
              if abs(damage(i,j))<0.1
                xx1=Coords(ii2,1);
                yy1=Coords(ii2,2);
                l=sqrt((xx-xx1)*(xx-xx1)+(yy-yy1)*(yy-yy1));
                weiyix1=uz1(ii2*2-1,1);           
                weiyiy1=uz1(ii2*2,1);
                detu=sqrt((weiyix+xx-weiyix1-xx1)*(weiyix+xx-weiyix1-xx1)+(weiyiy+yy-weiyiy1-yy1)*(weiyiy+yy-weiyiy1-yy1));
                panju=(detu-l)/l;
                if panju>maxx
                    maxx=panju;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %判断损伤
                s0=s_juzhen(i,j);
                if panju>s0
                    damage1=damage1+1;
                    damage(i,j)=1;
                end
              end
       damage2=damage2+1;
        end
    end
sunshang(ii1,1)=damage1/damage2;
end







if mod(itime,1)==0
sunshang33=zeros(100,100);
uux=zeros(100,100);
uuy=zeros(100,100);
ii=1;
for i=1:100
    for j=1:100
        sunshang33(i,j)=sunshang(ii);
        uux(i,j)=uz1(ii*2-1);
        uuy(i,j)=uz1(ii*2);
        ii=ii+1;
    end
end



qqq1='sunshang';
qqq2=num2str(itime);
qqq3=strcat(qqq1,qqq2);
fname=[num2str(qqq3), '.xls'];
xlswrite(fname,sunshang33);

qqq1='ux';
qqq3=strcat(qqq1,qqq2);
fname=[num2str(qqq3), '.xls'];
xlswrite(fname,uux);

qqq1='uy';
qqq3=strcat(qqq1,qqq2);
fname=[num2str(qqq3), '.xls'];
xlswrite(fname,uuy);
end

end











