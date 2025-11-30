clear
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PD粒子
ii=0;
for i=1:106
    for j=1:100
        ii=ii+1;
        node(ii,1)=i-50.5-3;
        node(ii,2)=j-50.5;
    end
end
node=node*0.5;
size_node=ii;
maxx=0;
nt=1000;
ux0=0.033;
dux=0.0002;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%材料参数
E11=159960;
E22=8960;
MU12=0.35;
MU21=0.35*E22/E11;
G12=MU12*E22/(1.0-MU12*MU21);
Q11=E11/(1.0-MU12*MU21);
Q22=E22/(1.0-MU12*MU21);
Q12=MU12*E22/(1.0-MU12*MU21);
Q21=Q12;
Q66=G12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PD参数
tt=1;                   %厚度
dettx=0.5;            %单位长度
vx=dettx*dettx*tt;  %体积
detla=3.015*dettx;     %delta
yuneilizi=28;
dd=2/(pi*tt*detla*detla*detla);
aa=(Q12-Q66)/2;
bbf=(Q11-Q12-2*Q66)/(2*detla*12);
bbt=(Q22-Q12-2*Q66)/(2*detla*12);
bbft=6*Q66/(pi*tt*detla*detla*detla*detla);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%影响域
yunei=zeros(size_node,yuneilizi);
for i=1:size_node
    xx1=node(i,1);
    yy1=node(i,2);
     ii=1;
    for j=1:size_node
        if i~=j
            xx2=node(j,1);
            yy2=node(j,2);
            l=sqrt((xx1-xx2)^2+(yy1-yy2)^2);
            if l<=detla
                yunei(i,ii)=j;
                ii=ii+1;
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 键的分类
fenlei=zeros(size_node,yuneilizi);
for i=1:size_node
    xx1=node(i,1);
    yy1=node(i,2);
    for j=1:yuneilizi
        jj=yunei(i,j);
        if jj>0
            xx2=node(jj,1);
            yy2=node(jj,2);
            if abs((xx1-xx2)-(yy1-yy2))<0.1*dettx % 45度方向
                fenlei(i,j)=1;
            end
            if abs((xx1-xx2)+(yy1-yy2))<0.1*dettx % -45度方向
                fenlei(i,j)=2;
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%预制裂纹
damage=zeros(size_node,yuneilizi);
yyshuju=0;
for i=1:size_node
    yy=node(i,1);
    xx=node(i,2);
    for j=1:yuneilizi
        if yunei(i,j)~=0
            ii2=yunei(i,j);
            yy1=node(ii2,1);
            xx1=node(ii2,2);
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
                        if abs(detx2)<15
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
                        if abs(detx2)<15
                            damage(i,j)=1;
                        end
                 end  
        %    end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%边界条件
ii_left=0;
ii_right=0;
boundary_left=[];
boundary_right=[];
for i=1:size_node
    xx=node(i,1);
    if xx<-25
        ii_left=ii_left+1;
        boundary_left(ii_left,1)=i;
    end
    if xx>25
        ii_right=ii_right+1;
        boundary_right(ii_right,1)=i;
    end
end
size_left=ii_left;
size_right=ii_right;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xishu11=2*aa*dd*dd*detla*detla*vx*vx*vx;
for itime=1:nt
    itime
    kz1=zeros(size_node*2,size_node*2);
    pz1=zeros(size_node*2,1);
for i=1:size_node
    xx=node(i,1);
    yy=node(i,2);
    for j=1:yuneilizi
        if yunei(i,j)~=0
            jj1=yunei(i,j);
            xx1=node(jj1,1);
            yy1=node(jj1,2);
            fenmu1=sqrt((xx1-xx)*(xx1-xx)+(yy1-yy)*(yy1-yy));
            detx=(xx1-xx);
            dety=(yy1-yy);
            um1=1;
            if damage(i,j)>0.5
                um1=0.00001;
            end
            for kk1=1:yuneilizi
                if yunei(i,kk1)~=0
                    kk2=yunei(i,kk1);
                    xx2=node(kk2,1);
                    yy2=node(kk2,2);
                    fenmu2=sqrt((xx2-xx)*(xx2-xx)+(yy2-yy)*(yy2-yy));
                    um2=1;
                    if damage(i,kk1)>0.5
                        um2=0.00001;
                    end
                    xishu1=um1*um2*xishu11/fenmu1/fenmu2;
                    kk_PD(1,1)=xishu1*detx*detx;
                    kk_PD(1,2)=xishu1*detx*dety;
                    kk_PD(2,1)=xishu1*dety*detx;
                    kk_PD(2,2)=xishu1*dety*dety;
                    kk_PD(1,3)=-xishu1*detx*detx;
                    kk_PD(1,4)=-xishu1*detx*dety;
                    kk_PD(2,3)=-xishu1*dety*detx;
                    kk_PD(2,4)=-xishu1*dety*dety;
                    
                    kz1(i*2-1,i*2-1)=kz1(i*2-1,i*2-1)+kk_PD(1,1);
                    kz1(i*2-1,i*2)=kz1(i*2-1,i*2)+kk_PD(1,2);
                    kz1(i*2-1,kk2*2-1)=kz1(i*2-1,kk2*2-1)+kk_PD(1,3);
                    kz1(i*2-1,kk2*2)=kz1(i*2-1,kk2*2)+kk_PD(1,4);
                    kz1(i*2,i*2-1)=kz1(i*2,i*2-1)+kk_PD(2,1);
                    kz1(i*2,i*2)=kz1(i*2,i*2)+kk_PD(2,2);
                    kz1(i*2,kk2*2-1)=kz1(i*2,kk2*2-1)+kk_PD(2,3);
                    kz1(i*2,kk2*2)=kz1(i*2,kk2*2)+kk_PD(2,4);
                end
            end

            for kk1=1:yuneilizi
                if yunei(jj1,kk1)~=0
                    kk2=yunei(jj1,kk1);
                    xx2=node(kk2,1);
                    yy2=node(kk2,2);
                    um2=1;
                    if damage(i,kk1)>0.5
                        um2=0.00001;
                    end
                    fenmu2=sqrt((xx2-xx1)*(xx2-xx1)+(yy2-yy1)*(yy2-yy1));
                    xishu1=um1*um2*xishu11/fenmu1/fenmu2;
                    kk_PD(1,1)=xishu1*detx*detx;
                    kk_PD(1,2)=xishu1*detx*dety;
                    kk_PD(2,1)=xishu1*dety*detx;
                    kk_PD(2,2)=xishu1*dety*dety;
                    kk_PD(1,3)=-xishu1*detx*detx;
                    kk_PD(1,4)=-xishu1*detx*dety;
                    kk_PD(2,3)=-xishu1*dety*detx;
                    kk_PD(2,4)=-xishu1*dety*dety;
                
                    kz1(i*2-1,kk2*2-1)=kz1(i*2-1,kk2*2-1)+kk_PD(1,1);
                    kz1(i*2-1,kk2*2)=kz1(i*2-1,kk2*2)+kk_PD(1,2);
                    kz1(i*2-1,jj1*2-1)=kz1(i*2-1,jj1*2-1)+kk_PD(1,3);
                    kz1(i*2-1,jj1*2)=kz1(i*2-1,jj1*2)+kk_PD(1,4);
                    kz1(i*2,kk2*2-1)=kz1(i*2,kk2*2-1)+kk_PD(2,1);
                    kz1(i*2,kk2*2)=kz1(i*2,kk2*2)+kk_PD(2,2);
                    kz1(i*2,jj1*2-1)=kz1(i*2,jj1*2-1)+kk_PD(2,3);
                    kz1(i*2,jj1*2)=kz1(i*2,jj1*2)+kk_PD(2,4);
                end
            end
            
            bb=bbft;
            if fenlei(i,j)==1
                bb=bbft+bbf;
            end
            if fenlei(i,j)==2
                bb=bbft+bbt;
            end
            xishu22=4.0*bb*detla*vx*vx;
            xishu2=um1*xishu22/fenmu1;
            kk_PD(1,1)=xishu2*detx*detx;
            kk_PD(1,2)=xishu2*detx*dety;
            kk_PD(2,1)=xishu2*dety*detx;
            kk_PD(2,2)=xishu2*dety*dety;
            kk_PD(1,3)=-xishu2*detx*detx;
            kk_PD(1,4)=-xishu2*detx*dety;
            kk_PD(2,3)=-xishu2*dety*detx;
            kk_PD(2,4)=-xishu2*dety*dety;
                
            kz1(i*2-1,i*2-1)=kz1(i*2-1,i*2-1)+kk_PD(1,1);
            kz1(i*2-1,i*2)=kz1(i*2-1,i*2)+kk_PD(1,2);
            kz1(i*2-1,jj1*2-1)=kz1(i*2-1,jj1*2-1)+kk_PD(1,3);
            kz1(i*2-1,jj1*2)=kz1(i*2-1,jj1*2)+kk_PD(1,4);
            kz1(i*2,i*2-1)=kz1(i*2,i*2-1)+kk_PD(2,1);
            kz1(i*2,i*2)=kz1(i*2,i*2)+kk_PD(2,2);
            kz1(i*2,jj1*2-1)=kz1(i*2,jj1*2-1)+kk_PD(2,3);
            kz1(i*2,jj1*2)=kz1(i*2,jj1*2)+kk_PD(2,4);
        end
    end
 end
 
for i=1:size_left
	ii=boundary_left(i,1);
    kz1(ii*2-1,ii*2-1)=kz1(ii*2-1,ii*2-1)*1.0e8;
    pz1(ii*2-1,1)=kz1(ii*2-1,ii*2-1)*(-ux0-itime*dux);
  
end
for i=1:size_right
	ii=boundary_right(i,1);
    kz1(ii*2-1,ii*2-1)=kz1(ii*2-1,ii*2-1)*1.0e8;
    pz1(ii*2-1,1)=kz1(ii*2-1,ii*2-1)*(ux0+itime*dux);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uz1=kz1\pz1;
sunshang=zeros(size_node,1);
for i=1:size_node
    xx=node(i,1);
    yy=node(i,2);
    weiyix=uz1(i*2-1,1);
    weiyiy=uz1(i*2,1);
    damage1=0;
    damage2=0;
    for j=1:yuneilizi
        if yunei(i,j)~=0
            ii2=yunei(i,j);
            if abs(damage(i,j))>0.9
                damage1=damage1+1;
            end
            if abs(damage(i,j))<0.1
                xx1=node(ii2,1);
                yy1=node(ii2,2);
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
                s0=0.0067077;
                if fenlei(i,j)==1
                    s0=0.01775862;
                end
                if panju>s0
                    damage1=damage1+1;
                    damage(i,j)=1;
                end
              end
       damage2=damage2+1;
        end
    end
    sunshang(i,1)=damage1/damage2;
end

ux=zeros(100,106);
ii=1;
for i=1:106
    for j=1:100
        ux(j,i)=uz1(ii*2-1,1);
        ii=ii+1;
    end
end

uy=zeros(100,106);
ii=1;
for i=1:106
    for j=1:100
        uy(j,i)=uz1(ii*2,1);
        ii=ii+1;
    end
end

sun=zeros(100,106);
ii=1;
for i=1:106
    for j=1:100
        sun(j,i)=sunshang(ii,1);
        ii=ii+1;
    end
end
qqq1='sunshang';
qqq2=num2str(itime);
qqq3=strcat(qqq1,qqq2);
fname=[num2str(qqq3), '.xls'];
xlswrite(fname,sun);

qqq1='ux';
qqq3=strcat(qqq1,qqq2);
fname=[num2str(qqq3), '.xls'];
xlswrite(fname,ux);

qqq1='uy';
qqq3=strcat(qqq1,qqq2);
fname=[num2str(qqq3), '.xls'];
xlswrite(fname,uy);
end













