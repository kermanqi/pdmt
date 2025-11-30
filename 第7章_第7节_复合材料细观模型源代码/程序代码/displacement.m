%体积分数35.0%
clear
tic
ntime=10000;
tt=1;
radius=0.102;
dettx=0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%节点坐标
size_node=256*250;
ii=0;
node=zeros(size_node,2);
for i=1:256
    for j=1:250
        ii=ii+1;
        node(ii,1)=(i-128-0.5)*dettx;
        node(ii,2)=(j-125-0.5)*dettx;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
size_node=size(node,1);
set_circle=load('circle.txt');
size_circle=size(set_circle,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
number_horizon_particle=28;
detla=3.001*dettx;
vx=dettx*dettx*tt;
%材料参数
%基体
tt=2;                   %厚度
Em=3500;               %弹性模量
cm=9*Em/(pi*tt*detla*detla*detla);                %微模量系数
%纤维束
Ef=28000;               %弹性模量
cf=9*Ef/(pi*tt*detla*detla*detla);                %微模量系数
%界面
Ej=2*Em*Ef/(Em+Ef);
cj=9*Ej/(pi*tt*detla*detla*detla);                %微模量系数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%域内粒子
horizon_particle=zeros(size_node,number_horizon_particle);
property_bond=zeros(size_node,number_horizon_particle);
damage_particle=zeros(size_node,number_horizon_particle);
for i=1:size_node
    xx1=node(i,1);
    yy1=node(i,2);
    property_i=1; %基体
    for ij=1:size_circle
        xx_circle=set_circle(ij,1);
        yy_circle=set_circle(ij,2);
        length=sqrt((xx1-xx_circle)*(xx1-xx_circle)+(yy1-yy_circle)*(yy1-yy_circle));
        if length<radius
            property_i=2; %纤维
            break;
        end
    end
    jj=0;
    for j=1:size_node
        xx2=node(j,1);
        yy2=node(j,2);
        l=sqrt((xx2-xx1)*(xx2-xx1)+(yy2-yy1)*(yy2-yy1));
        if l<detla && l>0
            jj=jj+1;
            horizon_particle(i,jj)=j;
            property_j=1; %基体
            for ij=1:size_circle
                xx_circle=set_circle(ij,1);
                yy_circle=set_circle(ij,2);
                length=sqrt((xx2-xx_circle)*(xx2-xx_circle)+(yy2-yy_circle)*(yy2-yy_circle));
                if length<radius
                    property_j=2; %纤维
                    break;
                end
            end
            if property_i==1 && property_j==1 %基体
                property_bond(i,jj)=1;
            end
            if property_i==2 && property_j==2 %纤维
                property_bond(i,jj)=2;
            end
            if property_i==1 && property_j==2 %混合
                property_bond(i,jj)=3;
            end
            if property_i==2 && property_j==1 %混合
                property_bond(i,jj)=3;
            end
            
        end
        if jj>=number_horizon_particle
            break
        end
    end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%边界粒子
bianjie_zuo=[];
ii=0;
for i=1:size_node
    xx=node(i,1);
    if xx<-1.25
        ii=ii+1;
        bianjie_zuo(ii,1)=i;
    end
end

bianjie_you=[];
ii=0;
for i=1:size_node
    xx=node(i,1);
    if xx>1.25
        ii=ii+1;
        bianjie_you(ii,1)=i;
    end
end

size_zuo=size(bianjie_zuo,1);
size_you=size(bianjie_you,1);
global_load_vector=zeros(size_node*2,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
size_vector=0;
    for i=1:size_node                       % 对粒子循环
        xx_i=node(i,1);                     % 获取第i个粒子的x方向坐标
        yy_i=node(i,2);                     % 获取第i个粒子的y方向坐标
        for j=1:number_horizon_particle           % 对第i个粒子域内的粒子进行循环
            if damage_particle(i,j)<0.5         % 如果第i个粒子与第jj个粒子组成的键没有发生破坏，则执行条件语句中的内容
                jj=horizon_particle(i,j);           % 第i个粒子域内的第j个粒子的粒子编号为jj
                if jj>0                             % 如果第i个粒子域内的第j个粒子存在，则执行条件语句中的内容
                    xx_j=node(jj,1);            % 获取第jj个粒子的x方向坐标
                    yy_j=node(jj,2);            % 获取第jj个粒子的y方向坐标
                    length_i_jj=sqrt((xx_j-xx_i)*(xx_j-xx_i)+(yy_j-yy_i)*(yy_j-yy_i));   % 初始构型下粒子i与粒子jj之间的距离为length_i_jj
                    detl_x=(xx_j-xx_i)/length_i_jj;                                     % 计算式子()中的detx
                    detl_y=(yy_j-yy_i)/length_i_jj;                                     % 计算式子()中的dety
                    property=property_bond(i,j);
                    if property==1
                        c0=cm;
                    end
                    if property==2
                        c0=cf;
                    end
                    if property==3
                        c0=cj;
                    end
                    kk_PD=c0/length_i_jj*[detl_x*detl_x detl_x*detl_y -detl_x*detl_x -detl_x*detl_y; detl_y*detl_x detl_y*detl_y -detl_y*detl_x -detl_y*detl_y]*vx*vx; % 粒子i与粒子jj组成键的刚度矩阵
                    i1=[i*2-1,i*2];                 % 需要组装到总体刚度矩阵中的行向量
                    j1=[i*2-1,i*2,jj*2-1,jj*2];     % 需要组装到总体刚度矩阵中的列向量
                    for ii3=1:2                     % 对行向量进行循环
                        for jj3=1:4                 % 对列向量进行循环
                            size_vector=size_vector+1;  % 变量+1
                            vector_x(size_vector)=i1(ii3);  % 行向量元素赋值
                            vector_y(size_vector)=j1(jj3);  % 列向量元素赋值
                            stiffness_matrix(size_vector)=kk_PD(ii3,jj3);  % 总体刚度矩阵元素赋值
                        end
                    end
                end
            end
        end
    end
    global_stiffness_matrix=sparse(vector_x,vector_y,stiffness_matrix,size_node*2,size_node*2); % 合材总体刚度矩阵
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 左边界

        for i=1:size_zuo
            ii=bianjie_zuo(i,1);
            global_stiffness_matrix(ii*2-1,ii*2-1)=global_stiffness_matrix(ii*2-1,ii*2-1)*1.0e8;
            global_load_vector(ii*2-1,1)=global_stiffness_matrix(ii*2-1,ii*2-1)*(-0.5);
        end
        for i=1:size_you
            ii=bianjie_you(i,1);
            global_stiffness_matrix(ii*2-1,ii*2-1)=global_stiffness_matrix(ii*2-1,ii*2-1)*1.0e8;
            global_load_vector(ii*2-1,1)=global_stiffness_matrix(ii*2-1,ii*2-1)*(0.5);
        end

global_displacement=global_stiffness_matrix\global_load_vector;


displacement_x=zeros(250,256);
ii=0;
for j=1:256
    for i=1:250
        ii=ii+1;
        displacement_x(i,j)=global_displacement(ii*2-1,1);   
    end
end

displacement_y=zeros(250,256);
ii=0;
for j=1:256
    for i=1:250
        ii=ii+1;
        displacement_y(i,j)=global_displacement(ii*2,1);   
    end
end
