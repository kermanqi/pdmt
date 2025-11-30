function[dimension,dettx,m,length_x,length_y,length_z,xx_min,yy_min,zz_min,xx_max,yy_max,zz_max,detla,c0,volume_particle,size_horizon_particle, ...
    boundary_condition_left,boundary_condition_right,boundary_condition_bottom,boundary_condition_top,boundary_condition_behind,boundary_condition_front,temperature_left, ...
    temperature_right,temperature_bottom,temperature_top,temperature_behind,temperature_front,heat_flux_left,heat_flux_right,heat_flux_bottom, ...
    heat_flux_top,heat_flux_behind,heat_flux_front,body_heat_source,dt,time,output_interval,density,cv,temperature_initial,nt]=function_input()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                             
% 需输入的模型参数
k0=1;          % 导热系数
dimension=1;        % 模型维度 =1 一维问题；=2 二维问题；=3 三维问题
dettx=0.5;          % 相邻粒子的距离
m=3;                % 粒子作用范围(包含几层粒子)
length_x=50;    % x方向长度
length_y=1;    % y方向长度；对于一维问题，length_y取值为截面积A
length_z=1;    % z方向长度；对于二维平面应力问题，length_z取厚度t；对于一维问题，length_z取值为1
xx_min=-25;     % x方向最小坐标
yy_min=0;     % y方向最小坐标；对于一维问题，length_y取0
zz_min=0;     % z方向最小坐标；对于一维问题，二维平面应力问题，length_z取0
dt=1.0e-9;              % 每一步的时间增量
time=3.0e-5;              % 计算总时长
output_interval=1000;     % 输出间隔，即每经过output_interval个计算步输出一次结果
density=8.0e-9;              % 材料密度
cv=1.0e2;              % 比热容
temperature_initial=100;              % 初始温度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nt=int32(time/dt);                 % 总增量步数
xx_max=length_x+xx_min;     % x方向最大坐标
yy_max=length_y+yy_min;     % y方向最大坐标
zz_max=length_z+zz_min;     % z方向最大坐标
detla=(m+0.015)*dettx;  % 粒子作用范围
if dimension==1
    area=length_y*length_z; % 截面积
	c0=2*k0/(area*detla*detla);                %微模量系数
	volume_particle=dettx*area;   %粒子体积
    if m==2
        size_horizon_particle=4;    % 如果m取2，粒子作用范围内最多包含4个粒子
    end
    if m==3
        size_horizon_particle=6;    % 如果m取3，粒子作用范围内最多包含6个粒子
    end
    if m==4
        size_horizon_particle=8;    % 如果m取4，粒子作用范围内最多包含8个粒子
    end
    if m==5
        size_horizon_particle=10;    % 如果m取5，粒子作用范围内最多包含10个粒子
    end
    if m==6
        size_horizon_particle=12;    % 如果m取6，粒子作用范围内最多包含12个粒子
    end
end
if dimension==2
    thick=length_z;
	c0=6*k0/(pi*thick*detla*detla*detla);                %微模量系数
	volume_particle=dettx*dettx*thick;   %粒子体积
    if m==2
        size_horizon_particle=12;    % 如果m取2，粒子作用范围内最多包含12个粒子
    end
    if m==3
        size_horizon_particle=28;    % 如果m取3，粒子作用范围内最多包含28个粒子
    end
    if m==4
        size_horizon_particle=48;    % 如果m取4，粒子作用范围内最多包含48个粒子
    end
    if m==5
        size_horizon_particle=80;    % 如果m取5，粒子作用范围内最多包含80个粒子
    end
    if m==6
        size_horizon_particle=112;    % 如果m取6，粒子作用范围内最多包含112个粒子
    end
end
if dimension==3
	c0=6*k0/(pi*detla*detla*detla*detla);                %微模量系数
	volume_particle=dettx*dettx*dettx;   %粒子体积
    if m==2
        size_horizon_particle=32;    % 如果m取2，粒子作用范围内最多包含32个粒子
    end
    if m==3
        size_horizon_particle=122;    % 如果m取3，粒子作用范围内最多包含122个粒子
    end
    if m==4
        size_horizon_particle=256;    % 如果m取4，粒子作用范围内最多包含256个粒子
    end
    if m==5
        size_horizon_particle=514;    % 如果m取5，粒子作用范围内最多包含514个粒子
    end
    if m==6
        size_horizon_particle=924;    % 如果m取6，粒子作用范围内最多包含924个粒子
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 需输入的边界条件
% 对于一维情况，只需要处理与left，right相关的边界条件；对于二维情况，需要处理与left，right，bottom，top相关的边界条件；对于三维情况，需要处理所有相关的边界条件
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%边界条件，=0 无边界条件施加；=1 温度边界；=2 热流边界
boundary_condition_left=2;	% 左边界
boundary_condition_right=2;	% 右边界
boundary_condition_bottom=0; 	% 下边界
boundary_condition_top=0;     % 上边界
boundary_condition_behind=0;      % 后边界
boundary_condition_front=0;      % 前边界
% 温度边界赋值，如果指定了对应的温度边界条件，需在下面给出温度数值
temperature_left=0;      % 左边界
temperature_right=0;      % 右边界
temperature_bottom=0;     % 下边界
temperature_top=0;      % 上边界
temperature_behind=0;      % 后边界
temperature_front=0;      % 前边界
% 热流边界赋值，如果指定了对应的热流条件，需在下面给出热流数值
heat_flux_left=-30;      % 左边界
heat_flux_right=-30;      % 右边界
heat_flux_bottom=0;      % 下边界
heat_flux_top=0;      % 上边界
heat_flux_behind=0;      % 后边界
heat_flux_front=0;      % 前边界
% 内热源密度赋值
body_heat_source=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


