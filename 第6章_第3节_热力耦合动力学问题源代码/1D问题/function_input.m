function[dimension,dettx,m,length_x,length_y,length_z,xx_min,yy_min,zz_min,xx_max,yy_max,zz_max,detla,c0,volume_particle,size_horizon_particle, ...
    boundary_condition_left,boundary_condition_right,boundary_condition_bottom,boundary_condition_top,boundary_condition_behind,boundary_condition_front,velocity_left, ...
    velocity_right,velocity_bottom,velocity_top,velocity_behind,velocity_front,uniform_load_left,uniform_load_right,uniform_load_bottom, ...
    uniform_load_top,uniform_load_behind,uniform_load_front,body_force,c0_temperature,boundary_condition_left_temperature,boundary_condition_right_temperature, ...
    boundary_condition_bottom_temperature,boundary_condition_top_temperature,boundary_condition_behind_temperature,boundary_condition_front_temperature,temperature_left, ...
    temperature_right,temperature_bottom,temperature_top,temperature_behind,temperature_front,heat_flux_left,heat_flux_right,heat_flux_bottom, ...
    heat_flux_top,heat_flux_behind,heat_flux_front,body_heat_source,alf,temperature_initial,cv,dt,output_interval,density,nt]=function_input()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                             
% 需输入的模型参数
E0=210000;          % 弹性模量
k0=1;          % 导热系数
alf=1e-6;   % 热膨胀系数
dimension=1;        % 模型维度 =1 一维问题；=2 二维平面应力问题；=3 三维问题
dettx=0.5;          % 相邻粒子的距离
m=3;                % 粒子作用范围(包含几层粒子)
length_x=50;    % x方向长度
length_y=1;    % y方向长度；对于一维问题，length_y取值为截面积A
length_z=1;    % z方向长度；对于二维平面应力问题，length_z取厚度t；对于一维问题，length_z取值为1
xx_min=-25;     % x方向最小坐标
yy_min=0;     % y方向最小坐标；对于一维问题，length_y取0
zz_min=0;     % z方向最小坐标；对于一维问题，二维平面应力问题，length_z取0
temperature_initial=100; % 模型初始温度
cv=1.0e2;              % 比热容
dt=1.0e-8;              % 每一步的时间增量
time=5.0e-4;              % 计算总时长
output_interval=10000;              % 输出间隔，即每经过output_interval个计算步输出一次结果
density=8.0e-9;              % 材料密度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nt=round(time/dt);                 % 总增量步数
xx_max=length_x+xx_min;     % x方向最大坐标
yy_max=length_y+yy_min;     % y方向最大坐标
zz_max=length_z+zz_min;     % z方向最大坐标
detla=(m+0.015)*dettx;  % 粒子作用范围
if dimension==1
    area=length_y*length_z; % 截面积
	c0=2*E0/(area*detla*detla);                % 微模量系数
    c0_temperature=2*k0/(area*detla*detla);    % 微导热系数
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
    thick=length_z;  % 模型厚度
	c0=9*E0/(pi*thick*detla*detla*detla);                % 微模量系数
    c0_temperature=6*k0/(pi*thick*detla*detla*detla);    % 微导热系数
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
	c0=12*E0/(pi*detla*detla*detla*detla);                % 微模量系数
    c0_temperature=6*k0/(pi*detla*detla*detla*detla);     % 微导热系数
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
% 需输入的力学边界条件
% 对于一维情况，只需要处理与left，right相关的边界条件；对于二维情况，需要处理与left，right，bottom，top相关的边界条件；对于三维情况，需要处理所有相关的边界条件
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%边界条件，其中三个分量分别表示x,y,z方向，=0 无边界条件施加；=1 位移边界条件；=2 均匀分布荷载；对于一维情况，y,z方向分量不做处理；对于二维情况，z方向分量不做处理
boundary_condition_left=[1,0,0];	% 左边界
boundary_condition_right=[0,0,0];	% 右边界
boundary_condition_bottom=[0,0,0]; 	% 下边界
boundary_condition_top=[0,0,0];     % 上边界
boundary_condition_behind=[0,0,0];      % 后边界
boundary_condition_front=[0,0,0];      % 前边界
% 位移边界赋值，其中三个分量分别表示x,y,z方向，如果指定了对应方向的位移边界条件，需在下面给出位移数值
velocity_left=[0.0,0,0];      % 左边界
velocity_right=[0,0,0];      % 右边界
velocity_bottom=[0,0,0];     % 下边界
velocity_top=[0,0,0];      % 上边界
velocity_behind=[0,0,0];      % 后边界
velocity_front=[0,0,0];      % 前边界
% 均布荷载密度赋值，如果指定了对应方向的荷载条件，需在下面给出荷载密度数值
uniform_load_left=[0,0,0];      % 左边界
uniform_load_right=[0,0,0];      % 右边界
uniform_load_bottom=[0,0,0];      % 下边界
uniform_load_top=[0,0,0];      % 上边界
uniform_load_behind=[0,0,0];      % 后边界
uniform_load_front=[0,0,0];      % 前边界
% 体力密度赋初值，如果指定了对应方向的荷载条件，需在下面给体力在各个方向的分量
body_force=[0,0,0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 需输入的热传导边界条件
% 对于一维情况，只需要处理与left，right相关的边界条件；对于二维情况，需要处理与left，right，bottom，top相关的边界条件；对于三维情况，需要处理所有相关的边界条件
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%边界条件，=0 无边界条件施加；=1 温度边界；=2 热流边界
boundary_condition_left_temperature=1;	% 左边界
boundary_condition_right_temperature=2;	% 右边界
boundary_condition_bottom_temperature=0; 	% 下边界
boundary_condition_top_temperature=0;     % 上边界
boundary_condition_behind_temperature=0;      % 后边界
boundary_condition_front_temperature=0;      % 前边界
% 温度边界赋值，如果指定了对应的温度边界条件，需在下面给出温度数值
temperature_left=100;      % 左边界
temperature_right=0;      % 右边界
temperature_bottom=0;     % 下边界
temperature_top=0;      % 上边界
temperature_behind=0;      % 后边界
temperature_front=0;      % 前边界
% 热流边界赋值，如果指定了对应的热流条件，需在下面给出热流数值
heat_flux_left=0;      % 左边界
heat_flux_right=5;      % 右边界
heat_flux_bottom=0;      % 下边界
heat_flux_top=0;      % 上边界
heat_flux_behind=0;      % 后边界
heat_flux_front=0;      % 前边界
% 内热源密度赋值
body_heat_source=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


