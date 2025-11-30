function[dimension,dettx,m,length_x,length_y,length_z,xx_min,yy_min,zz_min,xx_max,yy_max,zz_max,detla,cf,cm,volume_particle,size_horizon_particle, ...
    boundary_condition_left,boundary_condition_right,boundary_condition_bottom,boundary_condition_top,boundary_condition_behind,boundary_condition_front,displacement_left, ...
    displacement_right,displacement_bottom,displacement_top,displacement_behind,displacement_front,uniform_load_left,uniform_load_right,uniform_load_bottom, ...
    uniform_load_top,uniform_load_behind,uniform_load_front,body_force,theta]=function_input()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                             
% 需输入的模型参数
E1=195960;          % 纤维方向弹性模量
E2=8960;            % 垂直纤维方向弹性模量
theta=pi/4;         % 纤维方向与x方向的夹角，theta可以取0，pi/4，pi/2
dimension=2;        % 对于复合材料单层板，dimension取2
dettx=0.5;          % 相邻粒子的距离
m=3;                % 粒子作用范围(包含几层粒子)
length_x=50;    % x方向长度
length_y=50;    % y方向长度
length_z=1;    % 对于复合材料单层板，length_z取厚度thick
thick=length_z;    % 单层板厚度厚
xx_min=-25;     % x方向最小坐标
yy_min=-25;     % y方向最小坐标
zz_min=0;     % 对于复合材料单层板，zz_min取0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xx_max=length_x+xx_min;     % x方向最大坐标
yy_max=length_y+yy_min;     % y方向最大坐标
zz_max=length_z+zz_min;     % z方向最大坐标
detla=(m+0.015)*dettx;  % 粒子作用范围
cm=8*E1*E2/((E1-E2/9)*pi*thick*detla*detla*detla);	% 基体键微模量系数
volume_particle=dettx*dettx*thick;   %粒子体积
cf_parameter1=2*E1*(E1-E2)/(E1-E2/9);	% 纤维键微模量系数的参数1
cf_parameter2=0;	% 纤维键微模量系数的参数2赋初值为0
size_bond_0_90=[2,3,4,5,6];	% 对于角度theta=0和theta=pi/2，m取不同数值时，模型中包含不同长度的纤维键的数量
size_bond_45=[1,2,2,3,4];	% 对于角度theta=pi/4，m取不同数值时，模型中包含不同长度的纤维键的数量
if theta==pi/4  % theta=pi/4，则执行以下内容
	for i=1:size_bond_45(m-1) % m取不同数值时，对模型中包含不同长度的键的数量进行循环
        cf_parameter2=cf_parameter2+2*sqrt(2)*dettx*i*volume_particle; % 对纤维键微模量系数的参数2进行累加，模型包含每种长度键的个数均为2个，因此为2*sqrt(2)*dettx*i*volume_particle
    end
    cf=cf_parameter1/cf_parameter2;	% 纤维键微模量系数
else  % 如果theta=0或theta=pi/2，则执行以下内容
	for i=1:size_bond_0_90(m-1) % m取不同数值时，对模型中包含不同长度的键的数量进行循环
        cf_parameter2=cf_parameter2+2*dettx*i*volume_particle;  % 对纤维键微模量系数的参数2进行累加，模型包含每种长度键的个数均为2个，因此为2*dettx*i*volume_particle
	end
    cf=cf_parameter1/cf_parameter2;	% 纤维键微模量系数
end 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 需输入的边界条件
% 对于复合材料单层板，需要处理与left，right，bottom，top相关的边界条件
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 边界条件，其中三个分量分别表示x,y,z方向，=0 无边界条件施加；=1 位移边界条件；=2 均匀分布荷载；对于复合材料单层板，z方向分量不做处理，boundary_condition_behind 和boundary_condition_front不做处理
boundary_condition_left=[1,1,0];	% 左边界
boundary_condition_right=[1,1,0];	% 右边界
boundary_condition_bottom=[0,0,0]; 	% 下边界
boundary_condition_top=[0,0,0];     % 上边界
boundary_condition_behind=[0,0,0];      % 后边界
boundary_condition_front=[0,0,0];      % 前边界
% 位移边界赋值，其中三个分量分别表示x,y,z方向，如果指定了对应方向的位移边界条件，需在下面给出位移数值；对于复合材料单层板，z方向分量不做处理，displacement_behind和displacement_front不做处理
displacement_left=[-0.5,0,0];      % 左边界
displacement_right=[0.5,0,0];      % 右边界
displacement_bottom=[0,0,0];     % 下边界
displacement_top=[0,0,0];      % 上边界
displacement_behind=[0,0,0];      % 后边界
displacement_front=[0,0,0];      % 前边界
% 均布荷载密度赋值，如果指定了对应方向的荷载条件，需在下面给出荷载密度数值；对于复合材料单层板，z方向分量不做处理，uniform_load_behind和uniform_load_front不做处理
uniform_load_left=[0,0,0];      % 左边界
uniform_load_right=[0,0,0];      % 右边界
uniform_load_bottom=[0,0,0];      % 下边界
uniform_load_top=[0,0,0];      % 上边界
uniform_load_behind=[0,0,0];      % 后边界
uniform_load_front=[0,0,0];      % 前边界
% 体力密度赋初值，如果指定了对应方向的荷载条件，需在下面给体力在各个方向的分量；对于复合材料单层板，z方向分量不做处理
body_force=[0,0,0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


