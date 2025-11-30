 function[dimension,dettx,m,length_x,length_y,length_z,thick,size_ply,xx_min,yy_min,zz_min,xx_max,yy_max,zz_max,detla,cf,cm,bn,bs,volume_particle,size_horizon_particle,size_horizon_particle_n,size_horizon_particle_s, ...
    boundary_condition_left,boundary_condition_right,boundary_condition_bottom,boundary_condition_top,boundary_condition_behind,boundary_condition_front,displacement_left, ...
    displacement_right,displacement_bottom,displacement_top,displacement_behind,displacement_front,uniform_load_left,uniform_load_right,uniform_load_bottom, ...
    uniform_load_top,uniform_load_behind,uniform_load_front,body_force,theta,detla_n,detla_s]=function_input()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                             
% 需输入的模型参数
E1=195960;          % 纤维方向弹性模量
E2=8960;            % 垂直纤维方向弹性模量
Em=8960;            % 基体弹性模量
Gm=3500;            % 基体剪切模量
dimension=3;        % 对于复合材料层合板，dimension取3
dettx=1;          % 层内相邻粒子的距离
m=3;                % 粒子作用范围(包含几层粒子)
length_x=50;    % x方向长度
length_y=50;    % y方向长度
thick=0.1;    % 单层度厚
size_ply=3;    % 层合板总层数
theta=[0,pi/4,pi/2];         % 各个铺层纤维方向与x方向的夹角，theta可以取0，pi/4，pi/2
length_z=thick*size_ply;    % 层合板总厚度
xx_min=-25;     % x方向最小坐标
yy_min=-25;     % y方向最小坐标
zz_min=0;     % z方向最小坐标
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xx_max=length_x+xx_min;     % x方向最大坐标
yy_max=length_y+yy_min;     % y方向最大坐标
zz_max=length_z+zz_min;     % z方向最大坐标
detla=(m+0.015)*dettx;  % 粒子层内作用的范围
detla_n=1.015*thick;  % 粒子面外法向作用的范围，只与相邻层粒子相互作用，通常取为厚度thick
detla_s=sqrt(detla*detla+detla_n*detla_n);  % 粒子层间剪切作用的范围，=sqrt(detla*detla+detla_n*detla_n)
cm=8*E1*E2/((E1-E2/9)*pi*thick*detla*detla*detla);	% 基体键微模量系数
volume_particle=dettx*dettx*thick;   %粒子体积
cf_parameter1=2*E1*(E1-E2)/(E1-E2/9);	% 纤维键微模量系数的参数1
size_bond_0_90=[2,3,4,5,6];	% 对于角度theta=0和theta=pi/2，m取不同数值时，模型中包含不同长度的纤维键的数量
size_bond_45=[1,2,2,3,4];	% 对于角度theta=pi/4，m取不同数值时，模型中包含不同长度的纤维键的数量
cf=zeros(1,size_ply);	% 所有铺层纤维键微模量系数赋初值为0
for ii=1:size_ply	% 对所有铺层进行循环
    cf_parameter2=0;	% 纤维键微模量系数的参数2赋初值为0
    if theta(ii)==pi/4  % 如果第ii层纤维方向与x方向的夹角为pi/4，则执行以下内容
        for i=1:size_bond_45(m-1) % m取不同数值时，对模型中包含不同长度的键的数量进行循环
            cf_parameter2=cf_parameter2+2*sqrt(2)*dettx*i*volume_particle; % 对纤维键微模量系数的参数2进行累加，模型包含每种长度键的个数均为2个，因此为2*sqrt(2)*dettx*i*volume_particle
        end
        cf(ii)=cf_parameter1/cf_parameter2;	% 第ii层纤维键微模量系数
    else  % 如果第ii层纤维方向与x方向的夹角为0或pi/2，则执行以下内容
        for i=1:size_bond_0_90(m-1) % m取不同数值时，对模型中包含不同长度的键的数量进行循环
            cf_parameter2=cf_parameter2+2*dettx*i*volume_particle;  % 对纤维键微模量系数的参数2进行累加，模型包含每种长度键的个数均为2个，因此为2*dettx*i*volume_particle
        end
        cf(ii)=cf_parameter1/cf_parameter2;	% 第ii层纤维键微模量系数
    end
end
bn=Em/(4*detla_n*thick*volume_particle);    % 横向键微模量系数
bs=Gm/(16*pi*detla_s*thick*thick*thick*((detla*detla+2*thick*thick)/sqrt(detla*detla+thick*thick)-2*thick));    % 面外剪切键微模量系数

size_horizon_particle_n=2;  % 每个粒子最多包含2个横向键
if m==2
	size_horizon_particle=12;    % 如果m取2，每个粒子在层内最多包含12个键
    size_horizon_particle_s=12;    % 如果m取2，每个粒子与上一层或下一层最多包含12个面外剪切键
end
if m==3
	size_horizon_particle=28;    % 如果m取3，每个粒子在层内最多包含28个键
    size_horizon_particle_s=28;    % 如果m取3，每个粒子与上一层或下一层最多包含28个面外剪切键
end
if m==4
	size_horizon_particle=48;    % 如果m取4，每个粒子在层内最多包含48个键
    size_horizon_particle_s=48;    % 如果m取4，每个粒子与上一层或下一层最多包含48个面外剪切键
end
if m==5
	size_horizon_particle=80;    % 如果m取5，每个粒子在层内最多包含80个键
    size_horizon_particle_s=80;    % 如果m取5，每个粒子与上一层或下一层最多包含80个面外剪切键
end
if m==6
	size_horizon_particle=112;    % 如果m取6，每个粒子在层内最多包含112个键
    size_horizon_particle_s=112;    % 如果m取6，每个粒子与上一层或下一层最多包含112个面外剪切键
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 需输入的边界条件
% 对于复合材料层合板，需要处理与left，right，bottom，top相关的边界条件
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 边界条件，其中三个分量分别表示x,y,z方向，=0 无边界条件施加；=1 位移边界条件；=2 均匀分布荷载；对于复合材料层合板，boundary_condition_behind 和boundary_condition_front不做处理
boundary_condition_left=[1,1,1];	% 左边界
boundary_condition_right=[1,1,1];	% 右边界
boundary_condition_bottom=[0,0,0]; 	% 下边界
boundary_condition_top=[0,0,0];     % 上边界
boundary_condition_behind=[0,0,0];      % 后边界
boundary_condition_front=[0,0,0];      % 前边界
% 位移边界赋值，其中三个分量分别表示x,y,z方向，如果指定了对应方向的位移边界条件，需在下面给出位移数值；对于复合材料层合板，displacement_behind和displacement_front不做处理
displacement_left=[-0.5,0,0];      % 左边界
displacement_right=[0.5,0,0];      % 右边界
displacement_bottom=[0,0,0];     % 下边界
displacement_top=[0,0,0];      % 上边界
displacement_behind=[0,0,0];      % 后边界
displacement_front=[0,0,0];      % 前边界
% 均布荷载密度赋值，如果指定了对应方向的荷载条件，需在下面给出荷载密度数值；对于复合材料层合板，uniform_load_behind和uniform_load_front不做处理
uniform_load_left=[0,0,0];      % 左边界
uniform_load_right=[0,0,0];      % 右边界
uniform_load_bottom=[0,0,0];      % 下边界
uniform_load_top=[0,0,0];      % 上边界
uniform_load_behind=[0,0,0];      % 后边界
uniform_load_front=[0,0,0];      % 前边界
% 体力密度赋初值，如果指定了对应方向的荷载条件，需在下面给体力在各个方向的分量
body_force=[0,0,0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


