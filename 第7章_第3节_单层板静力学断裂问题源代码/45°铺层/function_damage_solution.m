function[damage,max_sf,max_sm,damage_particle]=function_damage_solution(particle,size_particle,horizon_particle, ...
    size_horizon_particle,global_displacement_vector,damage_particle,sf,sm,itime,output_interval,size_particle_damage,particle_damage,classify_bond)

damage=zeros(size_particle,1); % 向量damage赋初值为0，向量damage表示每个粒子的损伤程度
max_sf=0;
max_sm=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%二维情况
for ii=1:size_particle_damage                       % 对所有可以发生损伤的粒子循环
	i=particle_damage(ii);                  % 可以发生损伤的粒子集合中的第ii个粒子编号为粒子i
	xx_i=particle(i,1);                     % 获取第i个粒子的x方向坐标
	yy_i=particle(i,2);                     % 获取第i个粒子的y方向坐标
	uu_i=global_displacement_vector(i*2-1,1);   % 获取第i个粒子的x方向的位移
	vv_i=global_displacement_vector(i*2,1);   % 获取第i个粒子的y方向的位移
	number_damage_bond=0;                       % 用于统计第i个粒子域内发生破坏的键的总数，在此赋初值为0
	number_bond=0;                              % 用于统计第i个粒子域内键的总数，在此赋初值为0
	for j=1:size_horizon_particle           % 对第i个粒子域内的粒子进行循环
        jj=horizon_particle(i,j);           % 第i个粒子域内的第j个粒子的粒子编号为jj
        if jj>0                             % 如果第i个粒子域内的第j个粒子存在，则执行条件语句中的内容
            number_bond=number_bond+1;      % 如果第i个粒子域内的第j个粒子存在，第i个粒子域内键的数量加1
            if damage_particle(i,j)>0.5     % 如果第i个粒子域内的第j个键的损伤程度大于0.5，则该键已经被破坏，执行条件语句中的内容
                number_damage_bond=number_damage_bond+1;    % 第i个粒子域内发生破坏的键的数量加1
            end
            if damage_particle(i,j)<0.5         % 如果第i个粒子与第jj个粒子组成的键没有发生破坏，则执行条件语句中的内容
                xx_j=particle(jj,1);            % 获取第jj个粒子的x方向坐标
                yy_j=particle(jj,2);            % 获取第jj个粒子的y方向坐标
                uu_j=global_displacement_vector(jj*2-1,1);   % 获取第jj个粒子的x方向的位移
                vv_j=global_displacement_vector(jj*2,1);   % 获取第jj个粒子的y方向的位移
                length1_i_jj=sqrt((xx_j-xx_i)*(xx_j-xx_i)+(yy_j-yy_i)*(yy_j-yy_i));   % 初始构型下粒子i与粒子jj之间的距离为length_i_jj
                length2_i_jj=sqrt((xx_j+uu_j-xx_i-uu_i)*(xx_j+uu_j-xx_i-uu_i)+(yy_j+vv_j-yy_i-vv_i)*(yy_j+vv_j-yy_i-vv_i));   % 当前构型下粒子i与粒子jj之间的距离为length2_i_jj
                s0=(length2_i_jj-length1_i_jj)/length1_i_jj;   % 粒子i与粒子jj之间键的伸长率为s0
                if classify_bond(i,j)>0.5   % 粒子i与粒子jj组成的键为纤维方向的键
                    if s0>sf                                            % 如果粒子i与粒子jj之间键的伸长率s0大于纤维方向的临界伸长率sf，则执行条件语句中的内容
                        damage_particle(i,j)=1;                         % 粒子i与粒子jj之间键发生破坏，damage_particle(i,j)赋值为1
                        number_damage_bond=number_damage_bond+1;        % 第i个粒子域内发生破坏的键的数量加1
                    end
                    if s0>max_sf    % 如果s0大于max_sf，则执行条件语句中的内容
                        max_sf=s0; % 更新max_sf，保持max_sf取纤维方向键伸长率中的最大值
                    end
                else   % 粒子i与粒子jj组成的键不是纤维方向的键
                    if s0>sm                                            % 如果粒子i与粒子jj之间键的伸长率s0大于临界伸长率sm，则执行条件语句中的内容
                        damage_particle(i,j)=1;                         % 粒子i与粒子jj之间键发生破坏，damage_particle(i,j)赋值为1
                        number_damage_bond=number_damage_bond+1;        % 第i个粒子域内发生破坏的键的数量加1
                    end
                    if s0>max_sm    % 如果s0大于max_sm，则执行条件语句中的内容
                        max_sm=s0; % 更新max_sm，保持max_sm取其他方向键伸长率中的最大值
                    end 
                end
            end
        end
    end
	damage(i,1)=number_damage_bond/number_bond;                 % 第i个粒子的损伤度为其域内发生破坏的键的数量与其域内所有键的数量的比值
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mod(itime,output_interval)==0  % 如果itime对output_interval取余数为0，则执行条件语句中的内容
    str1='damage_';           % 定义输出文件名称中的第1部分 damage_
    str2=num2str(itime);           % 定义输出文件名称中的第2部分 itime
    str3='.txt';                    % 定义输出文件名称中的第3部分 .txt，即文件拓展名
    str4=strcat(str1,str2,str3);    % 将文件名称写到一个变量中，即damage_itime.txt，其中itime指具体计算步数，如果itime=1,则输出文件名为damage_1.txt
    dlmwrite(str4,damage);          % 将damage写入文本文件damage_itime.txt
end