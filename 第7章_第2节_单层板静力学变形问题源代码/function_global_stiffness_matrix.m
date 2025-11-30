function[global_stiffness_matrix]=function_global_stiffness_matrix(particle,size_particle,horizon_particle,size_horizon_particle,cf,cm,volume_particle,classify_bond)

size_vector=0; % size_vector赋初值为0
for i=1:size_particle                       % 对粒子循环
	xx_i=particle(i,1);                     % 获取第i个粒子的x方向坐标
	yy_i=particle(i,2);                     % 获取第i个粒子的y方向坐标
	for j=1:size_horizon_particle           % 对第i个粒子域内的粒子进行循环
        jj=horizon_particle(i,j);           % 第i个粒子域内的第j个粒子的粒子编号为jj
        if jj>0                             % 如果第i个粒子域内的第j个粒子存在，则执行条件语句中的内容
            xx_j=particle(jj,1);            % 获取第jj个粒子的x方向坐标
            yy_j=particle(jj,2);            % 获取第jj个粒子的y方向坐标
            length_i_jj=sqrt((xx_j-xx_i)*(xx_j-xx_i)+(yy_j-yy_i)*(yy_j-yy_i));   % 粒子i与粒子jj之间的距离为length_i_jj
            detl_x=(xx_j-xx_i)/length_i_jj;                                     % 计算式子()中的detx
            detl_y=(yy_j-yy_i)/length_i_jj;                                     % 计算式子()中的dety
            if classify_bond(i,j)==1	% classify_bond(i,j)等于1表示纤维方向的键
                c0=cf+cm;   % 微模量系数取纤维键微模量系数与基体键微模量系数的和
            else	% classify_bond(i,j)不等于1表示其他方向的键
                c0=cm;   % 微模量系数取基体键微模量系数
            end
            kk_PD=c0/length_i_jj*[detl_x*detl_x detl_x*detl_y -detl_x*detl_x -detl_x*detl_y; detl_y*detl_x detl_y*detl_y -detl_y*detl_x -detl_y*detl_y]*volume_particle*volume_particle; % 粒子i与粒子jj组成键的刚度矩阵
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
global_stiffness_matrix=sparse(vector_x,vector_y,stiffness_matrix,size_particle*2,size_particle*2); % 合成总体刚度矩阵