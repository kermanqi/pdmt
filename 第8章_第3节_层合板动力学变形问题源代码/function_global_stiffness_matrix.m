function[global_stiffness_matrix]=function_global_stiffness_matrix(particle,size_particle,horizon_particle,size_horizon_particle,horizon_particle_n,size_horizon_particle_n, ...
    horizon_particle_s_lower,horizon_particle_s_upper,cf,cm,bn,bs,volume_particle,classify_bond,detla_n,detla_s,particle_ply)

size_vector=0; % size_vector赋初值为0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 面内作用
for i=1:size_particle                       % 对所有粒子进行循环
	xx_i=particle(i,1);                     % 获取第i个粒子的x坐标
	yy_i=particle(i,2);                     % 获取第i个粒子的y坐标
	for j=1:size_horizon_particle           % 对所有与粒子i进行面内作用的粒子进行循环
        jj=horizon_particle(i,j);           % 与粒子i进行面内作用的第j个粒子的粒子编号为jj
        if jj>0                             % 如果第i个粒子域内的第j个粒子存在，则执行条件语句中的内容
            xx_j=particle(jj,1);            % 获取第jj个粒子的x方向坐标
            yy_j=particle(jj,2);            % 获取第jj个粒子的y方向坐标
            length_i_jj=sqrt((xx_j-xx_i)*(xx_j-xx_i)+(yy_j-yy_i)*(yy_j-yy_i));   % 粒子i与粒子jj之间的距离为length_i_jj
            detl_x=(xx_j-xx_i)/length_i_jj;                                     % 计算式子()中的detx
            detl_y=(yy_j-yy_i)/length_i_jj;                                     % 计算式子()中的dety
            if classify_bond(i,j)==1	% classify_bond(i,j)等于1表示纤维方向的键
                ii=particle_ply(i);     % 获取粒子i所在层数
                c0=cf(ii)+cm;   % 微模量系数取纤维键微模量系数与基体键微模量系数的和
            else	% classify_bond(i,j)不等于1表示其他方向的键
                c0=cm;   % 微模量系数取基体键微模量系数
            end
            kk_PD=c0/length_i_jj*[detl_x*detl_x detl_x*detl_y -detl_x*detl_x -detl_x*detl_y; detl_y*detl_x detl_y*detl_y -detl_y*detl_x -detl_y*detl_y]*volume_particle*volume_particle; % 面内作用的刚度矩阵
            i1=[i*3-2,i*3-1];                 % 需要组装到总体刚度矩阵中的行向量
            j1=[i*3-2,i*3-1,jj*3-2,jj*3-1];     % 需要组装到总体刚度矩阵中的列向量
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 面外法向作用
for i=1:size_particle       % 对所有粒子进行循环
    zz_i=particle(i,3);     % 获取第i个粒子的z坐标
	for j=1:size_horizon_particle_n           % 对所有与粒子i进行面外法向作用的粒子进行循环
        jj=horizon_particle_n(i,j);           % 与粒子i进行面外法向作用的第j个粒子的粒子编号为jj
        if jj>0                             % 如果第i个粒子域内的第j个粒子存在，则执行条件语句中的内容
            zz_j=particle(jj,3);    % 获取第jj个粒子的z坐标
            length_i_jj=abs(zz_j-zz_i);    % 粒子i与粒子jj之间的距离为length_i_jj
            detl_z=(zz_j-zz_i)/length_i_jj;     % 计算式子()中的detz                                
            kk_PD_n=4*bn*detla_n/length_i_jj*[detl_z*detl_z -detl_z*detl_z]*volume_particle*volume_particle; % 面外法向作用的刚度矩阵 
            i1=i*3;                 % 需要组装到总体刚度矩阵中的行向量
            j1=[i*3,jj*3];     % 需要组装到总体刚度矩阵中的列向量
            for ii3=1:1                     % 对行向量进行循环
                for jj3=1:2                 % 对列向量进行循环
                    size_vector=size_vector+1;  % 变量+1
                    vector_x(size_vector)=i1(ii3);  % 行向量元素赋值
                    vector_y(size_vector)=j1(jj3);  % 列向量元素赋值
                    stiffness_matrix(size_vector)=kk_PD_n(ii3,jj3);  % 总体刚度矩阵元素赋值
                end
            end
        end
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 层间剪切作用
for i=1:size_particle       % 对所有粒子进行循环
	xx_i=particle(i,1);     % 获取第i个粒子的x坐标
	yy_i=particle(i,2);     % 获取第i个粒子的y坐标
    zz_i=particle(i,3);     % 获取第i个粒子的z坐标
    for ii=1:size_horizon_particle_n    % 对所有与粒子i进行面外法向作用的粒子进行循环
        iin=horizon_particle_n(i,ii);   % 与粒子i进行面外法向作用的第ii个粒子的粒子编号为iin
        if iin>0                        % 如果第i个粒子域内的第ii个粒子存在，则执行条件语句中的内容
            xx_iin=particle(iin,1);  % 获取第iin个粒子的x坐标                   
            yy_iin=particle(iin,2);  % 获取第iin个粒子的y坐标                   
            zz_iin=particle(iin,3);  % 获取第iin个粒子的z坐标
            for j=1:size_horizon_particle           % 对所有与粒子i进行面内作用的粒子进行循环
                jj=horizon_particle(i,j);           % 与粒子i进行面内作用的第j个粒子的粒子编号为jj
                if jj>0                             % 如果第i个粒子域内的第j个粒子存在，则执行条件语句中的内容
                    xx_j=particle(jj,1);            % 获取第jj个粒子的x坐标
                    yy_j=particle(jj,2);            % 获取第jj个粒子的y坐标
                    zz_j=particle(jj,3);            % 获取第jj个粒子的z坐标
                    if ii==1                        % 如果ii等于1，则表示粒子iin为粒子i下面一层的粒子
                        jjs=horizon_particle_s_lower(i,j); % 获取位于粒子jj下面一层的粒子jjs
                    end
                    if ii==2                        % 如果ii等于2，则表示粒子iin为粒子i上面一层的粒子
                        jjs=horizon_particle_s_upper(i,j); % 获取位于粒子jj上面一层的粒子jjs
                    end
                    xx_jjs=particle(jjs,1);     % 获取第jjs个粒子的x坐标                
                    yy_jjs=particle(jjs,2);     % 获取第jjs个粒子的y坐标               
                    zz_jjs=particle(jjs,3);     % 获取第jjs个粒子的z坐标
                    length_i_jjs=sqrt((xx_jjs-xx_i)*(xx_jjs-xx_i)+(yy_jjs-yy_i)*(yy_jjs-yy_i)+(zz_jjs-zz_i)*(zz_jjs-zz_i)); % 粒子i与粒子jjs之间的距离为length_i_jjs
                    detl_x_i_jjs=(xx_jjs-xx_i)/length_i_jjs;   % 计算式子()中的detx_i_jjs
                    detl_y_i_jjs=(yy_jjs-yy_i)/length_i_jjs;   % 计算式子()中的dety_i_jjs 
                    detl_z_i_jjs=(zz_jjs-zz_i)/length_i_jjs;   % 计算式子()中的detz_i_jjs 
                    length_j_iin=sqrt((xx_iin-xx_j)*(xx_iin-xx_j)+(yy_iin-yy_j)*(yy_iin-yy_j)+(zz_iin-zz_j)*(zz_iin-zz_j));  % 粒子j与粒子iin之间的距离为length_j_iin 
                    detl_x_j_iin=(xx_iin-xx_j)/length_j_iin;   % 计算式子()中的detx_j_iin 
                    detl_y_j_iin=(yy_iin-yy_j)/length_j_iin;   % 计算式子()中的dety_j_iin 
                    detl_z_j_iin=(zz_iin-zz_j)/length_j_iin;   % 计算式子()中的detz_j_iin 
                    kk_PD_s1=8*bs*detla_s/length_i_jjs*[detl_x_i_jjs*detl_x_i_jjs detl_x_i_jjs*detl_y_i_jjs detl_x_i_jjs*detl_z_i_jjs -detl_x_i_jjs*detl_x_i_jjs -detl_x_i_jjs*detl_y_i_jjs -detl_x_i_jjs*detl_z_i_jjs;
                                                        detl_y_i_jjs*detl_x_i_jjs detl_y_i_jjs*detl_y_i_jjs detl_y_i_jjs*detl_z_i_jjs -detl_y_i_jjs*detl_x_i_jjs -detl_y_i_jjs*detl_y_i_jjs -detl_y_i_jjs*detl_z_i_jjs;
                                                        detl_z_i_jjs*detl_x_i_jjs detl_z_i_jjs*detl_y_i_jjs detl_z_i_jjs*detl_z_i_jjs -detl_z_i_jjs*detl_x_i_jjs -detl_z_i_jjs*detl_y_i_jjs -detl_z_i_jjs*detl_z_i_jjs]*volume_particle*volume_particle;% 层间剪切作用的刚度矩阵1 
                    i1=[i*3-2,i*3-1,i*3];                 % 需要组装到总体刚度矩阵中的行向量
                    j1=[i*3-2,i*3-1,i*3,jjs*3-2,jjs*3-1,jjs*3];     % 需要组装到总体刚度矩阵中的列向量
                    for ii3=1:3                     % 对行向量进行循环
                        for jj3=1:6                 % 对列向量进行循环
                            size_vector=size_vector+1;  % 变量+1
                            vector_x(size_vector)=i1(ii3);  % 行向量元素赋值
                            vector_y(size_vector)=j1(jj3);  % 列向量元素赋值
                            stiffness_matrix(size_vector)=kk_PD_s1(ii3,jj3);  % 总体刚度矩阵元素赋值
                        end
                    end
                    kk_PD_s2=8*bs*detla_s/length_j_iin*[detl_x_i_jjs*detl_x_j_iin detl_x_i_jjs*detl_y_j_iin detl_x_i_jjs*detl_z_j_iin -detl_x_i_jjs*detl_x_j_iin -detl_x_i_jjs*detl_y_j_iin -detl_x_i_jjs*detl_z_j_iin;
                                                        detl_y_i_jjs*detl_x_j_iin detl_y_i_jjs*detl_y_j_iin detl_y_i_jjs*detl_z_j_iin -detl_y_i_jjs*detl_x_j_iin -detl_y_i_jjs*detl_y_j_iin -detl_y_i_jjs*detl_z_j_iin;
                                                        detl_z_i_jjs*detl_x_j_iin detl_z_i_jjs*detl_y_j_iin detl_z_i_jjs*detl_z_j_iin -detl_z_i_jjs*detl_x_j_iin -detl_z_i_jjs*detl_y_j_iin -detl_z_i_jjs*detl_z_j_iin]*volume_particle*volume_particle;% 层间剪切作用的刚度矩阵2
                    i1=[i*3-2,i*3-1,i*3];                 % 需要组装到总体刚度矩阵中的行向量
                    j1=[iin*3-2,iin*3-1,iin*3,jj*3-2,jj*3-1,jj*3];     % 需要组装到总体刚度矩阵中的列向量
                    for ii3=1:3                     % 对行向量进行循环
                        for jj3=1:6                 % 对列向量进行循环
                            size_vector=size_vector+1;  % 变量+1
                            vector_x(size_vector)=i1(ii3);  % 行向量元素赋值
                            vector_y(size_vector)=j1(jj3);  % 列向量元素赋值
                            stiffness_matrix(size_vector)=kk_PD_s2(ii3,jj3);  % 总体刚度矩阵元素赋值
                        end
                    end
                end
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global_stiffness_matrix=sparse(vector_x,vector_y,stiffness_matrix,size_particle*3,size_particle*3); % 合成总体刚度矩阵