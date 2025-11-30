function[damage_particle]=function_initial_crack(particle,size_particle,horizon_particle,size_horizon_particle,dimension,initial_crack, ...
xx1_crack,yy1_crack,zz1_crack,xx2_crack,yy2_crack,zz2_crack,dettx)

damage_particle=zeros(size_particle,size_horizon_particle);    % damage_particle矩阵清零
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 如果存在初始裂纹，则需找出被初始裂纹切断的键，并将键的损伤变量设置为1
if initial_crack==1    % initial_crack=1表示模型中存在初始裂纹
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 二维情况
    if dimension==2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 竖直裂纹
        if abs(xx1_crack-xx2_crack)<0.01*dettx    % 如果裂纹起始点的x坐标与裂纹结束点的x坐标差值的绝对值小于0.01*dettx，则认定该裂纹为竖直裂纹
            for i=1:size_particle %对粒子循环
                xx_i=particle(i,1); % 获取第i个粒子的x方向坐标
                yy_i=particle(i,2); % 获取第i个粒子的y方向坐标
                for j=1:size_horizon_particle           % 对第i个粒子域内的粒子进行循环
                    jj=horizon_particle(i,j);           % 第i个粒子域内的第j个粒子的粒子编号为jj
                    if jj>0                             % 如果第i个粒子域内的第j个粒子存在，则执行条件语句中的内容
                        xx_j=particle(jj,1);            % 获取第jj个粒子的x方向坐标
                        yy_j=particle(jj,2);            % 获取第jj个粒子的y方向坐标
                        if xx_i>xx1_crack && xx_j<xx1_crack     % 如果第i个粒子在裂纹线右侧，第j个粒子在裂纹线左侧
                            detl_y1=abs(yy_j-yy_i);             % detl_y1表示第i个粒子和第jj个粒子在y方向的坐标差
                            if yy_j>yy_i                        % 如果第jj个粒子的y坐标大于第i个粒子的y坐标
                                detl_y2=detl_y1*abs(xx1_crack-xx_i)/abs(xx_j-xx_i); %按case1计算长度detl_y2
                                yy_crack=yy_i+detl_y2;      %由第i个粒子和第jj个粒子组成的键与裂纹线所在直线的交点的y坐标
                            end
                            if yy_j<=yy_i                    % 如果第jj个粒子的y坐标小于等于第i个粒子的y坐标
                                detl_y2=detl_y1*abs(xx1_crack-xx_j)/abs(xx_j-xx_i); %按case2计算长度detl_y2
                                yy_crack=yy_j+detl_y2;      %由第i个粒子和第jj个粒子组成的键与裂纹线所在直线的交点的y坐标
                            end
                            if yy_crack>=yy1_crack && yy_crack<=yy2_crack % 如果键与裂纹线所在直线的交点的y坐标位于裂纹起始点的y坐标与裂纹结束点的y坐标之间，则该键被裂纹线切断
                                damage_particle(i,j)=1;    % damage_particle(i,jj)赋值为1代表由第i个粒子和第jj个粒子组成的键被裂纹线切断
                            end
                        end
                        if xx_i<xx1_crack && xx_j>xx1_crack  % 如果第i个粒子在裂纹线左侧，第j个粒子在裂纹线右侧
                            detl_y1=abs(yy_j-yy_i);  % detl_x1表示第i个粒子和第jj个粒子在y方向的坐标差
                            if yy_j>yy_i                    % 如果第jj个粒子的y坐标大于第i个粒子的y坐标
                                detl_y2=detl_y1*abs(xx1_crack-xx_i)/abs(xx_j-xx_i); %按case3计算长度detl_y2
                                yy_crack=yy_i+detl_y2;      %由第i个粒子和第jj个粒子组成的键与裂纹线所在直线的交点的y坐标
                            end
                            if yy_j<=yy_i                    % 如果第jj个粒子的yx坐标小于等于第i个粒子的y坐标
                                detl_y2=detl_y1*abs(xx1_crack-xx_j)/abs(xx_j-xx_i); %按case2计算长度detl_y2
                                yy_crack=yy_j+detl_y2;      %由第i个粒子和第jj个粒子组成的键与裂纹线所在直线的交点的y坐标
                            end
                            if yy_crack>=yy1_crack && yy_crack<=yy2_crack % 如果键与裂纹线所在直线的交点的y坐标位于裂纹起始点的y坐标与裂纹结束点的y坐标之间，则该键被裂纹线切断
                                damage_particle(i,j)=1;    % damage_particle(i,jj)赋值为1代表由第i个粒子和第jj个粒子组成的键被裂纹线切断
                            end
                        end
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 斜裂纹
        if abs(xx1_crack-xx2_crack)>0.01*dettx    % 如果裂纹起始点的x坐标与裂纹结束点的x坐标差值的绝对值大于0.01*dettx，则认定该裂纹不是竖直裂纹
            kk_crack=(yy2_crack-yy1_crack)/(xx2_crack-xx1_crack);    % 由裂纹起始点与裂纹结束点组成的直线方程 y=kx+b 中的参数 k
            bb_crack=yy2_crack-xx2_crack*kk_crack;                   % 由裂纹起始点与裂纹结束点组成的直线方程 y=kx+b 中的参数 b
            for i=1:size_particle                                   %对粒子循环
                xx_i=particle(i,1);                                 % 获取第i个粒子的x方向坐标
                yy_i=particle(i,2);                                 % 获取第i个粒子的y方向坐标
                for j=1:size_horizon_particle                       % 对第i个粒子域内的粒子进行循环
                    jj=horizon_particle(i,j);                       % 第i个粒子域内的第j个粒子的粒子编号为jj
                    if jj>0                                         % 如果第i个粒子域内的第j个粒子存在，则执行条件语句中的内容
                        xx_j=particle(jj,1);                        % 获取第jj个粒子的x方向坐标
                        yy_j=particle(jj,2);                        % 获取第jj个粒子的y方向坐标
                        if abs(xx_i-xx_j)<0.01*dettx                                      % 如果第i个粒子的x坐标第jj个粒子的x坐标相等，则执行条件语句中的内容
                            if xx_i>xx1_crack && xx_i<xx2_crack             % 如果第i个粒子的x坐标位于由裂纹起始点的x坐标与裂纹结束点的x坐标之间，则执行条件语句中的内容
                                yy_i_crack=xx_i*kk_crack+bb_crack;          % 由第i个粒子与第jj个粒子组成的键所在直线与裂纹的交点
                                det_yy_i_crack=yy_i-yy_i_crack;             % 由第i个粒子的y坐标与上述交点的y坐标的差值，差值为正表示粒子i位于裂纹线上侧，差值为负表示粒子i位于裂纹线下侧
                                det_yy_j_crack=yy_j-yy_i_crack;             % 由第jj个粒子的y坐标与上述交点的y坐标的差值，差值为正表示粒子jj位于裂纹线上侧，差值为负表示粒子jj位于裂纹线下侧
                                if det_yy_i_crack*det_yy_j_crack<0          % 如果第i个粒子的y坐标与上述交点的y坐标的差值以及第jj个粒子的y坐标与上述交点的y坐标的差值的乘积为负，表示粒子i与粒子jj位于裂纹的两侧，则执行条件语句中的内容
                                    damage_particle(i,j)=1;                 % damage_particle(i,jj)赋值为1代表由第i个粒子和第jj个粒子组成的键被裂纹线切断
                                end
                            end
                        else                                                            % 如果第i个粒子的x坐标不位于由裂纹起始点的x坐标与裂纹结束点的x坐标之间，则执行else语句中的内容
                            kk_bond=(yy_j-yy_i)/(xx_j-xx_i);                            % 由第i个粒子与第jj个粒子组成的直线方程 y=kx+b 中的参数 k
                            bb_bond=yy_i-xx_i*kk_bond;                                  % 由第i个粒子与第jj个粒子组成的直线方程 y=kx+b 中的参数 b
                            if kk_crack~=kk_bond                                        % 如果经过裂纹起始点与裂纹结束点的直线斜率与经过第i个粒子与第jj个粒子的直线斜率不相等，则两条直线一定存在交点，则执行条件语句中的内容
                                xx_point=(bb_bond-bb_crack)/(kk_crack-kk_bond);         % xx_point表示以上两条直线交点的x坐标            
                                if xx_point<xx2_crack && xx_point>xx1_crack             % 如果xx_point的x坐标位于由裂纹起始点的x坐标与裂纹结束点的x坐标之间，则执行条件语句中的内容
                                    yy_i_crack=xx_i*kk_crack+bb_crack;                  % 经过第i个粒子的y方向的直线与裂纹线的交点
                                    det_yy_i_crack=yy_i-yy_i_crack;                     % 由第i个粒子的y坐标与上述交点的y坐标的差值，差值为正表示粒子i位于裂纹线上侧，差值为负表示粒子i位于裂纹线下侧
                                    yy_j_crack=xx_j*kk_crack+bb_crack;                  % 经过第jj个粒子的y方向的直线与裂纹线的交点
                                    det_yy_j_crack=yy_j-yy_j_crack;                     % 由第jj个粒子的y坐标与上述交点的y坐标的差值，差值为正表示粒子jj位于裂纹线上侧，差值为负表示粒子jj位于裂纹线下侧
                                    if det_yy_i_crack*det_yy_j_crack<0                  % 如果第i个粒子的y坐标与上述交点的y坐标的差值以及第jj个粒子的y坐标与上述交点的y坐标的差值的乘积为负，表示粒子i与粒子jj位于裂纹的两侧，则执行条件语句中的内容
                                        damage_particle(i,j)=1;                         % damage_particle(i,jj)赋值为1代表由第i个粒子和第jj个粒子组成的键被裂纹线切断
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
