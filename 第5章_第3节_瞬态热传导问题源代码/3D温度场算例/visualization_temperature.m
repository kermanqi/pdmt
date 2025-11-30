size_x=10;
size_y=10;
size_z=26;
global_temperature_vector=load('temperature_30000.txt');
temperature=zeros(size_x,size_y,size_z);    % 温度写为矩阵形式
ii=0;                                   % 变量初始值取0
for i=1:size_z                          % 按z方向粒子总数循环
	for j=1:size_y                      % 按y方向粒子总数循环
        for k=1:size_x                  % 按x方向粒子总数循环
            ii=ii+1;                        % ii取值+1，从而遍历每个粒子
            temperature(k,j,i)=global_temperature_vector(ii,1);   % 总体温度矢量写入温度矩阵
        end
	end
end
figure(1)
slice(temperature,[1,size_y],[1,size_x],[1,size_z]);   % 温度结果可视化
set(gca,'DataAspectRatio',[1 1 1]);  % 设置坐标轴比例
set(colorbar,'FontSize',14,'FontName','Times New Roman');  % 设置图例的字号字体
axis off
colormap Jet  % 图像颜色映射设置为Jet