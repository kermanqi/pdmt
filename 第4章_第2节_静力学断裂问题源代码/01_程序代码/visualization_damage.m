size_x=106;
size_y=100;
damage=load('damage_240.txt');
damage_matrix=zeros(size_y,size_x);    % 粒子损伤度写为矩阵形式
ii=0;                                   % 变量初始值取0
for i=1:size_y                          % 按y方向粒子总数循环
	for j=1:size_x                      % 按x方向粒子总数循环
	ii=ii+1;                        % ii取值+1，从而遍历每个粒子
	damage_matrix(i,j)=damage(ii,1);     % 粒子损伤度写入矩阵damage_matrix
	end
end
pcolor(damage_matrix);shading interp;  % 粒子损伤度结果可视化
set(colorbar,'FontSize',14,'FontName','Times New Roman');  % 设置图例的字号字体
colormap Jet  % 图像颜色映射设置为Jet