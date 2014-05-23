clear;
%% 程序作者：
% 孙晓博（2010011732）
% 工物系核01班
% sunxb10@gmail.com
% Nov.16th, 2013
%% 输入参数：
N = 8;  % S_N近似的阶数
K = 50*N; % 沿r方向的网格数
tol = 1e-8;  % 判断S_N近似收敛的阈值
%% 计算初值：
sigma_t = 0.050; % 宏观截面数据(cm^(-1))
sigma_s = 0.030;
nusigma_f = 0.0225;
a0 = 66.0053;  % 平板半厚度(cm)
delta = a0/K;  % 差分步长
x = [0:delta:a0];
[mu,wt] = getGaussianQuadTuple(N);  % 获取求积组，节点mu与权重wt的元素成对出现，其中mu的元素按照“负正负正负正”的顺序交错排列，偶数位上的元素皆为正
l = zeros(N,1);
keff = zeros(2*K+1,1);  % 有效增值系数
Phi = zeros(2*K+1,1);  % 中子注量率
old_Phi = 100*ones(2*K+1,1);  % 迭代前的中子注量率初值设为10
src = 1e12*ones(2*K+1,N); % 源项初值设为1e10
Psi = zeros(2*K+1,N);  % 中子角注量率
bad_grid = 1;  % bad_grid是迭代前后相差超过tol的节点，初始置为1
cycle_count = 0; % 统计循环迭代次数

% f和f_mu用于处理各向异性散射问题（需要输入介质的等效原子质量MA）
% % 散射函数的勒让德展开
% f_mu =  @(mu) (sqrt(mu.^2+MA^2-1)+mu).^2./(2*MA*sqrt(mu.^2+MA^2-1));
% f = zeros(N,1,'double');
% for i= 1:N   
%     f(i) = myRomberg(@(mu) f_mu(mu).*myLegendre((i-1),mu),-1,1);
% end
% MatLegendre = getLegendreList(N,mu);  % 各离散角度的勒让德函数值，先提前算好并储存于MatLegendre中，如此可避免重复计算

%% 迭代计算：
t_Start = tic; % MATLAB计时语句
while(bad_grid > 0)
    % mu<0的部分
    for i = 1:(N/2)
        l(i) = delta*sigma_t/mu(i);
        Psi(2*K+1,i) = 0;  % 右侧真空边界
        for k = K:(-1):1  % 向前迭代
           Psi(2*k-1,i) = (2+l(i))/(2-l(i))*Psi(2*k+1,i)-2*l(i)/(2-l(i))*src(2*k)/sigma_t;
           Psi(2*k,i) = (Psi(2*k+1,i)+Psi(2*k-1,i))/2;
        end
    end
    % mu>0的部分
    for i = (N/2+1):N
        l(i) = delta*sigma_t/mu(i);
        Psi(1,i) = Psi(1,N+1-i);  % 左侧对称条件
        for k = 1:K  % 向后迭代
           Psi(2*k+1,i) = (2-l(i))/(2+l(i))*Psi(2*k-1,i)+2*l(i)/(2+l(i))*src(2*k)/sigma_t;
           Psi(2*k,i) = (Psi(2*k+1,i)+Psi(2*k-1,i))/2;
        end
    end
    % 更新源项
    for i = 1:K
        src(2*i) = (sigma_s+nusigma_f)/4*wt*(Psi(2*i-1,:)+Psi(2*i+1,:))';
%         % 使用各向异性散射函数更新源项
%         for p = 1:N
%             src(2*i,p) = 0;
%             for q = 1:N
%                 src(2*i,p) = src(2*i,p) + (2*q-1)/4*sigma_s*f(q)*MatLegendre(p,q)*(wt*((psi(2*i-1,:)+psi(2*i+1,:))'.*MatLegendre(:,q)));
%             end
%             src(2*i,p) = src(2*i,p) + nusigma_f/4*wt*(psi(2*i-1,:)+psi(2*i+1,:))';
%         end
    end
    % 检查是否达到收敛要求
    bad_grid = 0;
    for i = 1:(2*K+1)
        Phi(i) = wt*Psi(i,:)';  % 根据角注量率计算注量率（高斯积分）
        keff_temp = Phi(i)/old_Phi(i);  % 计算当前格点的有效增值系数
        if(abs((keff_temp-keff(i))/keff(i))>=tol)  % 根据有效增值系数来判断是否达到要求
           bad_grid = bad_grid + 1;
        end
        keff(i) = keff_temp;
    end
    old_Phi = Phi;
    cycle_count = cycle_count + 1;
end
elapsed_time = toc(t_Start); % MATLAB计时语句

%% 输出计算结果：
fp = fopen('1D_plane_outputs.txt','wt');
fprintf(fp,'有效增值系数（平均值）： %.8g\n',mean(keff));
fprintf(fp,'有效增值系数（最大值）： %.8g\n',max(keff));
fprintf(fp,'有效增值系数（最小值）： %.8g\n',min(keff));
fprintf(fp,'循环迭代总次数： %d\n',cycle_count);
fprintf(fp,'循环迭代总时长： %.3g s\n\n\n',elapsed_time);
fprintf(fp,'************************************\n');
fprintf(fp,'*           中子通量校验           *\n');
fprintf(fp,'************************************\n');
phi_center = Phi(1); % 中心通量
fprintf(fp,'   校验点位置           模拟结果\n');
fprintf(fp,'（a为平板半厚度）   （中心通量归一）\n');
fprintf(fp,'------------------------------------\n');
fprintf(fp,'     0.25a             %.9g\n',Phi(uint32(0.25*a0/(delta/2)))/phi_center); % 对中心通量归一
fprintf(fp,'     0.50a             %.9g\n',Phi(uint32(0.50*a0/(delta/2)))/phi_center);
fprintf(fp,'     0.75a             %.9g\n',Phi(uint32(0.75*a0/(delta/2)))/phi_center);
fprintf(fp,'     1.00a             %.9g\n',Phi(uint32(1.00*a0/(delta/2)))/phi_center);
fclose(fp);

%% 绘制中子通量分布曲线：
plot(x,Phi(1:2:(2*K+1)),'-r','linewidth',1.5);
xlabel('$x$','FontSize',16,'Interpreter','latex','FontWeight','bold');
ylabel('Fluence Rate $\phi$','FontSize',16,'Interpreter','latex','FontWeight','bold');
print('-depsc2','1D_plane_neutron_flux.eps'); % 将中子通量分布曲线存入PDF文件保存