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
R = 145.5436;  % 球的半径(cm)
delta = R/K;  % 差分步长
mu = zeros(2*N+1,1);
wt = zeros(2*N+1,1);
[mu(2:2:2*N),wt(2:2:2*N)] = getGaussianQuadTuple(N);  % 计算mu_{m}与wt_{2m}
mu(1:2*N:2*N+1) = [-1,1];
mu(3:2:2*N-1) = 0.5*(mu(2:2:2*N-2)+mu(4:2:2*N));
r = [0:delta/2:R]; % 计算r_{k}
A = zeros(2*K+1,1);
A(1:2:2*K+1)=r(1:2:2*K+1).^2; % 计算A_{2k+1}
V = zeros(2*K,1);
V(2:2:2*K)=(1/3)*(r(3:2:2*K+1).^3-r(1:2:2*K-1).^3); % 计算V_{2k}
alpha = zeros(2*N+1,1);
for i = 1:N
   alpha(2*i+1) = alpha(2*i-1) - mu(2*i)*wt(2*i);  % 计算alpha_{2m+1}
end
E = zeros(2*K,2*N);
E(2:2:2*K,2:2:2*N) = (A(1:2:2*K-1)+A(3:2:2*K+1))*abs(mu(2:2:2*N))'; % 计算E_{2k,2m}
G = zeros(2*K,2*N);
G(2:2:2*K,2:2:2*N) = (A(3:2:2*K+1)-A(1:2:2*K-1))*((alpha(3:2:2*N+1)+alpha(1:2:2*N-1))./wt(2:2:2*N))'; % 计算G_{2k,2m}
bad_grid = 1;  % bad_grid是迭代前后相差超过tol的节点，初始置为1
cycle_count = 0; % 统计循环迭代次数
Psi = zeros(2*K+1,2*N+1); % 中子角通量初值
src = 1e13*ones(2*K,1); % 源项
keff = ones(2*K,1); % 有效增殖系数
Phi = zeros(2*K,1); % 中子通量
old_Phi = 1e13*ones(2*K,1);

%% 迭代计算：
t_Start = tic; % MATLAB计时语句
% 迭代的具体公式请参阅本人所写的说明文档：单速一维均匀裸球堆S_N模拟的MATLAB实现.pdf
while(bad_grid > 0)
    % mu = -1
    for k = K:(-1):1
        Psi(2*k,1) = (V(2*k)*src(2*k)+(A(2*k+1)+A(2*k-1))*Psi(2*k+1,1))/(A(2*k+1)+A(2*k-1)+sigma_t*V(2*k));
        Psi(2*k-1,1) = 2*Psi(2*k,1) - Psi(2*k+1,1);
    end
    % mu < 0
    for m = 1:(N/2)
        for k = K:(-1):1
            Psi(2*k,2*m) = (V(2*k)*src(2*k)+E(2*k,2*m)*Psi(2*k+1,2*m)+G(2*k,2*m)*Psi(2*k,2*m-1))/(E(2*k,2*m)+G(2*k,2*m)+sigma_t*V(2*k));
            if(Psi(2*k,2*m) < 0)
                Psi(2*k,2*m) = 0;
            end
            Psi(2*k-1,2*m) = 2*Psi(2*k,2*m) - Psi(2*k+1,2*m);
            if(Psi(2*k-1,2*m) < 0)
                Psi(2*k-1,2*m) = 0;
            end
            Psi(2*k,2*m+1) = 2*Psi(2*k,2*m) - Psi(2*k,2*m-1);
            if(Psi(2*k,2*m+1) < 0)
                Psi(2*k,2*m+1) = 0;
            end
        end
    end
    % mu > 0
    for m = (N/2+1):N
        Psi(1,2*m) = Psi(1,2*N+2-2*m); % 球心对称条件
        for k = 1:K
            Psi(2*k,2*m) = (V(2*k)*src(2*k)+E(2*k,2*m)*Psi(2*k-1,2*m)+G(2*k,2*m)*Psi(2*k,2*m-1))/(E(2*k,2*m)+G(2*k,2*m)+sigma_t*V(2*k));
            if(Psi(2*k,2*m) < 0)
                Psi(2*k,2*m) = 0;
            end
            Psi(2*k+1,2*m) = 2*Psi(2*k,2*m) - Psi(2*k-1,2*m);
            if(Psi(2*k+1,2*m) < 0)
                Psi(2*k+1,2*m) = 0;
            end
            Psi(2*k,2*m+1) = 2*Psi(2*k,2*m) - Psi(2*k,2*m-1);
            if(Psi(2*k,2*m+1) < 0)
                Psi(2*k,2*m+1) = 0;
            end
        end
    end
    % 更新源项
    src(2:2:2*K) = (sigma_s + nusigma_f)/2*Psi(2:2:2*K,2:2:2*N)*wt(2:2:2*N);
    % 检查是否达到收敛要求
    bad_grid = 0;
    for k = 1:K
        Phi(2*k) = Psi(2*k,2:2:2*N)*wt(2:2:2*N);  % 根据角注量率计算注量率（高斯积分）
        keff_temp = Phi(2*k)/old_Phi(2*k);  % 计算当前格点的有效增值系数
        if(abs((keff_temp-keff(2*k))/keff(2*k))>=tol)  % 根据有效增值系数来判断是否达到要求
           bad_grid = bad_grid + 1;
        end
        keff(2*k) = keff_temp;
    end
    old_Phi = Phi;
    cycle_count = cycle_count + 1;
end
elapsed_time = toc(t_Start); % MATLAB计时语句

%% 输出计算结果：
fp = fopen('1D_sphere_outputs.txt','wt');
fprintf(fp,'有效增值系数（平均值）： %.8g\n',mean(keff(2:2:2*K)));
fprintf(fp,'有效增值系数（最大值）： %.8g\n',max(keff(2:2:2*K)));
fprintf(fp,'有效增值系数（最小值）： %.8g\n',min(keff(2:2:2*K)));
fprintf(fp,'循环迭代总次数： %d\n',cycle_count);
fprintf(fp,'循环迭代总时长： %.3g s\n\n\n',elapsed_time);
fprintf(fp,'************************************\n');
fprintf(fp,'*           中子通量校验           *\n');
fprintf(fp,'************************************\n');
Phi_center = Phi(2); % 中心通量
fprintf(fp,'   校验点位置           模拟结果\n');
fprintf(fp,'  （R为球半径）     （中心通量归一）\n');
fprintf(fp,'------------------------------------\n');
fprintf(fp,'     0.25R             %.9g\n',Phi(uint32(0.25*R/(delta/2)))/Phi_center); % 对中心通量归一
fprintf(fp,'     0.50R             %.9g\n',Phi(uint32(0.50*R/(delta/2)))/Phi_center);
fprintf(fp,'     0.75R             %.9g\n',Phi(uint32(0.75*R/(delta/2)))/Phi_center);
fprintf(fp,'     1.00R             %.9g\n',Phi(uint32(1.00*R/(delta/2)))/Phi_center);
fclose(fp);

%% 绘制中子通量分布曲线
plot(r(2:2:2*K),Phi(2:2:2*K),'-r','linewidth',1.5);
xlabel('$r$','FontSize',16,'Interpreter','latex','FontWeight','bold');
ylabel('Fluence Rate $\phi$','FontSize',16,'Interpreter','latex','FontWeight','bold');
print('-depsc2','1D_sphere_neutron_flux.eps'); % 将中子通量分布曲线存入PDF文件保存