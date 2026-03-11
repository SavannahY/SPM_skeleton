function [sol, t_eval, coeff] = run_first_principles_ecm_casadi(param, t_data, I_data, T_amb, solver_opts)
% RUN_FIRST_PRINCIPLES_ECM_CASADI - 终极核动力版 (God Mode)
% 采用 CasADi 的 CVODES 求解器 + mapaccum 技术，让 Padé 模型实现毫秒级求解！
    
    import casadi.*

    if nargin < 5
        solver_opts = struct();
    end

    param.T_amb = 273.15 + T_amb; 

    if ~isfield(solver_opts, 'linearization_soc') || isempty(solver_opts.linearization_soc)
        solver_opts.linearization_soc = param.SOC_IC;
    end

    % 1. 计算 Padé 系数矩阵
    coeff = pade_spm_coeffs(param, solver_opts.linearization_soc, param.T_amb);

    % 2. 初始状态
    [~, csn0, csp0] = conc_initial_sd(param, struct('nstates_electrode',1));
    x0 = zeros(8,1);
    x0(7) = csn0; 
    x0(8) = csp0;

    %% ========================================================
    %  CASADI 核心组装 (建造高铁) 
    %  ========================================================
    % 1. 定义符号变量 (Symbolic Variables)
    x = MX.sym('x', 8);
    I_sym = MX.sym('I');
    
    % 2. 定义导数 (RHS)
    dx = MX.zeros(8,1);
    dx(1:3) = coeff.n.A * x(1:3) + coeff.n.B * I_sym;
    dx(4:6) = coeff.p.A * x(4:6) + coeff.p.B * I_sym;
    dx(7) = -I_sym / (param.epsn * param.Acell * param.Ln * param.F);
    dx(8) =  I_sym / (param.epsp * param.Acell * param.Lp * param.F);
    
    % 3. 构建 DAE 系统
    ode = struct('x', x, 'p', I_sym, 'ode', dx);
    
% 4. 创建单步积分器 (适配最新版 CasADi 语法，消除 Warning)
    dt = t_data(2) - t_data(1);
    opts = struct(); % 清空被弃用的 'tf' 选项
    % 新语法：直接把时间起点 0 和终点 dt 作为参数传进函数里
    F = integrator('F', 'cvodes', ode, 0, dt, opts);
    
    % 5. 【核武器】使用 mapaccum 进行全序列映射
    % 它会直接在 C++ 层面完成整个 UDDS 的时间循环！
    N = length(t_data);
    F_mapaccum = F.mapaccum('F_map', N-1);
    
    % 6. 一键点火执行
    res = F_mapaccum('x0', x0, 'p', I_data(1:end-1)');
    
    % 7. 提取结果并把初始状态拼回去
    x_res = [x0, full(res.xf)];
    
    %% ========================================================
    %  极速后处理 (纯向量化计算电压)
    %  ========================================================
    x_res = x_res.'; % 转置为 [N x 8]
    t_eval = t_data(:);
    I_eval = I_data(:);

    delta_csn = (coeff.n.C * x_res(:, 1:3).').'; 
    delta_csp = (coeff.p.C * x_res(:, 4:6).').';
    
    csn_surf = csn0 + delta_csn;
    csp_surf = csp0 + delta_csp;

    theta_surf_n = min(max(csn_surf ./ param.csn_max, 1e-8), 1 - 1e-8);
    theta_surf_p = min(max(csp_surf ./ param.csp_max, 1e-8), 1 - 1e-8);

    sol.ocp_n = U_n(theta_surf_n, param.T_amb, param);
    sol.ocp_p = U_p(theta_surf_p, param.T_amb, param);
    sol.eta_n = eta_anode(theta_surf_n, param.T_amb, I_eval, param);
    sol.eta_p = eta_cathode(theta_surf_p, param.T_amb, I_eval, param);

    sol.V_oc = sol.ocp_p - sol.ocp_n;
    sol.V_cell = sol.V_oc + sol.eta_p - sol.eta_n - param.R0 .* I_eval;

    sol.theta_s_n_surf = theta_surf_n;
    sol.theta_s_p_surf = theta_surf_p;
    sol.soc_bulk_n = (x_res(:,7) - param.theta_n_0) ./ (param.theta_n_100 - param.theta_n_0);
    sol.soc_bulk_p = (x_res(:,8) - param.theta_p_0) ./ (param.theta_p_100 - param.theta_p_0);
end