function coeff = pade_spm_coeffs(param, linearization_soc, T_core)
%PADE_SPM_COEFFS Build the Padé-reduced single-particle model coefficients.
%
% This implementation follows the 3rd-order Padé reduction of the spherical
% diffusion transfer function used in
%   Prasad & Rahn, DSCC2012-MOVIC2012-8607.
%
% Practical note:
%   The paper is a local linear model about a chosen SOC. Here we reuse the
%   midterm skeleton's nonlinear OCP and Butler-Volmer output equations, but
%   keep the Padé diffusion reduction for a compact, fast model that works on
%   the supplied HPPC/UDDS/C/20 data.

    if nargin < 2 || isempty(linearization_soc)
        linearization_soc = param.SOC_IC;
    end
    if nargin < 3 || isempty(T_core)
        T_core = param.T_amb;
    end

    % Linearization stoichiometries from the same SOC mapping used by the
    % midterm skeleton.
    theta_n0 = linearization_soc * (param.theta_n_100 - param.theta_n_0) + param.theta_n_0;
    theta_p0 = param.theta_p_0 - linearization_soc * (param.theta_p_0 - param.theta_p_100);

    csn0 = theta_n0 * param.csn_max;
    csp0 = theta_p0 * param.csp_max;

    % Specific interfacial surface areas.
    asn = 3.0 * param.epsn / param.Rsn;
    asp = 3.0 * param.epsp / param.Rsp;

    % Third-order Padé coefficients for c_surf / I.
    coeff.n = local_particle_coeffs(param.Dsn, param.Rsn, asn, param.Acell, param.Ln, param.F, -1);
    coeff.p = local_particle_coeffs(param.Dsp, param.Rsp, asp, param.Acell, param.Lp, param.F, +1);

    coeff.csn0 = csn0;
    coeff.csp0 = csp0;
    coeff.theta_n0 = theta_n0;
    coeff.theta_p0 = theta_p0;
    coeff.asn = asn;
    coeff.asp = asp;

    % Optional paper-style linearized ECM quantities.
    dUp_dtheta = dUdtheta_fd(@U_p, theta_p0, T_core, param);
    dUn_dtheta = dUdtheta_fd(@U_n, theta_n0, T_core, param);
    dUp_dc = dUp_dtheta / param.csp_max;
    dUn_dc = dUn_dtheta / param.csn_max;

    % Exchange current densities at the linearization point.
    i0p = param.k0p * param.F * sqrt(param.ce * csp0 * (param.csp_max - csp0));
    i0n = param.k0n * param.F * sqrt(param.ce * csn0 * (param.csn_max - csn0));

    % With alpha_a + alpha_c = 1 in the provided parameter set.
    Rctp = param.Rg * T_core / (param.F * i0p);
    Rctn = param.Rg * T_core / (param.F * i0n);

    Cplus = 21 * dUp_dc;
    Cminus = 21 * dUn_dc;

    alpha1 = Cplus  / (param.Acell * param.F * asp * param.Lp * param.Rsp);
    alpha2 = param.Dsp / (param.Rsp^2);
    beta1  = Cminus / (param.Acell * param.F * asn * param.Ln * param.Rsn);
    beta2  = param.Dsn / (param.Rsn^2);

    K = Rctp / (asp * param.Acell * param.Lp) ...
      - Rctn / (asn * param.Acell * param.Ln) ...
      - param.R0;

    coeff.linearized.K = K;
    coeff.linearized.alpha1 = alpha1;
    coeff.linearized.alpha2 = alpha2;
    coeff.linearized.beta1 = beta1;
    coeff.linearized.beta2 = beta2;
    coeff.linearized.Cplus = Cplus;
    coeff.linearized.Cminus = Cminus;
    coeff.linearized.Rctp = Rctp;
    coeff.linearized.Rctn = Rctn;
    coeff.linearized.rc = local_rc_table(K, alpha1, alpha2, beta1, beta2);
end

function part = local_particle_coeffs(Ds, Rs, asj, Acell, Lj, F, signSurf)
    omega = Ds / (Rs^2);
    scale = signSurf * 21.0 / (asj * F * Acell * Lj * Rs);

    part.a1 = 3465.0 * omega^2;
    part.a2 = 189.0  * omega;
    part.b0 = scale * 495.0 * omega^2;
    part.b1 = scale * 60.0  * omega;
    part.b2 = scale;

    % Companion realization for (b2 s^2 + b1 s + b0)/(s^3 + a2 s^2 + a1 s)
    part.A = [0 1 0; 0 0 1; 0 -part.a1 -part.a2];
    part.B = [0; 0; 1];
    part.C = [part.b0 part.b1 part.b2];
    part.D = 0;
end

function rc = local_rc_table(K, alpha1, alpha2, beta1, beta2)
% Paper Table 2. These values can become non-passive depending on OCP slope
% sign and the chosen linearization point.
    rc.R1 = K;
    rc.C1 = 7.0 / (alpha1 + beta1);

    rc.C2 = 9.6246 / alpha1;
    rc.R2 = 0.0051 * alpha1 / alpha2;
    rc.C3 = 1.3277 / alpha1;
    rc.R3 = 0.0045 * alpha1 / alpha2;

    rc.C4 = 9.6246 / beta1;
    rc.R4 = 0.0051 * beta1 / beta2;
    rc.C5 = 1.3277 / beta1;
    rc.R5 = 0.0045 * beta1 / beta2;
end
