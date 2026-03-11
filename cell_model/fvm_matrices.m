function [A, B] = fvm_matrices(param, nstates_electrode)
%FVM_MATRICES Finite-volume matrices for SPM solid diffusion.
%
%   The FVM states are shell-averaged concentrations over Nr-1 control
%   volumes bounded by the radial interfaces stored in param.r_n / param.r_p.
%   The assembled model has the form
%
%       dx/dt = A*x - B*I
%
%   where I > 0 corresponds to discharge.
%
%   The discretization enforces lithium conservation by construction because
%   it is obtained from a shell-wise flux balance.

    arguments
        param struct
        nstates_electrode (1,1) {mustBeInteger,mustBePositive}
    end

    Ncv_expected = numel(param.r_n) - 1;
    if nstates_electrode ~= Ncv_expected
        error(['For FVM, nstates_electrode must equal Nr-1. ', ...
               'Expected %d, got %d.'], Ncv_expected, nstates_electrode);
    end

    [A_n, B_n] = build_fvm_electrode( ...
        param.Dsn, param.Rsn, param.epsn, param.Acell, param.Ln, ...
        param.F, param.r_n(:), param.dr_n, +1.0);

    [A_p, B_p] = build_fvm_electrode( ...
        param.Dsp, param.Rsp, param.epsp, param.Acell, param.Lp, ...
        param.F, param.r_p(:), param.dr_p, -1.0);

    A = blkdiag(A_n, A_p);
    B = [B_n; B_p];
end

function [A_j, B_j] = build_fvm_electrode(Ds, Rs, eps_s, Acell, Lj, F, r_edges, dr, g_sign)
%BUILD_FVM_ELECTRODE Build one-electrode FVM diffusion operator.
%
% r_edges are the normalized shell interfaces in [0, 1].

    Ncv = numel(r_edges) - 1;
    dV = r_edges(2:end).^3 - r_edges(1:end-1).^3;   % 3x normalized shell volumes
    alpha = Ds / Rs^2;

    A_j = spalloc(Ncv, Ncv, max(1, 3*Ncv));

    if Ncv == 1
        % Single shell: no internal diffusion term, only boundary flux input.
        A_j(1,1) = 0.0;
    else
        % First shell: left interface contribution is zero because r = 0.
        coeff_r = 3.0 * alpha * r_edges(2)^2 / (dV(1) * dr);
        A_j(1,1) = -coeff_r;
        A_j(1,2) =  coeff_r;

        % Interior shells.
        for i = 2:(Ncv-1)
            coeff_l = 3.0 * alpha * r_edges(i)^2   / (dV(i) * dr);
            coeff_r = 3.0 * alpha * r_edges(i+1)^2 / (dV(i) * dr);

            A_j(i,i-1) =  coeff_l;
            A_j(i,i)   = -(coeff_l + coeff_r);
            A_j(i,i+1) =  coeff_r;
        end

        % Last shell: right interface is handled through the boundary-flux input.
        coeff_l = 3.0 * alpha * r_edges(end-1)^2 / (dV(end) * dr);
        A_j(end,end-1) =  coeff_l;
        A_j(end,end)   = -coeff_l;
    end

    % Input coefficient from the surface flux balance.
    B_j = sparse(Ncv, 1);
    B_j(end) = g_sign / (eps_s * Acell * Lj * F * dV(end));
end
