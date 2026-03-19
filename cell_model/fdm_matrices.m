function [A, B] = fdm_matrices(param, nstates_electrode)
% FDM_MATRICES - Construct Finite Difference Method matrices for SPM
%
% DESCRIPTION:
%   Constructs the A and B matrices for the ODE system: dc/dt = A*c + B*I
%   for spherical diffusion using Finite Difference Method.
%   State variable is concentration c_s [mol/m³], not stoichiometry theta.
%
% PHYSICS:
%   PDE: ∂c_{s,j}/∂t = D_{s,j}/r² ∂/∂r(r² ∂c_{s,j}/∂r), j = n, p
%   BC at r=0: ∂c_{s,j}/∂r = 0 (symmetry)
%   BC at r=R: D_{s,j} ∂c_{s,j}/∂r = -I_app * g(I_app) / (a_{s,j} * A_cell * L_j * F)
%   where a_{s,j} = 3*eps_j/R_{s,j} and g(I_app) = -1 for cathode, +1 for anode
%
% INPUTS:
%   param             - Parameter structure with diffusivities and geometry
%   nstates_electrode - Number of spatial nodes per electrode [-]
%
% OUTPUTS:
%   A - State matrix [1/s] for dc/dt = A*c + B*I
%   B - Input matrix [mol/(m³·A·s)] for dc/dt = A*c + B*I

    % Construct diffusion matrix for anode
    j_values = (1:nstates_electrode)';

    above_diagonals = (1.0 + 0.5./(j_values-1)).^2;
    on_diagonals = (2.0 + 0.5./(j_values-1).^2);
    below_diagonals = (1.0 - 0.5./(j_values-1)).^2;

    D_n = spdiags([circshift(below_diagonals, -1), -on_diagonals, circshift(above_diagonals, 1)], [-1, 0, 1], nstates_electrode, nstates_electrode);
    D_n(end, end-1) = -D_n(end,end); D_n(1, 1) = -6.0; D_n(1, 2) = 6.0;
    D_n = (param.Dsn./(param.Rsn.^2)./(param.dr_n.^2)).*D_n;

    % Construct diffusion matrix for cathode
    D_p = spdiags([circshift(below_diagonals, -1), -on_diagonals, circshift(above_diagonals, 1)], [-1, 0, 1], nstates_electrode, nstates_electrode);
    D_p(end, end-1) = -D_p(end,end); D_p(1, 1) = -6.0; D_p(1, 2) = 6.0;
    D_p = (param.Dsp./(param.Rsp.^2)./(param.dr_p.^2)).*D_p;
    
    % Complete A matrix (state matrix)
    A = blkdiag(D_n, D_p); 

    % Construct B vector (input matrix)
    % BC: D_s * dc/dr|_{r=R} = -I_app * g(I_app) / (a_s * A_cell * L * F)
    % where a_s = 3*eps/R_s and g = +1 for anode (n), -1 for cathode (p)
    % Discretized: incorporates the flux BC into the last node equation
    Nr = nstates_electrode;
    B = zeros(2*nstates_electrode, 1);
    % Anode: g(I_app) = +1, flux enters particle during discharge
    B(nstates_electrode) = 2.0.*((1.0 + 0.5./(Nr-1)).^2).*1.0./(param.epsn.*param.Ln.*param.Acell.*param.F)./(3.0.*param.dr_n);
    % Cathode: g(I_app) = -1, flux leaves particle during discharge
    B(end) = -2.0.*((1.0 + 0.5./(Nr-1)).^2).*1.0./(param.epsp.*param.Lp.*param.Acell.*param.F)./(3.0.*param.dr_p);

end
