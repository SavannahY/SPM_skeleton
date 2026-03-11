function [A, B] = fdm_matrices(param, nstates_electrode)
%FDM_MATRICES Build A and B matrices for solid diffusion (Alternative FDM)
%
% Target form used in your project:
%   xdot = A*x - B*Iapp
%
% A matches the Recitation Eq.(16)-(18) stencil (TA board).
% B is nonzero only at the surface node of each electrode.

N = nstates_electrode;
if N < 3
    error('nstates_electrode must be >= 2');
end

% Build dimensionless stencil M (same for both electrodes)
M = build_M_altFDM(N);

% Scale into each electrode block: (Ds/Rs^2)*(1/dr^2)*M
% An = (param.Dsn/(param.Rsn^2)) * (1/(param.dr_n^2)) * M;
% Ap = (param.Dsp/(param.Rsp^2)) * (1/(param.dr_p^2)) * M;


dr = 1/(N-1);

An = (param.Dsn/(param.Rsn^2)) * (1/(dr^2)) * M;
Ap = (param.Dsp/(param.Rsp^2)) * (1/(dr^2)) * M;



A = blkdiag(An, Ap);

% Build B (only surface node nonzero)
Bn = sparse(N,1);
Bp = sparse(N,1);

% Midterm sign convention for discharge I>0:
% g_n = +1 (anode), g_p = -1 (cathode)
g_n = +1;
g_p = -1;

% Bn(end) = (2/param.dr_n) * (1 + 1/(N-1)) * (g_n/(3*param.epsn*param.Acell*param.Ln*param.F));
% Bp(end) = (2/param.dr_p) * (1 + 1/(N-1)) * (g_p/(3*param.epsp*param.Acell*param.Lp*param.F));

Bn(end) = (2/dr) * (1 + 1/(N-1)) * (g_n/(3*param.epsn*param.Acell*param.Ln*param.F));
Bp(end) = (2/dr) * (1 + 1/(N-1)) * (g_p/(3*param.epsp*param.Acell*param.Lp*param.F));



B = [Bn; Bp];

end


function M = build_M_altFDM(N)
% Dimensionless stencil from Eq.(16)-(18)
% IMPORTANT: spdiags expects the +/-1 diagonal columns shifted/padded.

lower_des = zeros(N,1);     % desired A(i,i-1) coefficients stored at row i
main      = -2*ones(N,1);   % desired A(i,i)
upper_des = zeros(N,1);     % desired A(i,i+1) coefficients stored at row i

% Center node i=1: -6*c1 + 6*c2
main(1)      = -6;
upper_des(1) =  6;

% Interior i=2..N-1
if N > 2
    i = (2:N-1)';
    lower_des(i) = 1 - 1./(i-1);
    upper_des(i) = 1 + 1./(i-1);
end

% Surface i=N: 2*c_{N-1} - 2*c_N
lower_des(N) = 2;
main(N)      = -2;

% --- KEY FIX: align for spdiags ---
% In MATLAB spdiags construction, the column vectors are aligned with zeros
% at positions where the diagonal doesn't exist.
%
% For d = -1 (subdiagonal): element for (i,i-1) should be placed in row i-1 of the spdiags column.
% For d = +1 (superdiagonal): element for (i,i+1) should be placed in row i+1 of the spdiags column.
lower = [lower_des(2:N); 0];         % shift UP so row3->2 uses lower(2), etc.
upper = [0; upper_des(1:N-1)];       % shift DOWN so row1->2 uses upper(2), etc.

M = spdiags([lower, main, upper], [-1 0 1], N, N);

end
