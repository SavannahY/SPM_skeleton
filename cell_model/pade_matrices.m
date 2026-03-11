function pade = pade_matrices(param,method)
% PADE_MATRICES 2nd or 3rd order Padé reduced diffusion model
%
% States per electrode:
%   PADE2 : 1 avg + 2 deviation
%   PADE3 : 1 avg + 3 deviation

switch upper(method)
    case 'PADE2'
        order = 2;
    case 'PADE3'
        order = 3;
    otherwise
        error('Unknown Padé method')
end

%% Diffusion time constants
tau_n = param.Rsn^2 / param.Dsn;
tau_p = param.Rsp^2 / param.Dsp;

%% Average concentration gains
k_avg_n =  1/(param.epsn * param.Acell * param.Ln * param.F);
k_avg_p = -1/(param.epsp * param.Acell * param.Lp * param.F);

%% Surface deviation gains
% correct physical scaling from diffusion boundary condition
k_dev_n =  tau_n/(3*param.epsn*param.Acell*param.Ln*param.F);
k_dev_p = -tau_p/(3*param.epsp*param.Acell*param.Lp*param.F);

%% Padé transfer functions (dimensionless form)

switch order

    case 2

        % Padé coefficients for surface deviation transfer function
        num_base = [-1/5  -3/455  -1/45045];
        den_base = [1     4/65    3/5005];

    case 3

        num_base = [-1/5  -1/119  -2/26775  -1/11486475];
        den_base = [1     6/85     2/1785    4/1044225];

end

%% Convert to dimensional transfer functions
% num_base / den_base are expressed in the dimensionless variable q = tau*s.
% Re-scale coefficients so tf2ss sees a transfer function in physical time s.
scale_n = 1 ./ (tau_n .^ (0:order));
scale_p = 1 ./ (tau_p .^ (0:order));

num_n = k_dev_n * (num_base .* scale_n);
num_p = k_dev_p * (num_base .* scale_p);

den_n = den_base .* scale_n;
den_p = den_base .* scale_p;

%% Convert to state-space

[Adev_n,Bdev_n,Cdev_n,Ddev_n] = tf2ss(num_n,den_n);
[Adev_p,Bdev_p,Cdev_p,Ddev_p] = tf2ss(num_p,den_p);

%% Add average concentration state

A_n = blkdiag(0,Adev_n);
A_p = blkdiag(0,Adev_p);

B_n = [k_avg_n ; Bdev_n];
B_p = [k_avg_p ; Bdev_p];

%% Combine electrodes

pade.A_mat = blkdiag(A_n,A_p);
pade.Bvector = [B_n ; B_p];

pade.Csurf_n = [1 Cdev_n];
pade.Dsurf_n = Ddev_n;

pade.Csurf_p = [1 Cdev_p];
pade.Dsurf_p = Ddev_p;

pade.order = order;
pade.nstates_electrode = order + 1;

end