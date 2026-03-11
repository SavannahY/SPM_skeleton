function dUdtheta = dUdtheta_fd(ocp_fun, theta0, T_core, param)
%DUDTHETA_FD Central-difference derivative of an OCP curve with respect to stoichiometry.
    h = 1e-6;
    theta_m = max(1e-6, min(1 - 1e-6, theta0 - h));
    theta_p = max(1e-6, min(1 - 1e-6, theta0 + h));
    dUdtheta = (ocp_fun(theta_p, T_core, param) - ocp_fun(theta_m, T_core, param)) / (theta_p - theta_m);
end
