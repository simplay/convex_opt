% Convex Optimization - Project 2
% MICHAEL SINGLE
% 08-917-445
function [demosaicedImg] = demosaicing_michael_single(mosaiced, Omega, lambda, iterations)
% DEMOSAIC demosaices a given mosaiced image by solving a convex opt.
%          problem relying on primal dual solver 
%          and a fixed numer of iteration.
%
% @author: MICHAEL SINGLE
%          08-917-445
%
% @param mosaiced M x N x 3 Color Image encoded as a mosaiced image
%        according to bayer filter theory.
% @param Omega bayer filter tensor M x N x 3
%        Omega(:,:,1) is RED color mask
%        Omega(:,:,2) is GREEN color mask
%        Omega(:,:,3) is BLUE color mask
% @param lambda [Float] regularization factor to control exactness of result and the smoothness 
%        when solving the demosaicing convex minimization problem.
% @param iterations [Integer] number of iterations
% @return demosaicedImg M x N x 3 Color Image.

    % parameters 
    tau = 0.001; 
    sigma = 0.001;
    theta = 0.5;

    % initial guesses 
    Rx_n = mosaiced(:,:,1); 
    Gx_n = mosaiced(:,:,2);
    Bx_n = mosaiced(:,:,3);

    Ry_n = ones([size(Rx_n), 2]);
    Gy_n = ones([size(Gx_n), 2]);
    By_n = ones([size(Bx_n), 2]);
    
    Rx_tilde_n = Rx_n;
    Gx_tilde_n = Gx_n;
    Bx_tilde_n = Bx_n;
    
    h = waitbar(0,'Progress Gradient Descend');
    for i = 1:iterations
       waitbar(i/iterations) 
       [Ry_n, Rx_n, Rx_tilde_n] = primal_dual_solver(Ry_n, Rx_n, Rx_tilde_n, mosaiced(:,:,1), Omega(:,:,1), lambda, tau, sigma, theta);
       [Gy_n, Gx_n, Gx_tilde_n] = primal_dual_solver(Gy_n, Gx_n, Gx_tilde_n, mosaiced(:,:,2), Omega(:,:,2), lambda, tau, sigma, theta);
       [By_n, Bx_n, Bx_tilde_n] = primal_dual_solver(By_n, Bx_n, Bx_tilde_n, mosaiced(:,:,3), Omega(:,:,3), lambda, tau, sigma, theta);
    end
    close(h); 
    demosaicedImg = mat2Img(uR,uG,uB);
end

function [y_n_p_1, x_n_p_1, x_tilde_n_p_1] = primal_dual_solver(y_n, x_n, x_tilde_n, g, Omega, lambda, tau, sigma, theta)
% @param y_n current partial y iterative solution of primal dual formulation.
% @param x_n current partial x iterative solution of primal dual formulation.
% @param x_tilde_n current iterative solution of primal dual formulation.
% @param g initial mosaiced image of color channel corresponding to given u.
% @param Omega bayer filter mask M x N for a certain color channel
% @param lambda [Float] regularization factor to control exactness of result and the smoothness 
%        when solving the demosaicing convex minimization problem.
% @param tau learning rate for gradient descend.
% @param eps small perturbation to avoid division by zero issues.
% @return next iteration of color channel u.
	y_n_p_1 = y_n_plus_1_for(y_n, x_tilde_n, sigma);
    x_n_p_1 = x_n_plus_1_for(x_n, y_n_p_1, lambda, tau, Omega, g);
    x_tilde_n_p_1 = x_tilde_n_plus_1_for(x_n_p_1, x_n, theta);
end

function y_n_p_1 = y_n_plus_1_for(y_n, x_tilde, sigma)
    grad_tilde_x = grad_of(x_tilde);
    y_nominator = y_n + sigma*grad_tilde_x;
    norm_y_nominator = sqrt(y_nominator(:,:,1).^2 + y_nominator(:,:,2).^2);
    y_denominator = max(1, norm_y_nominator);

    y_n_p_1 = zeros(size(y_n));
    y_n_p_1(:,:,1) = y_nominator(:,:,1)./y_denominator;
    y_n_p_1(:,:,2) = y_nominator(:,:,2)./y_denominator;
end

function x_n_p_1 = x_n_plus_1_for(x_n, y_n_p_1, lambda, tau, Omega, g)
    div_y_n_p_1 = div_of(y_n_p_1);
    x_n_p_1_nominator = x_n + tau*div_y_n_p_1 + tau*lambda*(Omega.*g);
    x_n_p_1_denominator = 1 + tau*lambda*Omega;
    x_n_p_1 = x_n_p_1_nominator./x_n_p_1_denominator;
end

function x_tilde_n_p_1 = x_tilde_n_plus_1_for(x_n_p_1, x_n, theta)
    x_tilde_n_p_1 = x_n_p_1 + theta*(x_n_p_1 - x_n);
end

function grad_f = grad_of(f)
% foreward difference scheme

    grad_f = zeros([size(f), 2]);
    df_dx = f([2:end, end],:)-f(:,:);
    df_dy = f(:,[2:end, end])-f(:,:);
    grad_f(:,:,1) = df_dx;
    grad_f(:,:,2) = df_dy;
end

function div_v = div_of(v)
   grad_v_x = grad_of(v(:,:,1));
   grad_v_y = grad_of(v(:,:,2));
   div_v = grad_v_x(:,:,1) + grad_v_y(:,:,2);
end
