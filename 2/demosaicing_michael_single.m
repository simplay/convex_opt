% Convex Optimization - Project 2
% MICHAEL SINGLE
% 08-917-445
function [demosaicedImg] = demosaicing_michael_single(mosaiced, Omega, lambda, iterations, input)
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

    % parameter
    K_a = sqrt(4);
    tau = 10*0.2*1e-3; 
    sigma = 1/(tau*K_a);
    disp(['tau*sigma=',num2str(tau*sigma)]);
    theta = 0.0;

    % initial guesses 
    Rx_n = mosaiced(:,:,1); 
    Gx_n = mosaiced(:,:,2);
    Bx_n = mosaiced(:,:,3);

    Ry_n = zeros([size(Rx_n), 2]);
    Gy_n = zeros([size(Gx_n), 2]);
    By_n = zeros([size(Bx_n), 2]);
    
    Rx_tilde_n = Rx_n;
    Gx_tilde_n = Gx_n;
    Bx_tilde_n = Bx_n;
    %figure
    
    %h = waitbar(0,'Progress Gradient Descend');
    for i = 1:iterations
       %waitbar(i/iterations) 
       before_r_x = Rx_n;
       before_g_x = Gx_n;
       before_b_x = Bx_n;
       
       before_r_y = Ry_n;
       before_g_y = Gy_n;
       before_b_y = By_n;
       
       before_r_x_ = Rx_tilde_n;
       before_g_x_ = Gx_tilde_n;
       before_b_x_ = Bx_tilde_n;
       
       [Ry_n, Rx_n, Rx_tilde_n] = primal_dual_solver(before_r_y, before_r_x, before_r_x_, mosaiced(:,:,1), Omega(:,:,1), lambda, tau, sigma, theta);
       [Gy_n, Gx_n, Gx_tilde_n] = primal_dual_solver(before_g_y, before_g_x, before_g_x_, mosaiced(:,:,2), Omega(:,:,2), lambda, tau, sigma, theta);
       [By_n, Bx_n, Bx_tilde_n] = primal_dual_solver(before_b_y, before_b_x, before_b_x_, mosaiced(:,:,3), Omega(:,:,3), lambda, tau, sigma, theta);
       
       a = input(:,:,1);
       b = Rx_tilde_n;
       
      
       
       
       %plot(i,norm(a(:)-b(:)),'.');
       %hold on;
       %imagesc(Rx_n);
       %drawnow
    end
    %close(h); 
    demosaicedImg = mat2Img(Rx_tilde_n, Gx_tilde_n, Bx_tilde_n);
    %demosaicedImg = mat2Img(Rx_tilde_n, Rx_tilde_n, Rx_tilde_n);
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
    y_nominator = y_n - sigma*grad_tilde_x;
    norm_y_nominator = sqrt(y_nominator(:,:,1).^2 + y_nominator(:,:,2).^2);
    y_denominator = max(1, norm_y_nominator);

    y_n_p_1 = zeros(size(y_n));
    y_n_p_1(:,:,1) = y_nominator(:,:,1)./y_denominator;
    y_n_p_1(:,:,2) = y_nominator(:,:,2)./y_denominator;
end

function x_n_p_1 = x_n_plus_1_for(x_n, y_n_p_1, lambda, tau, Omega, g)
    div_y_n_p_1 = div_of(y_n_p_1);
    x_n_p_1_nominator = x_n - tau*div_y_n_p_1 + tau*lambda*(Omega.*g);
    x_n_p_1_denominator = 1 + tau*lambda*Omega;
    x_n_p_1 = x_n_p_1_nominator./x_n_p_1_denominator;
end

function x_tilde_n_p_1 = x_tilde_n_plus_1_for(x_n_p_1, x_n, theta)
    x_tilde_n_p_1 = x_n_p_1 + theta*(x_n_p_1 - x_n);
end

function grad_f = grad_of(f)

    grad_f = zeros([size(f), 2]);
    
    % foreward differences   
    df_dx = f([2:end, end],:)-f;
    df_dy = f(:,[2:end, end])-f;
    
    % backward differences
    %df_dx = f-f([1,1:end-1],:);
    %df_dy = f-f(:,[1,1:end-1]);
    
    % central differences
    %df_dx = ( f([2:end,end],:) - f([1,1:end-1],:) )/2;
    %df_dy = ( f(:,[2:end,end]) - f(:,[1,1:end-1]) )/2;
    
    grad_f(:,:,1) = df_dx;
    grad_f(:,:,2) = df_dy;

end

function grad_f = grad_of2(f)

    grad_f = zeros([size(f), 2]);
    
    % foreward differences   
    %df_dx = f([2:end, end],:)-f;
    %df_dy = f(:,[2:end, end])-f;
    
    % backward differences
    df_dx = f-f([1,1:end-1],:);
    df_dy = f-f(:,[1,1:end-1]);
    
    % central differences
    %df_dx = ( f([2:end,end],:) - f([1,1:end-1],:) )/2;
    %df_dy = ( f(:,[2:end,end]) - f(:,[1,1:end-1]) )/2;
    
    grad_f(:,:,1) = df_dx;
    grad_f(:,:,2) = df_dy;

end

function div_v = div_of(v)
   grad_v_x = grad_of2(v(:,:,1));
   grad_v_y = grad_of2(v(:,:,2));
   div_v = grad_v_x(:,:,1) + grad_v_y(:,:,2);
end
