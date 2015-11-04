%% Fixed-Point Q21.10 Linear filter Kalman (for navigation) 
% Example:
% x - current position in current sample time;
% vx - rate of current sample time;
% dt - sample time;
% R - value noise of measurement coordinate X;
% Rvx - value noise of measurement rate Vx;
%%
% All input data represent coordinate and rate at a fixed point Q21.10
%%
function [y] = KalmanC(x, vx, Rx, Rvx, dt)
    %#codegen
    Q = 1024;			% measurement noise
    persistent x_est x_est_vx p_est p_est_vx p_est_diag_0_1 k11 k12 k21 k22 x_tek vx_tek %p_est_diag_1_0            % Initial state conditions
    if isempty(x_est)
        p_est = int32(0);             % predicted error to X
        p_est_vx = int32(0);          % predicted error to Vx
        p_est_diag_0_1 = int32(0);    % diagonal element of a covariation matrix P_EST on the previous step
        k11 = int32(0);
        k12 = int32(0);
        k21 = int32(0);
        k22 = int32(0);
        x_est = x;
        x_est_vx = vx;
    end
    % time update - prediction
    % Prediction of a condition X and Vx
    x_prd = x_est + (dt * x_est_vx) / 1024;
    x_prd_vx = x_est_vx;

    % Calculation diagonal value of a covariation matrix P_EST on the previous step
    p_prd_dVx = p_est_diag_0_1 + (dt * p_est_vx) / 1024;

    % Prediction of an error X and Vx
    p_prd = p_est + (((dt * dt) / 1024) * p_est_vx) / 1024 + Q;
    p_prd_vx = p_est_vx + Q;

    % measurement update - correction
    % Addition of an error of measurement Rx to the predicted error (specially "noisering" on X and Vx).
    S = p_prd + Rx;							% Noiser of the predicted error on X error Rx
    Svx = p_prd_vx + Rvx;					% Noiser of the predicted error on Vx error Rvx

    % Reception of value of the predicted error on X and Vx without the added noise
    B = p_prd;								% Value of the predicted error on X without noise
    Bvx = p_prd_vx;							% Value of the predicted error on Vx without noise
	
        % Matrix Filter of Kalman (not transponation).
    k11 = (-(p_prd_dVx * p_prd_dVx - Svx * B) * 16) / ((S * Svx - p_prd_dVx * p_prd_dVx) / 64);		% KalmanGainX
    k12 = (-(p_prd_dVx * Bvx - Svx * p_prd_dVx) * 16) /((S * Svx - p_prd_dVx * p_prd_dVx) / 64);	% dt*Vx

    k21 = ((S * p_prd_dVx - p_prd_dVx * B) * 16) / ((S * Svx - p_prd_dVx * p_prd_dVx) / 64);		% xVx
    k22 = ((S * Bvx - p_prd_dVx * p_prd_dVx) * 16) / ((S * Svx - p_prd_dVx * p_prd_dVx) / 64);		% KalmanGainVx

    % Definition predicted condition on X and Vx. It is made on the basis of disclosing of matrixes (a linear kind).
    x_est = x_prd + (k11 * (x - x_prd) + k21 * (vx - x_prd_vx)) / 1024;
    x_est_vx = x_prd_vx + (k12 * (x - x_prd) + k22 * (vx - x_prd_vx)) / 1024;
    
    % Definition predicted errors on X and Vx. It is made on the basis of disclosing of matrixes (a linear kind).
    p_est = p_prd - (k11 * p_prd + k21 * p_prd_dVx) / 1024;
    p_est_vx = p_prd_vx - (k12 * p_prd_dVx + k22 * p_prd_vx) / 1024;
    
    % Calculation of a collateral diagonal of a covariation matrix condition errors (value of a collateral diagonal). They leave identical, but whether always?
    p_est_diag_0_1 = p_prd_dVx - (k12 * p_prd + k22 * p_prd_dVx) / 1024;							% p_est_diag_0_1 (0,1)
    %p_est_diag_1_0 = p_prd_dVx - (k11 * p_prd_dVx + k21 * p_prd_vx);						% p_est_diag_1_0 (1,0)
    
	%Filtered data
    y = x_est;
end
