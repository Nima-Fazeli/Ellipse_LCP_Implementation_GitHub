clc; clear; close all;

% Define end time (in seconds):
T_end = 6;

% Define simulation time step:
h = 1e-2;

% Define Coeff of friction
mu = 0.7;

% Define mass of block
m = 1; g = 9.81;

ellipse_res = 1e3; 
% Define Initial Conditions:
% The order of the vector is x,xdot,y,ydot,theta,thetadot
q_k = [  0;  5; 30*pi/180]; % configurations
v_k = [0.5;  0;         3]; % velocities

% q_k = [  0;  5; 60*pi/180];
% v_k = [3;  0;         1]; % Trace Plot

% Define Coeff of Rest:
e = 0.6;

% Define threshhold of COM velocity so that the cube stops bouncing
v_threshold = 0.1;

% Define size of square
a_e = 2;
b_e = 1;

% Collision Detection Resolution
epsilon = 1e-2;

% Define the time span for the ode45 function: (redundant)
tspan = [0 h];

% Start a while loop so that simulation runs till we hit the end time,
% initialize time here
t = 0;
x_vec = []; % gen coord
t_vec = []; % time
v_vec = []; % Vertex - not used here
f_vec = []; % Force Vec
l_vec = []; % lambda
b_vec = []; % Part of lcp
s_vec = []; % velocity
j_vec = []; % Jacobian
m_vec = []; % minimum detected point
i_vec = [];

E = [1;1];

% Pseudo Mass matrix:
M = [m,0,0;
    0,m,0;
    0,0,pi/4*m*(a_e*b_e^3+a_e^3*b_e)]; 

Minv  = inv(M); % inverse of mass matrix


while t < T_end
    
    v_kplus1 = v_k+Minv*(h*[0;-m*g;0]);
    q_kplus1 = q_k+h*v_kplus1;
    
    % Solve minimization problem to find closes point to the ground
    if q_kplus1(2)>max([a_e,b_e])*1.05
        min_vert = 1;
        beta_c = 0;
    else
        [min_vert,beta_c,X_min] = min_point_ellipse( q_kplus1(1:2), q_kplus1(3), [a_e,b_e], ellipse_res );
    end
    
    if (min_vert > 1e-4)
        % initialize the ICs for the next iteration
        q_k = q_kplus1;
        v_k = v_kplus1;
        t   = t+h;
        x_vec = [x_vec, q_k];
        s_vec = [s_vec, v_k];
        t_vec = [t_vec, t];
        v_vec = [v_vec, beta_c];
        f_vec = [f_vec, zeros(3,1)];
        b_vec = [b_vec, zeros(3,1)];
        l_vec = [l_vec, zeros(1,1)];
        j_vec = [j_vec, [zeros(3,1);zeros(3,1)]];
        m_vec = [m_vec, [q_kplus1(1:3); beta_c]];
        i_vec = [i_vec, min_vert];
    else
        
%         [min_vert,beta_c]=min_point_ellipse( q_kplus1(1:2), q_kplus1(3), [a_e,b_e], ellipse_res );
        
        % Calculate J_n and J_t: (contact jacobians)
        % J_n  - Normal
        n = calcN([X_min min_vert], q_kplus1); 
        
        % J_t  - Tangential
        d_col = calcD([X_min min_vert], q_kplus1); 
        
        D = [-d_col, d_col];
        
        % Define LCP:
        M_hat = [n'*Minv*n, n'*Minv*D,  0;
            D'*Minv*n, D'*Minv*D,  E;
            mu,       -E',  0];
        
        b  = h*Minv*[0;-m*g;0] + v_k;
        
        q_hat = [n'*b+n'*e*v_k; D'*b; 0];
        
        z = LCP(M_hat,q_hat);
        
        %         if abs(z(1))<0.05
        %             e = 0;
        % %             vlplus1 = zeros(3,1)+1e-3;
        %         end
        
        v_kplus1 = [Minv*n, Minv*D]*[z(1:3)]+b;
        q_k = q_k + v_kplus1*h;
        v_k = v_kplus1;
        
        t = t+h;
        v_vec = [v_vec, beta_c];
        x_vec = [x_vec, q_k];
        s_vec = [s_vec, v_k];
        t_vec = [t_vec, t];
        f_vec = [f_vec, [Minv*n, Minv*D]*z(1:3)];
        b_vec = [b_vec, z(1:3)];
        l_vec = [l_vec, z(4)];
        j_vec = [j_vec, [n;d_col]];
        m_vec = [m_vec, [q_kplus1(1:3); beta_c]];
        i_vec = [i_vec, min_vert];
    end
    
end

%%
figure(1)
clf
subplot(221)
plot(x_vec(1,:), x_vec(3,:))
xlabel 'x [m]'
ylabel 'y [m]'
grid on
subplot(222)
plot(t_vec, x_vec(1,:))
xlabel 'time [sec]'
ylabel 'x [m]'
grid on
subplot(223)
plot(t_vec, x_vec(2,:))
xlabel 'time [sec]'
ylabel 'y [m]'
grid on
subplot(224)
plot(t_vec, x_vec(3,:))
xlabel 'time [sec]'
ylabel '\theta [rad]'
grid on

%% Draw Ellipse

figure(2)
clf

% writerObj = VideoWriter('Ellipse_bounce_Full.avi');
% open(writerObj);

for i=1:10:size(x_vec,2)
    clf
    hold on
    
    ic=draw_ellipse_func(x_vec(1:2,i),x_vec(3,i), [a_e,b_e], ellipse_res);
    plot(ic(1,:),ic(2,:),'LineWidth',3);
    
    [x_l, y_l] = ellipse_min_line(x_vec(1:2,i), x_vec(3,i), m_vec(end,i), [a_e,b_e]);
    plot(x_l,y_l,'r','LineWidth',2);
    
    plot([-6,6],[0,0],'k')
    axis([-10,10,-1,8])
    axis equal
    grid on
    title(['Time: ',num2str(t_vec(i)),''])
    xlabel 'Horizontal displacement [m]'
    ylabel 'Vertical displacement [m]'
    pause(1e-2)
    
%         frame = getframe;
%         writeVideo(writerObj,frame);
end

% close(writerObj);

%%
figure(3)
subplot(221)
plot(t_vec,b_vec(1,:))
xlabel 'time [sec]'
ylabel 'c_n [Normal Impact Impulse]'
grid on
subplot(222)
plot(t_vec,b_vec(2,:),'--')
xlabel 'time [sec]'
ylabel '\beta_1'
grid on
subplot(223)
plot(t_vec,b_vec(3,:),'--')
xlabel 'time [sec]'
ylabel '\beta_2'
grid on
subplot(224)
plot(t_vec,l_vec(1,:))
xlabel 'time [sec]'
ylabel 'Lambda [Measure of Relative veclocity]'
grid on

figure(4)
clf
subplot(311)
plot(t_vec,f_vec(1,:),'--')
xlabel 'time [sec]'
ylabel 'Tangential Force'
grid on
subplot(312)
plot(t_vec,f_vec(2,:),'--')
xlabel 'time [sec]'
ylabel 'Normal Force'
grid on
subplot(313)
plot(t_vec,f_vec(3,:),'--')
xlabel 'time [sec]'
ylabel 'Torque'
grid on

%% Trace plot:

figure_ellipse = figure(10);
clf
% x_v = [v_vec(1,:)',v_vec(3,:)'];
% y_v = [v_vec(2,:)',v_vec(4,:)'];


for i=1:100:size(x_vec,2)
    hold on
    ic=draw_ellipse_func(x_vec(1:2,i),x_vec(3,i), [a_e,b_e], ellipse_res);
    plot(ic(1,:),ic(2,:),'--k','LineWidth',1);
end

    ic=draw_ellipse_func(x_vec(1:2,i),x_vec(3,i), [a_e,b_e], ellipse_res);
    plot(ic(1,:),ic(2,:),'k','LineWidth',2);

plot([-6,15],[0,0],'k')
axis([-3,9,-1,8])
% axis equal
grid on

xlabel 'Horizontal displacement [m]'
ylabel 'Vertical displacement [m]'
title 'Trace of Trajectory of the 2-D Ellipse'

%% Saving the data from the runs for sysid
% save('ellipse_full','t_vec','x_vec','a_e','b_e','ellipse_res')
