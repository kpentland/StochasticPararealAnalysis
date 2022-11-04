function [t,u] = RK(t,u0,f,method)
%This function is a Runge-Kutta ODE solver that uses a specified RK method
% depending on user preference.

%Inputs:
% t     : Time interval (i.e. t = (0:0.1:10))
% y0    : Initial condition (row vector)
% f     : Function handle of ODEs to be solved
% method: Runge-Kutta method to be used

%Outputs:
% t    : Time interval
% y    : Solution in matrix form (time steps by dimension of problem)


%Select the method to be used.
switch method
    case 'RK1' %Euler's method
        a = 0;
        b = 1; 
        c = 0;
    case 'RK2' %midpoint method
        a = [0,0;0.5,0];
        b = [0,1];
        c = [0,0.5];
    case 'RK3' %Kutta's third-order method
        a = [0,0,0;0.5,0,0;-1,2,0];
        b = [1/6,2/3,1/6];
        c = [0,0.5,1];
    case 'RK4' %classic fourth-order method
        a = [0,0,0,0;0.5,0,0,0;0,0.5,0,0;0,0,1,0];
        b = [1/6,1/3,1/3,1/6];
        c = [0,0.5,0.5,1];
    case 'RK5' %Butcher's fifth-order method (there are many - see his 2008 book)
        a = [0,0,0,0,0,0;0.25,0,0,0,0,0;0.125,0.125,0,0,0,0;0,0,0.5,0,0,0;3/16,-3/8,3/8,9/16,0,0;-3/7,8/7,6/7,-12/7,8/7,0];
        b = [7,0,32,12,32,7]./90;
        c = [0,0.25,0.25,0.5,0.75,1];
    case 'RK6' %Butcher's sixth-order method (again there are many)
        a = [0,0,0,0,0,0,0;1/3,0,0,0,0,0,0;0,2/3,0,0,0,0,0;1/12,1/3,-1/12,0,0,0,0;25/48,-55/24,35/48,15/8,0,0,0;3/20,-11/24,-1/8,1/2,1/10,0,0;-261/260,33/13,43/156,-118/39,32/195,80/39,0];
        b = [13/200,0,11/40,11/40,4/25,4/25,13/200];
        c = [0,1/3,2/3,1/3,5/6,1/6,1];
    case 'RK7' %Butcher's seventh-order method (again there are many)
        a = [0,0,0,0,0,0,0,0,0;1/6,0,0,0,0,0,0,0,0;0,1/3,0,0,0,0,0,0,0;1/8,0,3/8,0,0,0,0,0,0;148/1331,0,150/1331,-56/1331,0,0,0,0,0;-404/243,0,-170/27,4024/1701,10648/1701,0,0,0,0;2466/2401,0,1242/343,-19176/16807,-51909/16807,1053/2401,0,0,0;5/154,0,0,96/539,-1815/20384,-405/2464,49/1144,0,0;-113/32,0,-195/22,32/7,29403/3584,-729/512,1029/1408,21/16,0];
        b = [0,0,0,32/105,1771561/6289920,243/2560,16807/74880,77/1440,11/270];
        c = [0,1/6,1/3,1/2,2/11,2/3,6/7,0,1];
    case 'RK8' %Cooper-Verner eigth-order method (again there are many)
        s = sqrt(21);
        a = [0,0,0,0,0,0,0,0,0,0,0;1/2,0,0,0,0,0,0,0,0,0,0;1/4,1/4,0,0,0,0,0,0,0,0,0;1/7,(-7-3*s)/98,(21+5*s)/49,0,0,0,0,0,0,0,0;(11+s)/84,0,(18+4*s)/63,(21-s)/252,0,0,0,0,0,0,0;(5+s)/48,0,(9+s)/36,(-231+14*s)/360,(63-7*s)/80,0,0,0,0,0,0;(10-s)/42,0,(-432+92*s)/315,(633-145*s)/90,(-504+115*s)/70,(63-13*s)/35,0,0,0,0,0;1/14,0,0,0,(14-3*s)/126,(13-3*s)/63,1/9,0,0,0,0;1/32,0,0,0,(91-21*s)/576,11/72,(-385-75*s)/1152,(63+13*s)/128,0,0,0;1/14,0,0,0,1/9,(-733-147*s)/2205,(515+111*s)/504,(-51-11*s)/56,(132+28*s)/245,0,0;0,0,0,0,(-42+7*s)/18,(-18+28*s)/45,(-273-53*s)/72,(301+53*s)/72,(28-28*s)/45,(49-7*s)/18,0];
        b = [1/20,0,0,0,0,0,0,49/180,16/45,49/180,1/20];
        c = [0,1/2,1/2,(7+s)/14,(7+s)/14,1/2,(7-s)/14,(7-s)/14,1/2,(7+s)/14,1];
    case 'ImplicitRK1' %Implicit Euler
        
        u = zeros(length(u0),length(t));
        u(:,1) = u0;
        
        for n = 1:length(t)-1
           
            t_old = t(n);
            t_new = t(n+1);
            u_old = u(:,n);
            u_new = u_old + (t_new - t_old)*f(t_old,u_old);
            
            g = @(u_new,u_old,t_new,t_old,f)(u_new - u_old - (t_new - t_old)*f(t_new,u_new));
            
            
            u_next = fsolve ( @(u_new)g(u_new,u_old,t_new,t_old,f), u_new, optimoptions ( 'fsolve', 'Display', 'off' ));
            u(:,n+1) = u_next;
            
        end
        u = u';  %transpose solution matrix
        
        return
    case 'ExactLinear' %Exact Solution (only to be used for scalar linear ODE

        u = zeros(length(u0),length(t));
        u(:,1) = u0;
        
        lambda = f(1,u0)/u0;
        
        for n = 1:length(t)-1
            h = t(n+1)-t(n);            %length of time step
            u(:,n+1) = u(:,n)*exp(lambda*h);
        end
        u = u';
        return
    case 'ExactNonLinear' %Exact Solution (only to be used for nonlinear ODE f(u) = sqrt(u^2 + 2))

        u = zeros(length(u0),length(t));
        u(:,1) = u0;
                
        for n = 1:length(t)-1
            h = t(n+1)-t(n);            %length of time step
            u(:,n+1) = sqrt(2)*sinh(h + asinh(u(:,n)/sqrt(2)));
        end
        u = u';
        return
    case 'ExactLinearSystem' %Exact Solution (for linear system of ODEs)

        u = zeros(length(u0),length(t));
        u(:,1) = u0;       
        
        for n = 1:length(t)-1
            h = t(n+1)-t(n);            %length of time step
            u(:,n+1) = (u(:,n)'*expm(f(0,1)*h))';
        end
        u = u';
        return
    otherwise
        fprintf('Error: Please define the Runge Kutta method. \n')
end

%SOLVING THE EQUATIONS.

%Pre-define the solution matrix and specify initial conditions
u = zeros(length(u0),length(t));
u(:,1) = u0;

%Iterate over each time step explicitly
for n = 1:length(t)-1
    
    %This function carries out the iterative step of a general form
    % of the Runge-Kutta method with inputs: (time step, initial time,
    % intial condition, function, coefficient matrices).
    
    h = t(n+1)-t(n);            %length of time step
    dim = length(u0);           %dimension of the ODE problem
    S = length(b);              %order of the RK method (2nd/4th)
    k = zeros(dim,S);           %matrix for other k values
    k(:,1) = h*f(t(n),u(:,n));  %definition of k1
    
    %calculate the coefficients k
    for i = 2:S
        temp = zeros(dim,1);
        for j = 1:i-1
            temp = temp + a(i,j)*k(:,j);
        end
        k(:,i) = h*f(t(n) + c(i)*h,u(:,n) + temp);
    end
    
    %calculate the final solution  
    u(:,n+1) = u(:,n) + sum(b.*k,2);
    
end
u = u';  %transpose solution matrix
end



