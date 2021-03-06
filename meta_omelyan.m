%% Parameter choiche:
% The fundamental parameters, to be used in the following section are
% defined.

N = 600;

%Integration time;
T = 312;

%Omelyan parameter lambda;
lam = 0.1931833275;


%Number of measurements;
cycle = 2e4;


% S(x(t)) parameters;

e = (1*(pi^2))/N;     
cc = pi/(2*e);          %coefficient multiplying the force;
den = 1/(2*e);

dpi = 2*pi;             %some useful definition;
oml = 1-2*lam;

%metadynamic parameters;

Qtrh = 20;              %Threshold value of the charge
hgt = 1.03;             %height of the time dependent potential inside Qtrh;
dq = 1;              %width of the  time dependent potential inside Qtrh;
kk = 0.70;                 %strength of the potential outside Qtrh;
upd = 20;               % # of sweeps after which the update of the 
                        %time dependent potential is performed;

metad = 0;              
q = -Qtrh-dq:dq:Qtrh+dq;         %charge 'lattice';
td_pot = zeros(length(q),1);     %grid for the time dependent potential;
store = zeros(length(q),1);

%Lattice grid;
nu = linspace(2,N+1,N);
nd = linspace(0,N-1,N);
nu(N) = 1;
nd(1) = N;

%AcceptanceRate;
a_rate = 0;

%Path and distance between path points;
y = zeros(N,1);
y0 = zeros(N,1);
d = zeros(N,1);


%Useful vector for the minimization of the path;
z = zeros(N,1);
z_step = zeros(N,1);
Q_min = zeros(N,1);
d_min = zeros(N,1);
d0 = zeros(N,1);


%Strength vector;
Force = zeros(N,1);
%Charge vector;
Q = zeros(cycle,1);
Q0 = 0;


%Momenta taken with a normal distribution with mean 0 and variance 1;
p = normrnd(0,1,[cycle+1,N]);



%Initial and final kinetic values;
Kend = 0;
K = zeros(cycle+1,1);


for i = 1:cycle+1
    for n = 1:N
    K(i) = K(i) + (p(i,n)^2)/2;
    end
end

V = 0;

%% Main body:
% % The steps of the algorithm are performed as follows: 
% %1) Symplettic integration with Omelyan integrator;
% %2) computation of the final kinetic and potential energy;
% %3) Metropolis step;
% %4) Update of the time dependent potential;

for k = 1:cycle

    LP0 = -(K(k) + V);                %logarithm of the probability
                                      %associated with the last path;
    
  
    dt = (2*rand()-1)*0.005 + 0.025;    %time step of the integrator;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Omelyan integrator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %first step
        
        for n = 1:N
            y(n) = y(n) + lam*p(k,n)*dt;    %New path point;
        end
        
        for n = 1:N
            d(n) = y(nu(n))-y(n);           %new distance between points;
        end
        
        Q(k) =  sum(sin(dpi*d))/dpi;  %charge for the configuration;
        index = floor((Q(k)-q(1))/dq + 1.); %position of the charge on the lattice;
        
        if index > length(q)-1  
            metad = -2*kk*(Q(k)- Qtrh);
        elseif index <= 1
            metad = -2*kk*(Q(k)+Qtrh-dq);
        else

                %coefficient for the derivative
                %of the the time dependent potential;
                metad = (td_pot(index) - td_pot(index+1))/dq; 

        end
            
        
       
        for n = 1:N
            % time dependent force;
            f = metad*(-cos(dpi*d(n)) + cos(dpi*d(nd(n))));
            %total force;
            Force(n) =  -cc*(sin(dpi*d(nd(n)))-sin(dpi*d(n)))+ f;
        end
    
        for n = 1:N
            %Pn(i+1/2);
            p(k,n) = p(k,n) + (dt/2)*Force(n);
            %Yn(i+1);
            y(n) = y(n) + p(k,n)*dt*oml;
        end
        
         for n = 1:N
            d(n) = y(nu(n))-y(n);
            
         end
        
         Q(k) =  sum(sin(dpi*(d)))/dpi;
         index = floor((Q(k)-q(1))/dq + 1.);
         
         if index > length(q)-1  
            metad = -2*kk*(Q(k)- Qtrh);
         elseif index <= 1
            metad = -2*kk*(Q(k)+Qtrh-dq);
         else
                metad = (td_pot(index) - td_pot(index+1))/dq;

         end
            
         
         for n = 1:N
            f = metad*(-cos(dpi*d(n)) + cos(dpi*d(nd(n))));
            Force(n) =  -cc*(sin(dpi*d(nd(n)))-sin(dpi*d(n)))+ f; 
            p(k,n) = p(k,n) + (dt/2)*Force(n);
         end
         
         for n = 1:N
             y(n) = y(n) + 2*lam*dt*p(k,n);
         end
         
         for n = 1:N
            d(n) = y(nu(n))-y(n);
            
         end
         
         Q(k) =  sum(sin(dpi*d))/dpi;
         index = floor((Q(k)-q(1))/dq + 1.);
         if index > length(q)-1  
            metad = -2*kk*(Q(k)- Qtrh);
         elseif index <= 1
            metad = -2*kk*(Q(k)+Qtrh-dq);
         else
                metad = (td_pot(index) - td_pot(index+1))/dq;

         end
            
     
         %Intermediate steps
         for i = 2:T-1
         
            
             for n = 1:N
                
                f = metad*(-cos(dpi*d(n)) + cos(dpi*d(nd(n))));
                Force(n) =  -cc*(sin(dpi*d(nd(n)))-sin(dpi*d(n)))+ f;
                p(k,n) = p(k,n) + (dt/2)*Force(n);
             end

             for n = 1:N
                y(n) = y(n) + dt*oml*p(k,n); 
             end

             for n = 1:N
                d(n) = y(nu(n))-y(n);

             end

             Q(k) =  sum(sin(dpi*d))/dpi;
             index = floor((Q(k)-q(1))/dq + 1.);
             if index > length(q)-1  
                 metad = -2*kk*(Q(k)- Qtrh);
             elseif index <= 1
                 metad = -2*kk*(Q(k)+Qtrh-dq);
             else
                 metad = (td_pot(index) - td_pot(index+1))/dq;

             end


             for n = 1:N
              
                f = metad*(-cos(dpi*d(n)) + cos(dpi*d(nd(n))));
                Force(n) =  -cc*(sin(dpi*d(nd(n)))-sin(dpi*d(n)))+ f; 
                p(k,n) = p(k,n) + (dt/2)*Force(n);
             end

             for n = 1:N
                 y(n) = y(n) + 2*lam*dt*p(k,n);
             end

             for n = 1:N
                d(n) = y(nu(n))-y(n);

             end

             Q(k) =  sum(sin(dpi*d))/dpi;
             index = floor((Q(k)-q(1))/dq + 1);
             
             if index > length(q)-1  
                metad = -2*kk*(Q(k)- Qtrh);
             elseif index <= 1
                metad = -2*kk*(Q(k)+Qtrh-dq);
             else

                metad = (td_pot(index) - td_pot(index+1))/dq;

             end
            
         
        end
    
    %%%%%%%%%%%%%%%%
    %Last step

        for n = 1:N
           f = metad*(-cos(dpi*d(n)) + cos(dpi*d(nd(n))));
           Force(n) =  -cc*(sin(dpi*d(nd(n)))-sin(dpi*d(n)))+ f;
           p(k,n) = p(k,n) + (dt/2)*Force(n);
        end
        
        for n = 1:N
            y(n) = y(n) + dt*oml*p(k,n); 
        end
        
        for n = 1:N
            d(n) = y(nu(n))-y(n);           
        end
         
        Q(k) =  sum(sin(dpi*d))/dpi;
        index = floor((Q(k)-q(1))/dq + 1.);
        if index > length(q)-1  
            metad = -2*kk*(Q(k)- Qtrh);
        elseif index <= 1
            metad = -2*kk*(Q(k)+Qtrh-dq);
        else
                metad = (td_pot(index) - td_pot(index+1))/dq;

        end
            
         
        for n = 1:N
            f = metad*(-cos(dpi*d(n)) + cos(dpi*d(nd(n))));
            Force(n) =  -cc*(sin(dpi*d(nd(n)))-sin(dpi*d(n)))+ f; 
            p(k,n) = p(k,n) + (dt/2)*Force(n);
        end
         
        for n = 1:N
             y(n) = y(n) + lam*dt*p(k,n);
        end
        
         Q(k) =  sum(sin(dpi*d))/dpi;
         index = floor((Q(k)-q(1))/dq + 1.);
         if index > length(q)-1  
            metad = -2*kk*(Q(k)- Qtrh);
         elseif index <= 1
            metad = -2*kk*(Q(k)+Qtrh-dq);
         else
                metad = (td_pot(index) - td_pot(index+1))/dq;

         end
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metropolis step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    y = mod(y,1);
    
    for n = 1:N
        d(n) = y(nu(n))-y(n);
    end
    
    
    %Final potential energy;
    if  index > length(q)-1  %Update if Q > Qtrh;
     
      Vend = den*sum((sin(pi*d).^2)) + kk*(Q(k) -Qtrh)^2; 
    elseif index <= 1        %Update id Q < - Qtrh;
      Vend = den*sum((sin(pi*d).^2)) + kk*(Q(k) +Qtrh-dq)^2; 
    else                     %Update inside the barrier;
      Vend =  den*sum((sin(pi*d).^2))+td_pot(index)+(td_pot(index+1)-td_pot(index))*(Q(k)-q(index))/dq;
    end
    
    %Final Kinetic energy;
    Kend = sum((p(k,1:end).^2))/2;
    
    
    %Logarithm of the probability associated to the new path;
    LP1 = -Vend-Kend;
    
    %Metropolis step;
    
    if log(rand()) < LP1 - LP0
        %Update the stored path, charge and distances;
        y0 = y;
        d0 = d;
        Q0 = Q(k);
        %Update of the potential and the acceptance rate;
        a_rate = a_rate +1;      
        V = Vend;
        
    end
    
    y = y0;
    d = d0;    
    Q(k) =  Q0;
    
    index = floor((Q(k)-q(1))/dq + 1.);
    %Update of the time dependet potential every 'upd' steps
    if mod(k,upd) == 0
        if index > length(q)-1  
            metad = -2*kk*(Q(k)- Qtrh);
        elseif index <= 1
            metad = -2*kk*(Q(k)+Qtrh-dq);
        else
            %subtract the old value of the t-d-potential from the total
            %try updating also the first point
            %potential;
            V = V-td_pot(index)-(td_pot(index+1)-td_pot(index))*(Q(k)-q(index))/dq;
            %Update V(i);
            td_pot(index) = td_pot(index) + hgt*(1 - (Q(k) - q(index))/dq);           
            %Update V(i+1);
            td_pot(index+1) = td_pot(index+1) + hgt*((Q(k) - q(index))/dq);
            %Update coefficient of metadynamic;
            metad = (td_pot(index) - td_pot(index+1))/dq;
            %Update the last potential with the new value of the
            %t-d-potential;
            V = V+td_pot(index)+(td_pot(index+1)-td_pot(index))*(Q(k)-q(index))/dq;

        end
    end
    
    if k > 7e4 && mod(k,upd) == 0
        store = store + td_pot;
    end
    %minimization of the action;
    z = y;
    for j = 1:2
        for l = 1:N
            
            if abs(z(nu(l))-z(nd(l)))< 0.5 && abs(z(l)-z(nd(l))) < 0.5 && abs(z(nu(l))-z(l))<0.5
                
                z_step(l) = (z(nd(l)) + z(nu(l)))/2;
            else
                z_step(l) = z(l);
            end
            
            
        end
        z = z_step;
        
    end
    
    for l = 1:N
        d_min(l) = z(nu(l))-z(l);           
    end
    %charge for the minimized action;
    Q_min(k) = sum(sin(dpi*d_min))/dpi;
         
end

a_rate = a_rate/cycle;
store = store/1500;

%Write charge on file;
writematrix(Q,'Charge_300_0.5(3).txt','delimiter','tab');
%write final path on file;
writematrix(y,'path_300_0.5(3).txt','delimiter','tab');
%write acceptance rate on file;
writematrix(a_rate,'acceptance_rate_300_0.5(3).txt','delimiter','tab');
%write time dependet potential on file;
writematrix(td_pot,'time_d_potential_300_0.5(3).txt','delimiter','tab');

