%% Parameter choiche:
% The fundamental parameters, to be used in the following section are
% defined.

N = 300;

%Integration time;
T = 220;

%Omelyan parameter lambda;
lam = 0.1931833275;


%Number of measurements;
cycle = 5e4;


% S(x(t)) parameters;

e = (1*(pi^2))/N;     
cc = pi/(2*e);          %coefficient multiplying the force;
den = 1/(2*e);

dpi = 2*pi;             %some useful definition;
oml = 1-2*lam;

%metadynamic parameters;

Qtrh = 20;              %Threshold value of the charge
hgt = 0.73;             %height of the time dependent potential inside Qtrh;
dq = 0.82;              %width of the  time dependent potential inside Qtrh;
kk = 2;                 %strength of the potential outside Qtrh;
upd = 20;               % # of sweeps after which the update of the 
                        %time dependent potential is performed;

metad = 0;              
q = -Qtrh-dq:dq:Qtrh+dq;    %charge 'lattice';
v = zeros(length(q),1);     %grid for the time dependent potential;


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
% The steps of the algorithm are performed as follows: 
%1) Symplettic integration with Omelyan integrator;
%2) computation of the final kinetic and potential energy;
%3) Metropolis step;
%4) Update of the time dependent potential;

for k = 1:cycle

    LP0 = -(K(k) + V);                %logarithm of the probability
                                      %associated with the last path;
    
  
    dt = (rand()-1)*0.005 + 0.025;    %time step of the integrator;

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
        elseif index <= 0
            metad = -2*kk*(-Q(k)+Qtrh);
        else

                %coefficient for the derivative
                %of the the time dependent potential;
                metad = (v(index) - v(index+1))/dq; 

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
         elseif index <= 0
            metad = -2*kk*(-Q(k)+Qtrh);
         else
                metad = (v(index) - v(index+1))/dq;

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
         elseif index <= 0
            metad = -2*kk*(-Q(k)+Qtrh);
         else
                metad = (v(index) - v(index+1))/dq;

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
             elseif index <= 0
                 metad = -2*kk*(-Q(k)+Qtrh);
             else
                 metad = (v(index) - v(index+1))/dq;

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
             elseif index <= 0
                metad = -2*kk*(-Q(k)+Qtrh);
             else

                metad = (v(index) - v(index+1))/dq;

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
        elseif index <= 0
            metad = -2*kk*(-Q(k)+Qtrh);
        else
                metad = (v(index) - v(index+1))/dq;

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
         elseif index <= 0
            metad = -2*kk*(-Q(k)+Qtrh);
         else
                metad = (v(index) - v(index+1))/dq;

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
    elseif index <= 0        %Update id Q < - Qtrh;
      Vend = den*sum((sin(pi*d).^2)) + kk*(-Q(k) -Qtrh)^2; 
    else                     %Update inside the barrier;
      Vend =  den*sum((sin(pi*d).^2))+v(index)+(v(index+1)-v(index))*(Q(k)-q(index));
    end
    
    %Final Kinetic energy;
    Kend = sum((p(k,1:end).^2))/2;
    
    
    %Logarithm of the probability associated to the new path;
    LP1 = -Vend-Kend;
    
    %Metropolis step;
    
    if log(rand()) < LP1 - LP0
        %Update of the stored path, charge and distances;
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
        elseif index <= 0
            metad = -2*kk*(-Q(k)+Qtrh);
        else
            V = V-v(index) - (v(index+1)-v(index))*(Q(k)-q(index));
            v(index) = v(index) + hgt*(1 - (Q(k) - q(index))/dq);
            v(index+1) = v(index+1) + hgt*((Q(k) - q(index))/dq);
            metad = -(-v(index) + v(index+1))/dq;
            V = V +v(index)+(v(index+1)-v(index))*(Q(k)-q(index));

        end
    end
         
end

a_rate = a_rate/cycle;

filename1 =('Charge.txt');
filename2 = ('Path.txt');
filename3 = ('Parameters.txt');


%Paths
dlmwrite(filename1,Q,'-append','delimiter','\t');
dlmwrite(filename2,y,'-append','delimiter','\t');
dlmwrite(filename3,a_rate,'-append','delimiter','\t');

