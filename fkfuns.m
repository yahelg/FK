function f=fkfuns()
    f.fk = @fk;
    f.force = @force;
    f.plotconfig = @plotconfig;
    f.energy = @energy;
    f.spect = @spect;
    f.initx = @initx;
    f.fkgrad = @fkgrad;
    f.energies = @energies;
    f.plotEnergies = @plotEnergies;
    f.plotw = @plotw;
    f.plotground = @plotground;
    f.plotconfigfirst = @plotconfigfirst;
end

function x=fk(V,k,inx)
    n = length(inx);
    w = mean(diff(inx));
    r = random('unif',-w/10,w/10,1,n);
    y = inx + r;
%     if nargin < 3
%         V = pi;
%         k = 0.02;
%         n = 10;
%         f = fk(n,V,k);
%     end

%     R = 2 * (1/(k*n^2))^(1/3);
%     
%     y = sort(random('unif',-R*n/2,R*n/2,n,1));

%     x = [-0.2  -0.1   0.0   0.1 0.2 0.4 0.6 0.8 0.9];
%     n = length(x);
    options = odeset('RelTol',1e-9,'AbsTol',1e-9,'Events',@(t,x)...
     eventfun(V,k,t,x));
    [T,Y] = ode45(@(t,x) force(V,k,t,x),[0 Inf],y,options);
%      disp(Y(end,:));
%     plot(T,Y(:,1));
    x = Y(end,:);
%     disp(energy(V,k,x));
%     plotconfig(V,k,x);
%     w = mean(diff(x));
%     w = min(spect(V,k,x));
%     dx = force(V,k,0,x);
%     dx
end

function x=fkgrad(V,k,beta,inx)
    n = length(inx);
    y = inx;
    options = odeset('RelTol',1e-9,'AbsTol',1e-9,'Events',@(t,x)...
     eventfunmin(V,k,beta,t,x));
    [T,Y] = ode45(@(t,x) gradualforce(V,k,beta,t,x),[0 Inf],y,options);
    disp(max(T));
    x = Y(end,:);
end


function x=initx(n,k)
    R = 2 * (1/(k*n^2))^(1/3);
    
    y = sort(random('unif',-R*n/2,R*n/2,n,1));

%     x = [-0.2  -0.1   0.0   0.1 0.2 0.4 0.6 0.8 0.9];
%     n = length(x);
    options = odeset('RelTol',1e-9,'AbsTol',1e-9,'Events',@(t,x)...
     eventfun(0,k,t,x));
    [T,Y] = ode45(@(t,x) force(0,k,t,x),[0 Inf],y,options);
    x = Y(end,:);
end

function dx = force(V,k,t,x)
    n = length(x);

    dx = -k^2*x - V/(2*pi)*sin(2*pi*x);
    for i=1:n
        dx(i) = dx(i) - sum(1/2./(x(i+1:end) - x(i)).^2);
        dx(i) = dx(i) + sum(1/2./(x(1:i-1) - x(i)).^2);
    end
end

function dx = gradualforce(V,k,beta,t,x)
    n = length(x);

    Vt = V*t*beta;
    if Vt>V
        Vt = V;
    end
    
    dx = force(Vt,k,t,x);
end

function [x,isterm,dir] = eventfun(V,k,t,x)
    dx = force(V,k,0,x);
    x = norm(dx) - 1e-7;
    isterm = 1;
    dir = 0;  %or -1, doesn't matter
end

function [x,isterm,dir] = eventfunmin(V,k,beta,t,x)
    dx = force(V,k,0,x);
    tmin = 1/beta;
    x = norm(dx) - 1e-7 + (t<tmin);
    isterm = 1;
    dir = 0;  %or -1, doesn't matter
end

function plotconfig(V,k,x)
    s = min(x)-1/2; f = max(x)+1/2;
    seg = linspace(s,f,1000);
%     plot(seg,-V/(2*pi)^2*cos(2*pi*seg),seg,...
%         1/2*k^2*seg.^2,x,-V/(2*pi)^2*cos(2*pi*x),'*');
    plot(seg,-V/(2*pi)^2*cos(2*pi*seg)+1/2*k^2*seg.^2,...
        x,-V/(2*pi)^2*cos(2*pi*x)+1/2*k^2*x.^2,'*');
end

function E=energy(V,k,x)
    n = length(x);
    E = sum(1/2*k^2*x.^2 - V/(2*pi)^2 * cos(2*pi*x));
    for i=1:n
        E = E + sum(1/2./abs(x(i+1:end) - x(i)));
        E = E + sum(1/2./abs((x(1:i-1) - x(i))));
    end
end

function En=energies(V,k,cell)
    n = length(cell);
    En = zeros(1,n);
    for i=1:n
        En(i) = energy(V,k(i),cell(i,:));
    end
end

function plotEnergies(V,k,cell)
    [n,reps,~] = size(cell);
    En = zeros(n,10);
    for i=1:n
        e = zeros(reps,1);
        for j=1:reps
            e(j) = energy(V,k(i),cell(i,j,:));
        end
        d = sort(e);
        En(i,:) = d(1:10);
    end
    plot(k,En);
end

function plotw(V,k,cell)
    [n,reps,~] = size(cell);
    w = zeros(n,1);
    for i=1:n
        e = zeros(reps,1);
        for j=1:reps
            e(j) = energy(V,k(i),cell(i,j,:));
        end
        [~,f] = sort(e);
        x = cell(i,f(1),:);
        w(i) = mean(diff(x));
    end
    plot(k,w);
end

function plotground(V,k,cell)
    [n,reps,ions] = size(cell);
    w = zeros(n,1);
    kd = [];
    wd = [];
    for i=1:n
        e = zeros(reps,1);
        for j=1:reps
            e(j) = energy(V,k(i),cell(i,j,:));
        end
        [~,f] = sort(e);
        x = squeeze(cell(i,f(1),:));
        y = squeeze(cell(i,f(2),:));
        w(i) = mean(diff(x));
        if norm(x-y)>1e-4
            kd = [kd k(i)];
            wd = [wd w(i)];
        end
    end
    plot(k,w,'.-',kd,wd,'*');
    xlabel('k');
    ylabel('w');
    title(sprintf('w vs k for %d ions, %d repetitions V=%6.4g',...
        ions,reps,V));
end

function plotconfigfirst(V,k,cell)
    [n,reps,~] = size(cell);
    w = zeros(n,1);
    for i=1:n
        e = zeros(reps,1);
        for j=1:reps
            e(j) = energy(V,k(i),cell(i,j,:));
        end
        [~,f] = sort(e);
        x = squeeze(cell(i,f(1),:));
        y = squeeze(cell(i,f(2),:));
        if norm(x-y)<1e-4
            w(i) = mean(diff(x));
        else
            figure;
            subplot(2,1,1);
            plotconfig(V,k(i),x);
            title(sprintf('w=%6.4g, k=%6.4g, V=%6.4g, Energy=%6.4g',...
                mean(diff(x)),k(i),V,e(f(1))));
            subplot(2,1,2);
            plotconfig(V,k(i),y);
            title(sprintf('w=%6.4g, k=%6.4g, V=%6.4g, Energy=%6.4g',...
                mean(diff(y)),k(i),V,e(f(2))));
        end
    end
end

function spec=spect(V,k,x)
    n = length(x);
    spec = zeros(n);
    for i=1:n
        y = x;
        y(i) = [];
        spec(i,i) = -k^2-V*cos(2*pi*x(i))-sum(1./abs(y-x(i)).^3);
        for j=1:n
            if i~=j
                spec(i,j)=1/abs(x(i)-x(j))^3;
            end
        end
    end
    spec = -spec;
    spec = eig(spec);
end