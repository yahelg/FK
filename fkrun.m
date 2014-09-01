function cell=fkrun()

%     f = fkfuns;
%     k = 0.02;
%     x = f.initx(20,k);
%     f.plotconfig(0,k,x);
%     disp(x);
%     return;

    % Use parrallel computing
    if matlabpool('size') == 0 
        matlabpool open
    end
    ions = 20;
    reps = 50;
    
    k = linspace(0.1,1,100);
    n = length(k);
%     cell = zeros(n,20);
%     a = zeros(1,n);
    
%     parfor i=1:n
%         V = 4*pi;
%         f = fkfuns;
%         beta = 1/100;
%         inx = f.initx(20,k(i));
%         x = f.fkgrad(V,k(i),beta,inx);
%         cell(i,:) = x;
%         a(i) = mean(diff(x));
%     end
%     plot(k,a,'.-');
%     disp('ok');
%     return;

    cell = zeros(n,reps,ions);
    parfor i=1:n
        V = 4*pi;
        f = fkfuns;
        inx = f.initx(20,k(i));
%         r = [];
%         en = [];
        
        c = zeros(reps,ions);
        for j=1:reps
            x = f.fk(V,k(i),inx);
%             r = [r;x];
%             en = [en f.energy(V,k(i),x)];
            c(j,:) = x;
        end
        cell(i,:,:) = c;
        
%         [b,c] = sort(en);
%         r1 = r(c(1),:);
%         r2 = r(c(2),:);
%         cell(i,:) = r1;
%         if norm(r1-r2) < 1e-4
%             a(i) = mean(diff(r1));
%         else
%             a(i) = 0;
%         end
%         
%         a(i) = fk(20,4*pi,k(i));
    end
%     plot(k,a,'.-');
%     disp('ok');
end