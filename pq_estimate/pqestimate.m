% p,q estimate from single source

clear all;

load('DataFile1.mat');

function val = Ravg(A,i,j)
    val = ( A(i,j-1) + A(i-1,j) + A(i,j+1) + A(i+1,j) )/6 +(A(i-1,j-1)+ A(i-1,j+1)+A(i+1,j-1)+A(i+1,j+1))/12;
end

function val = Rf(f,g,fs,gs)
    val = (-2*f*( (4 - fs^2 - gs^2)*8 + 16*f*fs + 16*g*gs) + 16*fs*(4 + f^2 + g^2)) / ((4 + f^2 + g^2)^2*(4 + fs^2 + gs^2));
end

function val = Rfgval(f,g,s)
    val = ( (4 - f^2 - g^2)*(4 - s(1)^2 - s(2)^2) + 16*f*s(1) + 16*g*s(2) )/( (4 + f^2 + g^2) * (4 + s(1)^2 + s(2)^2) );
end

function val = Rg(f,g,fs,gs)
    val = (-2*g*((4 - fs^2 - gs^2)*8 + 16*f*fs + 16*g*gs) + 16*gs*(4 + f^2 + g^2))/((4 + f^2 + g^2)^2 * (4 + fs^2 + gs^2));
end

function val = Rp(p,q,ps,qs)
    val = (ps*(q*q +1) - p*(qs*q + 1))/(sqrt((qs*qs + ps*ps +1) * (p*p + q*q +1)) * (p*p+q*q+1));
end

function val = Rq(p,q,ps,qs)
    val = (qs*(p^2 +1) - q*(ps*p + 1))/(sqrt((qs^2 + ps^2 +1) * (p^2 + q^2 +1)) * (p^2+q^2+1));
end

function val = Rval(p,q,s)
    val = (s(1)*p + s(2)*q + 1)/sqrt( (s(2)*s(2) + s(1)*s(1) +1) * (p*p + q*q +1));
end

M = size(E,1);
N = size(E,2);

En = E;

% p_init, q_init contain the initial conditions (values on the boundary)
p_old = p_init;
q_old = q_init;

pn = zeros(size(En));
qn = zeros(size(En));

% Set the number of iterations, value of the regularization (lambda)
iters = 100;
lambda = 0.5 ;

% We use the iterative method derived from variational calculus to estimate the values of p,q
for kk = 1:iters,
    disp(kk)
    for i=2:(M-1),
        for j=2:(N-1),
            if(boundary(i,j)==0 && mask(i,j) ==1)
                % If the points lies not on the boundary, only then compute the value of p,q
                pn(i,j) = Ravg(p_old,i,j) + (1/lambda)*( En(i,j) - Rval(p_old(i,j), q_old(i,j), s))*Rp(p_old(i,j),q_old(i,j),s(1),s(2));
                qn(i,j) = Ravg(q_old,i,j) + (1/lambda)*( En(i,j) - Rval(p_old(i,j), q_old(i,j), s))*Rq(p_old(i,j),q_old(i,j),s(1),s(2));
            else 
                % Otherwise, since it is already on the boundary, we already have the correct value 
                pn(i,j) = p_old(i,j);
                qn(i,j) = q_old(i,j);
            end
        end
    end
    % Re assign the values for next iteration
    p_old = pn;
    q_old = qn;
end

imshow(pn);
imshow(qn);
