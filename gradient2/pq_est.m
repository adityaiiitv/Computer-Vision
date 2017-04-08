% p,q ESTIMATE FROM A SINGLE SOURCE
clear all;
load('ASSIGN2.mat');

function val = Rm(A,i,j)
    val = (A(i,j-1) + A(i-1,j) + A(i,j+1) + A(i+1,j))/6 + (A(i-1,j-1) + A(i-1,j+1) + A(i+1,j-1) + A(i+1,j+1))/12;
end

function val = Rf(f,g,fs,gs)
    val = (-2*f*((4-fs^2-gs^2)*8 + 16*f*fs + 16*g*gs) + 16*fs*(4+f^2+g^2)) / ((4+f^2+g^2)^2*(4+fs^2+gs^2));
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

function val=Rv(p,q,s)
    val = (s(1) * p+s(2)*q+1) / sqrt((s(2) * s(2)+s(1) * s(1)+1) * (p*p+q*q+1));
end

M=size(E,1);
N=size(E,2);

En=E;

% p_init, q_init HAVE INITIAL CONDITION
p_o=p_init;
q_o=q_init;

pn=zeros(size(En));
qn=zeros(size(En));

% ITERATIONS AND LAMBDA
iter=10;
lambda=0.5 ;

% WE USE ITERATIVE METHOD TO GET THE VALUES OF p,q
for kk=1:iter,
    disp(kk)
    for i=2:(M-1),
        for j=2:(N-1),
            if(boundary(i,j)==0&&mask(i,j)==1)
                % IF POINTS ON BOUNDARY, THEN COMPUTE p,q
                pn(i,j) = Rm(p_o,i,j) + (1/lambda)*( En(i,j) - Rv(p_o(i,j), q_o(i,j), s)) * Rp(p_o(i,j),q_o(i,j),s(1),s(2));
                qn(i,j) = Rm(q_o,i,j) + (1/lambda)*( En(i,j) - Rv(p_o(i,j), q_o(i,j), s)) * Rq(p_o(i,j),q_o(i,j),s(1),s(2));
            else 
                % IF ALREADY ON BOUNDARY, THEN CORRECT VALUE 
                pn(i,j)=p_o(i,j);
                qn(i,j)=q_o(i,j);
            end
        end
    end
    % VALUES ARE RE-ASSIGNED FOR THE NEXT ITERATION
    p_o=pn;
    q_o=qn;
end
imshow(pn);
imshow(qn);
