%sfs problem
clear all;
close all;
function val = Ravg(A,i,j)
    val = ( A(i,j-1) + A(i-1,j) + A(i,j+1) + A(i+1,j) )/6 +(A(i-1,j-1)+ A(i-1,j+1)+A(i+1,j-1)+A(i+1,j+1))/12;
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

M=256;
N=256;

radius=100; % Radius of sphere

Depth = zeros(M,N);
E = 0.2 * ones(M,N); % Background Light


s_orig=[0,0]; % p,q cordinate of source

% Initialization of different matrices
p_init = zeros(M,N);
q_init = zeros(M,N);
f_init = zeros(M,N);
g_init = zeros(M,N);
p_orig = zeros(M,N);
q_orig = zeros(M,N);
f_orig = zeros(M,N);

mask = zeros(M,N);
boundary = zeros(M,N);
boundaryfg = zeros(M,N);
g_orig = zeros(M,N);

for i=1:M,
    for j=1:N,
        current_radius = sqrt((i-M/2)^2 + (j-N/2)^2);
        if(current_radius < radius) % Assign Depth only for points that are part of sphere in image
            Depth(i,j) = round(sqrt(radius^2 - (i-M/2)^2 - (j-N/2)^2)); % estimate depth for each pixel of sphere
            p = (i-M/2)/Depth(i,j);
            q = (j-N/2)/Depth(i,j);
            mask(i,j)=1;
            
            % Calculate E from p,q & s
            temp = Rval(p, q, s_orig);
            
            % closed curve values for p,q
            if (round(current_radius) == round(radius/2))
                p_init(i,j)= p;
                q_init(i,j)= q;
                boundary(i,j)=1;
            end
            p_orig(i,j) = p;
            q_orig(i,j) = q;
            % closed curve values for f,g
            if (round(current_radius) >= radius -1 && round(current_radius) < radius+1)
                f_init(i,j) = 2*p/(1+sqrt(1+p^2+q^2));
                g_init(i,j) = 2*q/(1+sqrt(1+p^2+q^2));
                boundaryfg(i,j)=1;

            end
            f_orig(i,j) = 2*p/(1+sqrt(1+p^2+q^2));
            g_orig(i,j) = 2*q/(1+sqrt(1+p^2+q^2));
            
            % Ensure nonegativity of E
            if(temp>0)
                E(i,j) = temp;
            else
                E(i,j) = 0;

            end
            
        end
    end
end

%Source passed to algorithm used for reconstruction
s = s_orig;%[1,1];

figure;
imshow(mat2gray(E));

% Noisy image with gaussian noise
E_noise = imnoise(E,'gaussian',0,5);

save('DataFile1.mat', 'E', 's','radius','mask','boundary', 'p_init', 'q_init', 'f_init', 'g_init', 'E_noise', 'Depth','boundaryfg');
