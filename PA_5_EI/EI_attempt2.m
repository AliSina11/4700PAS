clear all;

nx = 50;
ny = 50;
V = zeros(nx,ny);
G = sparse(nx*ny, nx*ny);

Inclusion = 0;

% (Psi_(i+1, j) + 2*Psi_(i, j) + Psi_(i-1, j))/(dx)^2 +
%    + (Psi_(i, j+1) + 2*Psi_(i, j) + Psi_(i, j-1))/(dy)^2

% Finite Difference
%   Say we index G(P,Q):
%   Then P refers to the current row of the G matrix, which is multiplied
%   by V to get D. Q refers to the exact point in physical space. So for
%   any point, G(P,Q) is the center, G(P,Q+1) is the next (right) point,
%   G(P,Q-1) is the last (left) point, G(P,Q + nx) is a point that is a row
%   above (physically above, or next in y), G(P,Q - nx) is a point that is a row
%   above (physically below, or previous in y).

for i = 0:ny
    if (0 < i) && (i < ny - 1)
        for j = 0:nx
            if(0 < j) && (j < nx - 1)
                % center 
                G(i*ny + j + 1, i*ny + j + 1) = -4;
                % right
                G(i*ny + j + 1, i*ny + j + 1 + 1) = 1;
                % left
                G(i*ny + j + 1, i*ny + j + 1 - 1) = 1;
                % up
                G(i*ny + j + 1, i*ny + j + 1 + nx) = 1;
                % down
                G(i*ny + j + 1, i*ny + j + 1 - nx) = 1;
            end
        end
    end
end

% Inclusion
for i = 10:20
    if (0 < i) && (i < ny - 1)
        for j = 10:20
            if(0 < j) && (j < nx - 1)
                % center 
                G(i*ny + j + 1, i*ny + j + 1) = -2;
                % right
                G(i*ny + j + 1, i*ny + j + 1 + 1) = 1;
                % left
                G(i*ny + j + 1, i*ny + j + 1 - 1) = 1;
                % up
                G(i*ny + j + 1, i*ny + j + 1 + nx) = 1;
                % down
                G(i*ny + j + 1, i*ny + j + 1 - nx) = 1;
            end
        end
    end
end

% Background conditions
for i = 0:ny-1
    if i == 0 || i == ny - 1
        for j = 0:nx-1
            G(ny * i + j + 1, ny * i + j + 1) = 1;
        end
    else
        G(ny * i + 1, ny * i + 1) = 1;
        G(ny * i + nx, ny * i + nx) = 1;
    end
end

figure('name', 'G Matrix')
spy(G)

nmodes = 20;
[E, D] = eigs(G, nmodes, 'SM');

figure('name', 'Eigenvalues')
plot(diag(D), '*');

np = ceil(sqrt(nmodes))
figure('name', 'Modes')
for k = 1:nmodes
    M = E(:,k);
    for i = 1:nx
        for j = 1:ny
            n = i + (j-1)*nx;
            V(i,j) = M(n);
        end
        subplot(np,np,k), surf(V, 'Linestyle', 'none')
        title(['EV = ' num2str(D(k,k))])
    end
end