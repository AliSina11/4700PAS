% This example shows how to calculate and plot both the
% fundamental TE and TM eigenmodes of an example 3-layer ridge
% waveguide using the full-vector eigenmode solver.  

% Refractive indices:
n1 = 3.34;          % Lower cladding
n2 = 3.305;          % Core % ridge = 3.305 to *3.44* in 10 steps
n3 = 1.00;          % Upper cladding (air)

% Layer heights:
h1 = 2.0;           % Lower cladding
h2 = 1.3;           % Core thickness
h3 = 0.5;           % Upper cladding

% Horizontal dimensions:
rh = 1.1;           % Ridge height
rw = 1.0;           % Ridge half-width
side = 1.5;         % Space on side

% Grid size:
dx = 0.0125;        % grid size (horizontal)
dy = 0.0125;        % grid size (vertical)

lambda = 1.55;      % vacuum wavelength
nmodes = 1;         % number of modes to compute

rwArray = zeros(10,1);
nEffArray = zeros(10,1);
n2Array = zeros(10,1);
i = 1;
for n2 = 3.305:0.0135:3.44
    [x,y,xc,yc,nx,ny,eps,edges] = waveguidemesh([n1,n2,n3],[h1,h2,h3], ...
                                            rh,rw,side,dx,dy); 

    % First consider the fundamental TE mode:

    
    [Hx,Hy,neff] = wgmodes(lambda,n2,nmodes,dx,dy,eps,'000A');

    fprintf(1,'neff = %.6f\n',neff);
    

    figure(1);
    subplot(121);
    contourmode(x,y,Hx(:,:,1)); % contourmode edit
    title('Hx (TE mode)'); xlabel('x'); ylabel('y'); 
    for v = edges, line(v{:}); end

    subplot(122);
    contourmode(x,y,Hy(:,:,1)); % contourmode edit
    title('Hy (TE mode)'); xlabel('x'); ylabel('y'); 
    for v = edges, line(v{:}); end
    
    n2Array(i) = n2;
    %rwArray(i) = rw;
    nEffArray(i) = neff;
    i = i + 1;
end
figure(2);
plot(n2Array, nEffArray);
xlabel('ridge index n2');
ylabel('effective index');
    % Next consider the fundamental TM mode
    % (same calculation, but with opposite symmetry)

    %[Hx,Hy,neff] = wgmodes(lambda,n2,nmodes,dx,dy,eps,'000S');

    %fprintf(1,'neff = %.6f\n',neff);

    %figure(2);
    %subplot(121);
    %contourmode(x,y,Hx(:,:,i));
    %title('Hx (TM mode)'); xlabel('x'); ylabel('y'); 
    %for v = edges, line(v{:}); end

    %subplot(122);
    %contourmode(x,y,Hy(:,:,i));
    %title('Hy (TM mode)'); xlabel('x'); ylabel('y'); 
    %for v = edges, line(v{:}); end
