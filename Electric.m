function [varargout] = Electric(varargin)
% ELECTRIC outputs a variety of values based on the input information. If 
% the object type is circle, Electric will find the electric flux through 
% the surface with the specified electric field, radius, and incline angle.
% If the object type is a rectangle, the electric flux through the surface 
% will be found with the input electric field, width, length, and incline 
% angle. If the object type is a sphere, the electric flux and field will
% be output from the input enclosed charge and radius. For all cases, a
% visualization will be displayed. All input values must be in SI units and 
% theta must be in degrees.

if nargin==0
    error('Error: Must input specified values.')
else 
    obj = varargin{1};
    switch obj
    case 'circle' % For circle, 3 additional inputs must be entered (E,r,theta) 
        E = varargin{2};
        r = varargin{3};
        theta = varargin{4};
        % Graphing Circle
        crcl(r,theta);
        % Calculating flux
        [flux] = circlecalc(E,r,theta);
        varargout{1} = flux;
    case 'rectangle' % For rectangle, 4 additional inputs must be entered (E,w,l,theta) 
        E = varargin{2};
        w = varargin{3};
        l = varargin{4};
        theta = varargin{5};
        % Calculating flux
        [flux] = rectanglecalc(E,w,l,theta);
        varargout{1} = flux;
        % Graphing Rectangle
        rct(w,l,theta);
    case 'sphere' % For sphere, 2 additional inputs must be entered (qenc,r) 
        qenc = varargin{2};
        r = varargin{3};
        % Calculating flux and field
        [flux,E] = spherecalc(qenc,r);
        varargout{1} = flux;
        varargout{2} = E;
        % Graphing the Sphere
        sphr(r,qenc);
    case 'concentric spheres' % For concentric spheres, must input 4 additional values (q1,r1,q2,r2)
        q1 = varargin{2};
        r1 = varargin{3};
        q2 = varargin{4};
        r2 = varargin{5};
        % Calculating flux and field at each surface
        [flux1,E1,flux2,E2] = concentric(q1,r1,q2,r2);
        varargout{1} = flux1;
        varargout{2} = E1;
        varargout{3} = flux2;
        varargout{4} = E2;
        if r2>r1
            % Graphing Concentric Spheres
            sphr(r1,q1);
            hold on
            q2outer = q1+q2;
            [s] = sphr(r2,q2outer);
            alpha(s,.5);
            hold off
        else
            error('Error: r2 must be greater than r1'); % Maybe try choosing which ever value is larger here.
        end
    end
end
function [flux] = circlecalc(E,r,theta) 
    % circle outputs electric flux from the given parameters
    flux = E*(pi*(r^2))*cosd(theta);
end 
function [flux] = rectanglecalc(E,w,l,theta) 
    % rectangle outputs electric flux from the given parameters
    flux = E*w*l*(cosd(theta));
end 
function [flux,E] = spherecalc(qenc,r) 
    % sphere outputs electric flux and field from the given parameters
    e0 = 8.85*(10^(-12));  % permitivity of free space 
    flux = qenc/e0;
    E = flux/(4*pi*(r^2));
end
function [flux1,E1,flux2,E2] = concentric(q1,r1,q2,r2) 
    e0 = 8.85*(10^(-12));
    flux1 = q1/e0;             % Flux on inner surface
    E1 = flux1/(4*pi*(r1^2));  % Electric Field at inner surface
    flux2 = (q2+q1)/e0;        % Flux on outer surface
    E2 = flux2/(4*pi*(r2^2));  % Electric Field at outer surface
end 
function [s] = sphr(r,q) 
    [x1,y1,z1] = sphere(50);
     x = x1*r;
     y = y1*r;
     z = z1*r;
     s = surf(x,y,z);
     hold on
     set(s,'FaceColor','#77AC30','EdgeColor','none');
     camlight
     % Making field vectors point inward if surface charge is negative
     if q<0
        [x2,y2,z2] = sphere(50);
        x = x2*(r+2);
        y = y2*(r+2);
        z = z2*(r+2);
        s2 = surf(x,y,z); % Creating an invisible large sphere so inverted vectors appear
        hold on
        set(s2,'FaceColor','none','EdgeColor','none');
        [Nx,Ny,Nz] = surfnorm(x,y,z);
        vec = quiver3(x,y,z,-Nx,-Ny,-Nz);
        set(vec,'AutoScale','on', 'AutoScaleFactor',2)
        axis equal
     else % For positive surface charges
     [Nx,Ny,Nz] = surfnorm(x,y,z);
     vec = quiver3(x,y,z,Nx,Ny,Nz);
     set(vec,'AutoScale','on', 'AutoScaleFactor',2)
     axis equal
     end
end
function [c] = crcl(r,theta) 
    O = [0,0,0];   % center
    deg = 0:0.01:2*pi;
    x = r*cos(deg);
    y = r*sin(deg);
    z = zeros(size(x));
    c = patch(x,y,z,'k');
    hold on
    view(3)
    set(c,'FaceColor','#77AC30','EdgeColor','none');
    d = [0 1 0];
    ang = 90-theta;
    rotate(c,d,ang);
    % Creating an invisible vertical plane for constant E field
    y2 = -r:0.25:r;
    z2 = -r:0.25:r;
    [Y2,Z2] = meshgrid(y2,z2);
    X2 = -r*ones(size(Y2));
    r2 = surf(X2,Y2,Z2);
    set(r2,'FaceColor','none','EdgeColor','none');
    [Nx,Ny,Nz] = surfnorm(X2,Y2,Z2);
    vec = quiver3(X2,Y2,Z2,Nx,Ny,Nz);
    set(vec,'AutoScale','on', 'AutoScaleFactor',r*5)
    axis equal    
end
function [r] = rct(w,l,theta) 
    x = 0:0.25:l;
    y = 0:0.25:w;
    z = 0:0;
    [X,Y,Z] = meshgrid(x,y,z);
    r = surf(X,Y,Z);
    hold on
    set(r,'FaceColor','#77AC30','EdgeColor','none');
    ang = 90-theta;
    d = [0 1 0];
    rotate(r,d,ang);
    % Creating an invisible vertical plane for constant E field
    y2 = 0:0.25:w;
    z2 = -l/2:0.25:l/2;
    [Y2,Z2] = meshgrid(y2,z2);
    X2 = ones(size(Y2));
    r2 = surf(X2,Y2,Z2);
    set(r2,'FaceColor','none','EdgeColor','none');
    [Nx,Ny,Nz] = surfnorm(X2,Y2,Z2);
    vec = quiver3(X2,Y2,Z2,Nx,Ny,Nz);
    set(vec,'AutoScale','on', 'AutoScaleFactor',l*2.5)
    axis equal
end
end