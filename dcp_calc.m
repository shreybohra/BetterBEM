function [dcp, forces] = dcp_calc(twist, polars, z, R, tsr, B, chord, r0, V)

    % ==============
    % Algorithm 1 from Ning (2014)
    % ==============
    
    e = 1e-6;
    rho = 1.225;
    omega = (tsr*V)/R;
    
    if f(pi/2, twist, polars, z, R, tsr, B, chord, r0) > 0

        phi = fzero(@(phi) f(phi, twist, polars, z, R, ...
            tsr, B, chord, r0), [e, pi/2]);
        
    elseif fpb(-pi/4, twist, polars, z, R, tsr, B, chord, r0) < 0 ...
            && fpb(e, polars, z, R, tsr, B, chord, r0) > 0
        
        phi = fzero(@(phi) fpb(phi, twist, polars, z, R, ...
            tsr, B, chord, r0), [-pi/4, -e]);
 
    else        
        phi = fzero(@(phi) f(phi, twist, polars, z, R, ...
            tsr, B, chord, r0), [pi/2, pi]);
    end   
    
    % calculate  values from converged value of phi
    [a, ~, kdash, F] = afun(phi, twist, polars, z, R, B, chord, r0);
    adash = kdash/(1-kdash);
    
    aoa = phi - twist;
    
    if aoa<deg2rad(30) && aoa>deg2rad(-10)
        cla = polars{1};
        cda = polars{2};
        
    else
        cla = polars{3};
        cda = polars{4};
    end

    cl = cla(aoa);
    cd = cda(aoa);
    ld = cl/cd;
    
    % calculate dcp from Jamieson (2018)
    % Equation is split into multiple terms for readability
    
    t1 = 8*a*(1-a)*F*tsr*(z/R)^2;
    t2 = ld*(1-a) - (tsr*(z/R)*(1+adash));
    t3 = ld*tsr*(z/R)*(1+adash)+(1-a);
    
    dcp = (t1*t2)/t3;    
    
    cn = cl*cos(phi) + cd*sin(phi);
    ct = cl*sin(phi) - cd*cos(phi);
    vr = sqrt((V*(1-a))^2 + (omega*z*(1+adash))^2);
    
    Fn = 0.5*B*rho*vr^2*chord*cn;
    Ft = 0.5*B*rho*vr^2*chord*ct;
    
    forces = [Fn Ft];

    
end

function out = f(phi, twist, polars, r, R, tsr, B, chord, r0)

    [a, ~, kdash, ~] = afun(phi, twist, polars, r, R, B, chord, r0);
    local_tsr = tsr*(r/R);
    out = sin(phi)/(1 - a) - (cos(phi)*(1 - kdash))/(local_tsr);
   
end

function out = fpb(phi, twist, polars, r, R, tsr, B, chord, r0)

    [~, k, kdash, ~] = afun(phi, twist, polars, r, R, B, chord, r0);
    local_tsr = tsr*(r/R);    
    out = sin(phi)*(1 - k) - (cos(phi)*(1 - kdash))/(local_tsr);
    
end

function [a, k, kdash, F] = afun(phi, twist, polars, r, R, B, chord, r0)


    aoa = phi - twist;
    
    % Use the better fit in the range of data (-10 to 30 degrees)
    % Fall back on the less accurate 360 polar
    
    if aoa<deg2rad(30) && aoa>deg2rad(-10)
        cla = polars{1};
        cda = polars{2};
        
    else
        cla = polars{3};
        cda = polars{4};
    end


    sigma = (B*chord)/(2*pi*r); % (Lubbock, 2020)
            
    
    cl = cla(aoa);
    cd = cda(aoa);
    
    cn = cl*cos(phi) + cd*sin(phi);
    ct = cl*sin(phi) - cd*cos(phi);
    
    ftip = (B/2)*((R - r)/(r*abs(sin(phi))));
    Ftip = (2/pi) * acos(exp(-ftip));
    fhub = (B/2)*((r - r0)/(r0*abs(sin(phi))));
    Fhub = (2/pi) * acos(exp(-fhub));
    F = Ftip * Fhub;

    
    k = (sigma*cn)/(4*F*(sin(phi))^2);
    kdash = (sigma*ct)/(4*F*sin(phi)*cos(phi));
    
    if k < 2/3
        a = k/(1+k);
        
    else
        t1 = 2*F*k - (10/9 - F);
        t2 = 2*F*k - F*(4/3 - F);
        t3 = 2*F*k - (25/9 - 2*F);
        
        a = (t1 - sqrt(t2))/t3;
    end
end