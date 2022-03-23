function plotellipse(a,b,phi,x0,y0,format)
    th = 0:pi/100:2*pi;
    f_ellipse = x0+1j*y0 + exp(1i*phi)*(a*exp(1j*th) + b*exp(-1j*th));
    plot(real(f_ellipse),imag(f_ellipse),format);
end

