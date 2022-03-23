function plotcircle(r,x0,y0,format)
    th = 0:pi/100:2*pi;
    f_circle = r*exp(1j*th) + x0+1j*y0;
    plot(real(f_circle),imag(f_circle),format);
end

