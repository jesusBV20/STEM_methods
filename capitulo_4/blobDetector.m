function [x,y] = blobDetector(I, sgm_list, delta, n_suav)
    % Inicializamos el tensor S(x,y,sigma)
    N = length(sgm_list); 
    [n,m] = size(I);
    S = zeros(n*n_suav,m*n_suav,N);
    
    % Inicializamos la figura con la imagen original
    figure(1), clf, subplot(1,N+2,1);
    imagesc(I)
    title("ABF inv. FFT STEM")
    xlabel("x [px]"), ylabel("y [px]")
    
    % Generamos grid del filtro
    x = linspace(-m, m, m);
    y = linspace(-n, n, n);
    [xx,yy] = meshgrid(x,y);
    
    % Generamos grid de suavizado
    xq = linspace(-m, m, m*n_suav);
    yq = linspace(-n, n, n*n_suav);
    [xq,yq] = meshgrid(xq,yq);
    
    % Bucle para los N sigmas
    for i = 1:N
        sigma = sgm_list(i);
        
        % Calculamos NLoG y lo convolucionamos
        NLoG = - 1/pi/sigma^2*(1 - (xx.^2 + yy.^2)/2/sigma^2).*exp(- (xx.^2 + yy.^2)/2/sigma^2);
        Iconv = conv2(I,NLoG,"same");
        
        subplot(1,N+2,i+1), imagesc(Iconv), axis off, colormap gray;
        title("$\sigma$ = " + num2str(sigma),'interpreter','latex')
        
        % Suavizamos el resultado por interpolaci?n
        if n_suav
            Iconv = griddata(xx,yy,abs(Iconv),xq,yq,"cubic");
        end

        % Almacenamos resultado y mostramos imagen filtrada
        S(:,:,i) = Iconv;
    end
    
    % Inicializamos la imagen resultado
    if n_suav
        I = imresize(I, [n*n_suav, m*n_suav]);
    end
    subplot(1,N+2,N+2)
    imagesc(I), axis off, colormap gray, hold on;

    % - Buscamos los puntos m?ximos para el primer sigma -
    i = 1;
    
    % Imponemos umbral delta
    Z = zeros(size(S(:,:,i)));
    for in = 1:size(S(:,:,i),1)
        for im = 1:size(S(:,:,i),2)
            if S(in,im,i)>delta
               Z(in,im) = S(in,im,i);
            else
               Z(in,im) = 0;
            end
        end
    end
    
    % Buscamos todos los m?ximos locales y los marcamos
    isMax = islocalmax(Z,1) & islocalmax(Z,2);
    isMax = imregionalmax(Z);
    [y,x] = find(isMax == 1);    
    
    plot(x,y,'r+', 'MarkerSize', 5, 'LineWidth', 0.1);
    title("Resultado $\sigma$ = " + num2str(sgm_list(i)),'interpreter','latex')
    legend("Puntos de inter?s")
    
    % Visualizaci?n de S y Z (interesante para depurar)
%     figure(3), clf;
%     surf(Iconv)
%     shading interp
%     figure(4), clf;
%     surf(S(:,:,i))
%     shading interp
%     figure(5), clf;
%     surf(Z)
%     shading interp
end

