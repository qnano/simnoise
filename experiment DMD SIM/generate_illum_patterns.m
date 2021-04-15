function I = generate_illum_patterns(N, pattern_parameter1, pattern_parameter2, method)
% generation of the binary DMD illumination patterns for Fourier Ptychography simulation
%
% N is the size of the illumination pattern NxN
%
% pattern parameter1 is: pitch of the sine patterns in case the method is 'SIM'
%                       pitch of the multispot pattern in case the method is 'MSP'
%                       sparsity S of the pseudo-random pattern in case the method is 'PRP' 
%
% pattern parameter2 is the total number of patterns
%
% method is defined by a string: 'SIM' - for line patterns, ***the number of shifts and rotations can be changed*** 
%                                'MSP' - for multispot patterns
%                                'PRP' - for pseudo-random patterns   
% Author: Nadya Chakrova 
% September 2015


if method == 'SIM'
        
        pitch = pattern_parameter1;
        N_periods = ceil(N/pitch);
        N_ang = 9;
        N_shifts = 9;
        K = N_ang*N_shifts;
        I = zeros(N,N,K);
        angle = pi/180*linspace(0,180,N_ang+1);
        j=0;
        for i = 1:N_ang
            for sh = 1:N_shifts % 3 phase shifts of the line pattern, as in classical SIM 
            j = j+1;
            shift = (sh-1)/N_shifts;
            c1 = cos(angle(i));
            c2 = sin(angle(i));
            A = cos(N_periods.*(xx(N,N, 'radfreq')*c1 + yy(N,N, 'radfreq')*c2) + 2*pi*shift);
            A(A>0)=1;
            A(A<0)=0;
            w=tukeywin(N,0.1)*tukeywin(N,0.1)';
            A = im2mat(A).*w;
            I(:,:,j) = A./(K/2);
            end
        end

elseif method == 'PRP'
        
        S = pattern_parameter1;
        K = pattern_parameter2;
        indexes = zeros(N,N,S);
        for js = 1:S
          indexes(:,:,js) = ceil((K-js+1)*rand(N,N));
          if (js>1)
            shiftnums = zeros(N,N);
            for jss = 1:js-1
              shiftnums = shiftnums + double(indexes(:,:,js)>=indexes(:,:,jss));
            end
            indexes(:,:,js) = indexes(:,:,js) + shiftnums;
            flag = 1;
            while flag
              shiftnums = zeros(N,N);
              for jss = 1:js-1
                shiftnums = shiftnums + double(indexes(:,:,js)==indexes(:,:,jss));
              end
              indexes(:,:,js) = indexes(:,:,js) + shiftnums;
              if (sum(sum(shiftnums))==0)
                flag = 0;
              end
            end
          end
        end

        I = zeros(N,N,K);
        for kk = 1:N
          for ll = 1:N
            I(kk,ll,indexes(kk,ll,1:S)) = 1/S;
          end
        end

elseif method == 'MSP'

        block = 1;
        step = pattern_parameter1;
        K = (step/block)^2;
        nSizeX = N;
        nSizeY = N;

        pImageData = zeros(N, (step/block)^2*N, 'uint8');
        w=tukeywin(N,0.1)*tukeywin(N,0.1)';

        for l = 1:(step/block)

            for i = 1:block
                for j = 1:block
             pImageData(i:step:nSizeX, j+(step/block*(l-1))*nSizeY+block*(l-1):step:(step/block*(l-1)+1)*nSizeY+block*(l-1)) = 1;
                end
            end

            for k = 1:(step/block)-1
            pImageData(k*block+1:nSizeX, (k+(step/block)*(l-1))*nSizeY+1:(k+1+(step/block)*(l-1))*nSizeY) = pImageData(1:nSizeX-block*k, nSizeY*(step/block)*(l-1)+1:nSizeY+nSizeY*(step/block)*(l-1));
            end

        end

        ind1 = 0;
        ind2 = 0;
        for i = 1:(step/block)^2
            ind1 = ind2+1;
            ind2 = N*i;
            I(:,:,i) = pImageData(:, ind1:ind2);
        end

end

end

