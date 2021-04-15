function [] = check_patternestimaterepeats(allpatternpitch,allpatternangle,allpatternphases)
% This function is for analyzing illumination pattern estimates from 
% repeated acquisitions as for e.g. live-cell datasets

% extract parameters
[numsteps,numchannels,numframes,numangles] = size(allpatternphases);

% only put in the effort if there are multiple repeated acquisitions
if numframes>1
  for jchannel = 1:numchannels
    allpatternphases = reshape(allpatternphases(:,jchannel,:,:),[numsteps,numframes,numangles]);
    allpatternpitch = reshape(allpatternpitch(jchannel,:,:),[numframes,numangles]);
    allpatternangle = reshape(allpatternangle(jchannel,:,:),[numframes,numangles]);

    % analyze pitch and angle estimates
    meanpitch = mean(allpatternpitch,1);
    stdpitch = std(allpatternpitch,[],1);
    meanangle = mean(allpatternangle,1);
    stdangle = std(allpatternangle,[],1);

    fprintf('relative precision pitch = %3.2e \n',mean(stdpitch)/mean(meanpitch))
    fprintf('precision pitch = %3.2e nm \n',mean(stdpitch))
    fprintf('precision angle = %3.2e deg \n',mean(stdangle)*180/pi)

    figure
    subplot(1,2,1)
    plot(allpatternpitch,'-o')
    ylabel('pitch [nm]')
    xlabel('frame')
    xlim([0 numframes+1])
    subplot(1,2,2)
    plot(allpatternangle*180/pi,'-o')
    ylabel('angle [deg]')
    xlabel('frame')
    xlim([0 numframes+1])

    figure
    subplot(1,2,1)
    plot(allpatternpitch-meanpitch,'-o')
    ylabel('\Delta pitch [nm]')
    xlabel('frame')
    xlim([0 numframes+1])
    subplot(1,2,2)
    plot(allpatternangle*180/pi-meanangle*180/pi,'-o')
    ylabel('\Delta angle [deg]')
    xlabel('frame')
    xlim([0 numframes+1])

    % analyze phase estimates
    meanphase = squeeze(mean(allpatternphases,2));
    stdphase = squeeze(std(allpatternphases,[],2));

    % extra computation to avoid phase wrapping induced biases
    stdphase_add = squeeze(std(mod(allpatternphases+pi,2*pi),[],2));
    stdphase = min(stdphase,stdphase_add);
    
    fprintf('precision phases = %3.2f deg \n',mean(stdphase,2)*180/pi)  

    figure
    for jangle = 1:numangles
      subplot(1,numangles,jangle)
      plot(squeeze(allpatternphases(:,:,jangle))'*180/pi,'-o')
      ylabel('phase [deg]')
      xlabel('frame')
      xlim([0 numframes+1])
      ylim([0 360])
      title(strcat('jangle=',num2str(jangle)))
    end

    figure
    scrsz = [1 1 1536 864];
    set(gcf,'Position',round([0.15*scrsz(3) 0.10*scrsz(4) 0.6*scrsz(3) 0.6*scrsz(4)]));
    for jstep = 1:numsteps
      for jangle = 1:numangles
        subplot(numangles,numsteps,numsteps*(jangle-1)+jstep)
        plot(squeeze(allpatternphases(jstep,:,jangle)-meanphase(jstep,jangle))'*180/pi,'-o')
        ylabel('\Delta phase [deg]')
        xlabel('frame')
        xlim([0 numframes+1])
      %   ylim([0 360])
        title(strcat('jstep=',num2str(jstep),', jangle=',num2str(jangle)))
      end
    end

  end
end

end