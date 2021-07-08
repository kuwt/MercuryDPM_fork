%%% This function takes fig_handle, outputfile_name, width of figure (cm)
function PrintLaTeX(figno, filename, width,factor)

saveas(figno,[filename '.pdf'])

pause(.1)

filename=strcat(filename,'.matlab');

set(0, 'DefaultTextInterpreter', 'none');

if (~exist('factor','var'))
  factor=0.8;
end

laprint(figno,filename, 'keepfontprops', 'on', 'asonscreen', 'on', 'head', 'off', 'width', width,'factor',factor,'keepticklabels','off','figcopy','off');
