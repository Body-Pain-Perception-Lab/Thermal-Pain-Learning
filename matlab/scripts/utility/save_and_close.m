%model comparison
function c = save_and_close(filename, width, height)

  hFig = gcf;  % Get the current figure handle
  set(hFig, 'Position', [100, 100, width, height]);
  
  % Save the subplot as an image file (e.g., PNG or PDF)
  saveas(hFig, filename);
  
  % Close the figure (optional, depending on your preference)
  close(hFig);
end