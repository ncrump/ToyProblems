pro pi
  tic
  ; define variables
  n = 1000000000L
  s = 0.0d
  ; compute pi
  for i = 0, n-1 do begin
    s += 1.0d/((4.0d*i+1.0d)*(4.0d*i+3.0d))
  endfor
  pi = double(8.0d*s)
  ; print results
  print, "approximation to pi =",pi
  toc
end
