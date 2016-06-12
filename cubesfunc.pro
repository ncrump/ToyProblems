function cubesfunc, n, limit
; returns solutions to a^3 + b^3

; check for invalid input
if (abs(n) ne n or ceil(n) ne n) then begin
  return, 'input must be a positive integer'
endif else begin
  ; brute force method
  a = []
  b = []
  ; loop over a's and b's
  for i = 1, limit do begin
    for j = i, limit do begin
        ; check sum of cubes
        if (long(i)*i*i + long(j)*j*j eq n) then begin
          a = [a,i]
          b = [b,j]
        endif
    endfor
  endfor
  return, [a,b]
endelse

end
