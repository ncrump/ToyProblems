pro taxicab
; find smallest number that can be written as
; the sum of two positive cubes two different ways
; taxicab = a^3 + b^3 = c^3 + d^3

limit = 20  ; search limit

sums = []
a    = []
b    = []
; brute force method
for i = 1, limit do begin
  for j = i, limit do begin
    s    = i*i*i + j*j*j
    indx = where(s eq sums)
    ; if sum exists then done
    if (indx gt 0) then begin
      print, 'taxicab number =',s
      print, 'pairs:'
      print, [a[indx],b[indx]],[i,j]
      break
    ; if sum does not exist then store and continue
    endif else begin
      sums = [sums, s]
      a    = [a, i]
      b    = [b, j]
    endelse
  endfor
endfor

end
