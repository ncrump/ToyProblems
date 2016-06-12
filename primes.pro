pro primes
; print all primes less than 1000

p = [2,3]
; only check odds
for i = 5, 1000, 2 do begin
  prime = 1
  mx = n_elements(p)-1
  ; prime if not a multiple of other integers
  for j = 0, mx do begin
    if (i mod p[j] eq 0) then begin
      prime = 0
      break
    endif
  endfor
  if (prime eq 1) then p = [p,i]
endfor

; print primes
print, 'number of primes:'
print, n_elements(p)
print, 'list of primes:'
print, p

end
