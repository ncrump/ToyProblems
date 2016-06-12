pro sp500data
; analyze sp500 data file

; set file location
file = '../data/table.csv'

; read text
temp = rd_tfile(file, /hskip)

; get size
n = n_elements(temp)
m = n_elements(str2arr(temp,','))

; array of structures to store data
str  = {date:'',open:0.0D,high:0.0D,low:0.0D,close:0.0D,vol:0.0D,adjclos:0.0D}
data = replicate(str,n)

; parse text into structures
for i = 0, n-1 do begin
  s  = str2arr(temp[i], ',')
  for j = 0, m-1 do begin
    data[i].(j) = s[j]
  endfor
endfor

; get statistics
diff = data.close - data.open
momt = moment(diff)
mean = momt[0]
stdv = sqrt(momt[1])

; get histogram and gauss fit
y    = histogram(diff,locations=x)
yfit = mpfitpeak(x,y,coeff,nterms=3)

; set distribution width for plotting
width = 6
xmin  = mean - width*stdv
xmax  = mean + width*stdv

; print results
print, 'set mean:', mean
print, 'set stdv:', stdv
print, 'fit mean:', coeff[1]
print, 'fit stdv:', coeff[2]

; plot results
device, decomposed=0
!p.background=255
!p.color=0
plot,x,y,xrange=[xmin,xmax],charsize=2,$
     xtitle='S&P 500 Closing Difference (1950 to 2015)',$
     ytitle='Frequency'
oplot,x,yfit,linestyle=3
plt = tvrd(true=1)
write_png, 'Plots/SP500.png', plt

end
