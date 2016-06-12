pro sampleread
; reads unformatted text file

; set file location
file = '../data/sampletext.txt'
; set array to store text
text = []
line = ''

; open file
openr, lun, file, /get_lun

; loop over lines
while not eof(lun) do begin
  readf, lun, line
  text = [text,line]
  endwhile

; close file
free_lun, lun

; print lines
nlines = size(text, /n_elements)
print, 'lines:', nlines
for i = 0, nlines-1 do print, text[i]
end
