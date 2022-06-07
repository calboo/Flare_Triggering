pro fig9

; IDL script used to produce Fig9 in 
; The Effects of Oscillations & Collisions of Emerging Bipolar Regions
; on the Triggering of Solar Flares
; C.Boocock, K.Kusano, D. Tsiklauri

; For each timestep the magnetic field along the polarity inversion
; line B_x is loaded and plotted through the plane x=0 as a coloured
; contour. Additionally the vertical magnetic field B_z is loaded and
; plotted across the photospheric surface z=0 as a coloured contour.
  
; To produce Fig9 this script was used on the output from a simulation
; with two non-colliding injected fields

; Cycle through the relevant outputs
for i = 0,16,4 do begin
; Load the magnetic field data
ds = getdata(i,/bx,/bz)

; Plot the coloured contour of Bx
c0 = contour(reform(ds.bx[*,50,*]), /fill, n_levels=50, $
	XTICKINTERVAL = 100, RGB_TABLE=colortable(3),$
	XTITLE =  'Y / Mm', YTITLE = 'Z / Mm', $
        XTICKNAME = ['-40','-20','0','20','40'], $
        YTICKNAME = ['0','10','20','30','40'], $
        POSITION = [0.05,0.2,0.95,0.95], $
        DIMENSIONS  = [1100,600], $
        FONT_SIZE=16)
c0.Scale, 0.9, 0.9

; Add colour bar
cb0 = COLORBAR( $
  TARGET = c0, $
  TICKFORMAT = '(F5.2)', $
  FONT_SIZE = 16, MAJOR=6, $
  POSITION = [0.2, 0.10, 0.8, 0.13], $
  TITLE='Field Strength B!D x')

; Save image to a file
filename = 'cont1'+String(i, Format='(I3.3)')+'.png'
c0.save, filename

; Plot the coloured contour of Bz
c1 = contour(reform(ds.bz[*,*,0]), /fill, n_levels=50, $
	XTITLE =  'X / Mm', YTITLE = 'Y / Mm', $
        RGB_TABLE=colortable(3), $
	XTICKNAME = ['-40','-20','0','20','40'], $
        YTICKNAME = ['-10','5','0','5','10'], $
        DIMENSIONS=[1200,500],$
        FONT_SIZE=16)
c1.Scale, 1.0, 0.65

; Add colour bar
cb1 = COLORBAR( $
  TARGET = c1, $
  TICKFORMAT = '(F5.2)', $
  FONT_SIZE = 16, MAJOR=6, $
  POSITION = [0.2, 0.11, 0.8, 0.14], $
  TITLE='Field Strength B!D z')

; Save image to file
filename = 'cont2'+String(i, Format='(I3.3)')+'.png'
c1.save, filename

endfor

end

