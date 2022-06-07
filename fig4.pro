pro fig4 

; IDL script used to produce Fig4 in 
; The Effects of Oscillations & Collisions of Emerging Bipolar Regions
; on the Triggering of Solar Flares
; C.Boocock, K.Kusano, D. Tsiklauri

; The magnetic field along the polarity inversion line B_x is loaded at
; each timestep and B_x through the plane x=0 is  plotted as a
; coloured contour for each output. 

; To produce Fig4 this script was used on the output from a simulation
; with  Ax = 0.2 (4 Mm), Ay = 0.2 (4 Mm) and frequency = 1

; Set window size
w = window(dimensions=[1500,800])

; Set contour levels
n_levels = 50
levels = (FINDGEN(n_levels)/(n_levels-1)*2)-1

; Plot contour at t = 0 and label vertical axis
ds = getdata(0,/bx)
c0 = contour(reform(ds.bx[200,*,*]), /fill, n_levels=50, $
	XTICKINTERVAL = 25, RGB_TABLE=colortable(3),$
        C_VALUE = levels,$
	XTITLE =  'Y / Mm', YTITLE = 'Z / Mm', $
	XTICKNAME = ['-10','-5','0','5','10'], $
        YTICKNAME = ['0','10','20','30','40'], $
        POSITION= [0.05,0.2,0.20,0.95], /CURRENT)

; Plot contour at first output time
ds = getdata(1,/bx)
c1 = contour(reform(ds.bx[200,*,*]), /fill, n_levels=50, $
	XTICKINTERVAL = 25, RGB_TABLE=colortable(3),$
        C_VALUE = levels,$
	XTITLE =  'Y / Mm', $
	XTICKNAME = ['-10','-5','0','5','10'], $
        POSITION= [0.24,0.2,0.39,0.95], $
        YSHOWTEXT = 0, /CURRENT)

; Plot contour at second output time
ds = getdata(2,/bx)
c2 = contour(reform(ds.bx[200,*,*]), /fill, n_levels=50, $
	XTICKINTERVAL = 25, RGB_TABLE=colortable(3),$
        C_VALUE = levels,$
	XTITLE =  'Y / Mm', $
	XTICKNAME = ['-10','-5','0','5','10'], $
        POSITION= [0.43,0.2,0.58,0.95], $
        YSHOWTEXT = 0, /CURRENT)

; Plot contour at third output time
ds = getdata(3,/bx)
c3 = contour(reform(ds.bx[200,*,*]), /fill, n_levels=50, $
	XTICKINTERVAL = 25, RGB_TABLE=colortable(3),$
        C_VALUE = levels,$
	XTITLE =  'Y / Mm', $
	XTICKNAME = ['-10','-5','0','5','10'], $
        POSITION= [0.62,0.2,0.77,0.95], $
        YSHOWTEXT = 0, /CURRENT)

; Plot contour at fourth output time
ds = getdata(4,/bx)
c4 = contour(reform(ds.bx[200,*,*]), /fill, n_levels=50, $
	XTICKINTERVAL = 25, RGB_TABLE=colortable(3),$
        C_VALUE = levels,$
	XTITLE =  'Y / Mm', $
	XTICKNAME = ['-10','-5','0','5','10'], $
        POSITION= [0.81,0.2,0.96,0.95], $
        YSHOWTEXT=0, /CURRENT)

; Add colour bar
tick_labels = [STRTRIM(FIX(levels), 2), '']
cb = COLORBAR( $
  TARGET = c1, $
  TICKFORMAT = '(F5.2)', $
  FONT_SIZE = 12, $
  POSITION = [0.2, 0.09, 0.8, 0.12], $
  TITLE='Field Strength')

end



