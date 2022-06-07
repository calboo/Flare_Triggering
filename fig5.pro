pro fig5

; IDL script used to produce Fig5 in 
; The Effects of Oscillations & Collisions of Emerging Bipolar Regions
; on the Triggering of Solar Flares
; C.Boocock, K.Kusano, D. Tsiklauri

; Energy data is loaded from simulations in seperate directories
; The maximum kinetic energy from each run is stored in an array
; The maximum peak KE is from each run is then used to plot a coloured
; contour showing the variation of peak KE with A_x and frequency

; For these simulations Ay = 0.2 (4 Mm) was fixed
; Ax was varied between 0 and 20 Mm and
; the frequency was varied from 0.1 to 5

; Set up plot and array for peak KE values
SET_PLOT, 'X'
max_ke = dblarr(5,6)

; Read peak KE for simulations with A_x = 0 (0 Mm)

print, '0'

en = getenergy(wkdir="../0_w0.1/Data")
print, max(en.en_ke)
max_ke[0,0] = max(en.en_ke)
en = getenergy(wkdir="../0_w0.5/Data")
print, max(en.en_ke)
max_ke[1,0] = max(en.en_ke)
en = getenergy(wkdir="../0_w1/Data")
print, max(en.en_ke)
max_ke[2,0] = max(en.en_ke)
en = getenergy(wkdir="../0_w2/Data")
print, max(en.en_ke)
max_ke[3,0] = max(en.en_ke)
en = getenergy(wkdir="../0_w5/Data")
print, max(en.en_ke)
max_ke[4,0] = max(en.en_ke)

; Read peak KE for simulations with A_x = 0.2 (4 Mm)

print, 'A'

en = getenergy(wkdir="../A_w0.1/Data")
print, max(en.en_ke)
max_ke[0,1] = max(en.en_ke)
en = getenergy(wkdir="../A_w0.5/Data")
print, max(en.en_ke)
max_ke[1,1] = max(en.en_ke)
en = getenergy(wkdir="../A_w1/Data")
print, max(en.en_ke)
max_ke[2,1] = max(en.en_ke)
en = getenergy(wkdir="../A_w2/Data")
print, max(en.en_ke)
max_ke[3,1] = max(en.en_ke)
en = getenergy(wkdir="../A_w5/Data")
print, max(en.en_ke)
max_ke[4,1] = max(en.en_ke)

; Read peak KE for simulations with A_x = 0.4 (8 Mm)

print, 'B'

en = getenergy(wkdir="../B_w0.1/Data")
print, max(en.en_ke)
max_ke[0,2] = max(en.en_ke)
en = getenergy(wkdir="../B_w0.5/Data")
print, max(en.en_ke)
max_ke[1,2] = max(en.en_ke)
en = getenergy(wkdir="../B_w1/Data")
print, max(en.en_ke)
max_ke[2,2] = max(en.en_ke)
en = getenergy(wkdir="../B_w2/Data")
print, max(en.en_ke)
max_ke[3,2] = max(en.en_ke)
en = getenergy(wkdir="../B_w5/Data")
print, max(en.en_ke)
max_ke[4,2] = max(en.en_ke)

; Read peak KE for simulations with A_x = 0.6 (12 Mm)

print, 'C'

en = getenergy(wkdir="../C_w0.1/Data")
print, max(en.en_ke)
max_ke[0,3] = max(en.en_ke)
en = getenergy(wkdir="../C_w0.5/Data")
print, max(en.en_ke)
max_ke[1,3] = max(en.en_ke)
en = getenergy(wkdir="../C_w1/Data")
print, max(en.en_ke)
max_ke[2,3] = max(en.en_ke)
en = getenergy(wkdir="../C_w2/Data")
print, max(en.en_ke)
max_ke[3,3] = max(en.en_ke)
en = getenergy(wkdir="../C_w5/Data")
print, max(en.en_ke)
max_ke[4,3] = max(en.en_ke)

; Read peak KE for simulations with A_x = 0.6 (12 Mm)

print, 'D'

en = getenergy(wkdir="../D_w0.1/Data")
print, max(en.en_ke)
max_ke[0,4] = max(en.en_ke)
en = getenergy(wkdir="../D_w0.5/Data")
print, max(en.en_ke)
max_ke[1,4] = max(en.en_ke)
en = getenergy(wkdir="../D_w1/Data")
print, max(en.en_ke)
max_ke[2,4] = max(en.en_ke)
en = getenergy(wkdir="../D_w2/Data")
print, max(en.en_ke)
max_ke[3,4] = max(en.en_ke)
en = getenergy(wkdir="../D_w5/Data")
print, max(en.en_ke)
max_ke[4,4] = max(en.en_ke)

; Read peak KE for simulations with A_x = 1 (20 Mm)

print, 'E'

en = getenergy(wkdir="../E_w0.1/Data")
print, max(en.en_ke)
max_ke[0,5] = max(en.en_ke)
en = getenergy(wkdir="../E_w0.5/Data")
print, max(en.en_ke)
max_ke[1,5] = max(en.en_ke)
en = getenergy(wkdir="../E_w1/Data")
print, max(en.en_ke)
max_ke[2,5] = max(en.en_ke)
en = getenergy(wkdir="../E_w2/Data")
print, max(en.en_ke)
max_ke[3,5] = max(en.en_ke)
en = getenergy(wkdir="../E_w5/Data")
print, max(en.en_ke)
max_ke[4,5] = max(en.en_ke)

; Set up markers, axis and contour levels

X = [0.1,0.5,1.0,2.0,5.0,0.1,0.5,1.0,2.0,5.0,0.1,0.5,1.0,2.0,5.0,0.1,0.5,1.0,2.0,5.0,0.1,0.5,1.0,2.0,5.0,0.1,0.5,1.0,2.0,5.0]
Y = [0,0,0,0,0,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.6,0.6,0.6,0.6,0.6,0.8,0.8,0.8,0.8,0.8,1.0,1.0,1.0,1.0,1.0]*2e1

cx = [0.1,0.5,1.0,2.0,5.0]
cy = [0.0,0.2,0.4,0.6,0.8,1.0]*2e1

clevels = [0.0125,0.015,0.0175,0.02,0.0225,0.025,0.0275,0.03,0.0325,0.035,0.0375,0.04]

; Plot the coloured contour

c = CONTOUR(max_ke,cx,cy, /FILL, RGB_TABLE=ct,c_value=clevels,$
XTITLE='Frequency / $\tau$ !D A !N !U -1', YTITLE='Amplitude / Mm', $
POSITION=[0.1,0.22,0.95,0.95])

; Add markers to plot

dist_size = SCATTERPLOT(X,Y, /SYM_FILLED, RGB_TABLE=7, /OVERPLOT)

; Add colorbar

cb = COLORBAR(TARGET=c,POSITION=[0.1,0.09,0.9,0.12], $
TITLE='Peak Kinetic Energy')

; Print table of peak KE

print, 'Table'
print, max_ke

end
