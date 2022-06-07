pro fig6

; IDL script used to produce Fig6 in 
; The Effects of Oscillations & Collisions of Emerging Bipolar Regions
; on the Triggering of Solar Flares
; C.Boocock, K.Kusano, D. Tsiklauri

; Energy data is loaded from simulations in seperate directories
; The maximum kinetic energy from each run is stored in an array
; The maximum peak KE is from each run is then used to plot a coloured
; contour showing the variation of peak KE with A_x and A_y

; For these simulations frequency = 1 was fixed
; Ax was varied between 0 and 20 Mm and
; Ay was varied between 0 and 4 Mm

; Set up plot and array for peak KE values
SET_PLOT, 'X'
max_ke = dblarr(6,5)

; Read peak KE for simulations with A_y = 0 (0 Mm)

print, 'ay = 0.0'

en = getenergy(wkdir="../0.0_0.0/Data")
print, max(en.en_ke)
max_ke[0,0] = max(en.en_ke)
en = getenergy(wkdir="../0.0_0.2/Data")
print, max(en.en_ke)
max_ke[1,0] = max(en.en_ke)
en = getenergy(wkdir="../0.0_0.4/Data")
print, max(en.en_ke)
max_ke[2,0] = max(en.en_ke)
en = getenergy(wkdir="../0.0_0.6/Data")
print, max(en.en_ke)
max_ke[3,0] = max(en.en_ke)
en = getenergy(wkdir="../0.0_0.8/Data")
print, max(en.en_ke)
max_ke[4,0] = max(en.en_ke)
en = getenergy(wkdir="../0.0_1.0/Data")
print, max(en.en_ke)
max_ke[5,0] = max(en.en_ke)

; Read peak KE for simulations with A_y = 0.05 (1 Mm)

print, 'ay = 0.05'

en = getenergy(wkdir="../0.05_0.0/Data")
print, max(en.en_ke)
max_ke[0,1] = max(en.en_ke)
en = getenergy(wkdir="../0.05_0.2/Data")
print, max(en.en_ke)
max_ke[1,1] = max(en.en_ke)
en = getenergy(wkdir="../0.05_0.4/Data")
print, max(en.en_ke)
max_ke[2,1] = max(en.en_ke)
en = getenergy(wkdir="../0.05_0.6/Data")
print, max(en.en_ke)
max_ke[3,1] = max(en.en_ke)
en = getenergy(wkdir="../0.05_0.8/Data")
print, max(en.en_ke)
max_ke[4,1] = max(en.en_ke)
en = getenergy(wkdir="../0.05_1.0/Data")
print, max(en.en_ke)
max_ke[5,1] = max(en.en_ke)

; Read peak KE for simulations with A_y = 0.1 (2 Mm)

print, 'ay = 0.1'

en = getenergy(wkdir="../0.1_0.0/Data")
print, max(en.en_ke)
max_ke[0,2] = max(en.en_ke)
en = getenergy(wkdir="../0.1_0.2/Data")
print, max(en.en_ke)
max_ke[1,2] = max(en.en_ke)
en = getenergy(wkdir="../0.1_0.4/Data")
print, max(en.en_ke)
max_ke[2,2] = max(en.en_ke)
en = getenergy(wkdir="../0.1_0.6/Data")
print, max(en.en_ke)
max_ke[3,2] = max(en.en_ke)
en = getenergy(wkdir="../0.1_0.8/Data")
print, max(en.en_ke)
max_ke[4,2] = max(en.en_ke)
en = getenergy(wkdir="../0.1_1.0/Data")
print, max(en.en_ke)
max_ke[5,2] = max(en.en_ke)

; Read peak KE for simulations with A_y = 0.15 (3 Mm)

print, 'ay = 0.15'

en = getenergy(wkdir="../0.15_0.0/Data")
print, max(en.en_ke)
max_ke[0,3] = max(en.en_ke)
en = getenergy(wkdir="../0.15_0.2/Data")
print, max(en.en_ke)
max_ke[1,3] = max(en.en_ke)
en = getenergy(wkdir="../0.15_0.4/Data")
print, max(en.en_ke)
max_ke[2,3] = max(en.en_ke)
en = getenergy(wkdir="../0.15_0.6/Data")
print, max(en.en_ke)
max_ke[3,3] = max(en.en_ke)
en = getenergy(wkdir="../0.15_0.8/Data")
print, max(en.en_ke)
max_ke[4,3] = max(en.en_ke)
en = getenergy(wkdir="../0.15_1.0/Data")
print, max(en.en_ke)
max_ke[5,3] = max(en.en_ke)

; Read peak KE for simulations with A_y = 0.2 (4 Mm)

print, 'ay = 0.2'

en = getenergy(wkdir="../../production_1/0_w1/Data")
print, max(en.en_ke)
max_ke[0,4] = max(en.en_ke)
en = getenergy(wkdir="../../production_1/A_w1/Data")
print, max(en.en_ke)
max_ke[1,4] = max(en.en_ke)
en = getenergy(wkdir="../../production_1/B_w1/Data")
print, max(en.en_ke)
max_ke[2,4] = max(en.en_ke)
en = getenergy(wkdir="../../production_1/C_w1/Data")
print, max(en.en_ke)
max_ke[3,4] = max(en.en_ke)
en = getenergy(wkdir="../../production_1/D_w1/Data")
print, max(en.en_ke)
max_ke[4,4] = max(en.en_ke)
en = getenergy(wkdir="../../production_1/E_w1/Data")
print, max(en.en_ke)
max_ke[5,4] = max(en.en_ke)

; Set up markers, axis and contour levels

X = [0.0,0.2,0.4,0.6,0.8,1.0,0.0,0.2,0.4,0.6,0.8,1.0,0.0,0.2,0.4,0.6,0.8,1.0,0.0,0.2,0.4,0.6,0.8,1.0,0.0,0.2,0.4,0.6,0.8,1.0]*2e1
Y = [0,0,0,0,0,0,0.05,0.05,0.05,0.05,0.05,0.05,0.1,0.1,0.1,0.1,0.1,0.1,0.15,0.15,0.15,0.15,0.15,0.15,0.2,0.2,0.2,0.2,0.2,0.2]*2e1

cx = [0.0,0.2,0.4,0.6,0.8,1.0]*2e1
cy = [0.0,0.05,01,0.15,0.2]*2e1

clevels = [0.0125,0.015,0.0175,0.02,0.0225,0.025,0.0275,0.03,0.0325,0.035,0.0375,0.04]

; Plot the coloured contour

c = CONTOUR(max_ke,cx,cy, /FILL, RGB_TABLE=ct,c_value=clevels,$
XTITLE='A !D x !N / Mm', YTITLE='A !D y !N / Mm', $
POSITION=[0.1,0.22,0.95,0.95], YRANGE=[0,4])

; Add markers to plot

dist_size = SCATTERPLOT(X,Y, /SYM_FILLED, RGB_TABLE=7, /OVERPLOT)

; Add colorbar

cb = COLORBAR(TARGET=c,POSITION=[0.1,0.09,0.9,0.12], $
TITLE='Peak Kinetic Energy')

; Print table of peak KE

print, 'Table'
print, max_ke

end
