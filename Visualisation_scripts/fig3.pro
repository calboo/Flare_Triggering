pro fig3

; IDL script used to produce Fig3 in 
; The Effects of Oscillations & Collisions of Emerging Bipolar Regions
; on the Triggering of Solar Flares
; C.Boocock, K.Kusano, D. Tsiklauri

; Energy data is loaded from simulations in seperate directories
; The kinetic energies over each simulation are then plotted

; For these simulations Ay = 0.2 (4 Mm) and frequency = 1 were fixed
; Ax was varied between 0 and 20 Mm

; Load data from run with Ax = 1 (20 Mm), Ay = 0.2 (4 Mm) and frequency = 1
en = getenergy(WKDIR='../E_w1/Data')

; Set up graphing window
loadct,13
TVLCT, 255, 255, 255, 254 ; White color 
TVLCT, 0, 0, 0, 253       ; Black color 
!P.Color = 253  
!P.Background = 254 
plot, en.time,en.en_ke,THICK=2, XTITLE='Time / !7s!1 !I A!N', YTITLE='Kinetic Energy', XRANGE=[1,50]
oplot, en.time,en.en_ke,THICK=3,color=250

; Load data from run with Ax = 0.8 (16 Mm), Ay = 0.2 (4 Mm) and frequency = 1
en1 = getenergy(WKDIR='../D_w1/Data')
oplot, en1.time,en1.en_ke,THICK=2,line=0,color=220

; Load data from run with Ax = 0.6 (12 Mm), Ay = 0.2 (4 Mm) and frequency = 1
en2 = getenergy(WKDIR='../C_w1/Data')
oplot, en2.time,en2.en_ke,THICK=2,line=4,color=150

; Load data from run with Ax = 0.4 (8 Mm), Ay = 0.2 (4 Mm) and frequency = 1
en3 = getenergy(WKDIR='../B_w1/Data')
oplot, en3.time,en3.en_ke,THICK=2,line=3,color=100

; Load data from run with Ax = 0.2 (4 Mm), Ay = 0.2 (4 Mm) and frequency = 1
en4 = getenergy(WKDIR='../A_w1/Data')
oplot, en4.time,en4.en_ke,THICK=2,line=2,color=80

; Load data from run with Ax = 0 (0 Mm), Ay = 0.2 (4 Mm) and frequency = 1
en5 = getenergy(WKDIR='../0_w1/Data')
oplot, en5.time,en5.en_ke,THICK=2,line=6,color=0

end
