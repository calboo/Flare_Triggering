pro fig8

; IDL script used to produce Fig8 in 
; The Effects of Oscillations & Collisions of Emerging Bipolar Regions
; on the Triggering of Solar Flares
; C.Boocock, K.Kusano, D. Tsiklauri

; Energy data is loaded from simulations in seperate directories
; The kinetic energies over each simulation are then plotted

; Graphs of KE against time are plotted for simulations with a
; non-oscillating injected field, two stationary injected fields and
; two colliding injected fields
; For both simulations the field shear is theta_0 = 80 degrees and the
; emerging region has opposite polarity to the background phi_e = 180 degrees

; Load data from the run with two colliding injected fields
en = getenergy(WKDIR='../../Collision/Collide_static/Data')   

loadct,13
TVLCT, 255, 255, 255, 254 ; White color 
TVLCT, 0, 0, 0, 253       ; Black color 
!P.Color = 253  
!P.Background = 254 
plot, en.time,en.en_ke,THICK=3, XTITLE='Time / !7s!1 !I A!N', YTITLE='Kinetic Energy', XRANGE=[1,50]
oplot, en.time,en.en_ke,THICK=3,color=250

; Load data from the run with two non-colliding injected fields
en1 = getenergy(WKDIR='../../Collision/NoCollide_static/Data')
oplot, en1.time,en1.en_ke,THICK=2,line=3,color=80

; Load data from the run with a single stationary injected field
en2 = getenergy(WKDIR='../../Standard/80_180/Data')
oplot, en2.time,en2.en_ke,THICK=2,line=6,color=0

end
