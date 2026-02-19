speclab 0       #
		set me = 9.11e-28
		set c = 2.998e10
		set h = 6.63e-27
                set eV = 1.602e-12
		set ang = 1.e-8
		set mum = 1.e-4
		#
                set kev = 1.e3*eV
                set mev = 1.e6*eV
		#
                define EX   (lg(1.*kev/h))
		define Eg   (lg(1.*mev/h))
		define Eo   (lg(c/(5000*ang)))
		define Enir (lg(c/(2.*mum)))
		define Emir (lg(c/(10.*mum)))
		define Efir (lg(c/(100.*mum)))
		define Emm  (lg(c/(0.1)))
		define EHI  (lg(1420.4058e6))
		#
		limits $fx1 $fx2 0 1
		linemark $Eg "1 MeV"
		linemark $EX "1 keV"
		linemark $Eo "5000 \AA"
		linemark $Enir "2 \mu m"
		linemark $Emir "10 \mu m"
		linemark $Efir "100 \mu m"
		linemark $Emm "1 mm"
		linemark $EHI "21 cm"
		#
linemark 2	#
		if($1 > $fx1 && $1 < $fx2) {\
			angle 270
			expand 1.2
			#
			relocate $1 0.9 draw $1 0.95
			relocate $1 0.89
			putlabel 6 $2
			#
			expand 1.5
			angle 0
		}
		

