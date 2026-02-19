pl	0	#
		erase
		LOCATION 4500 31000 3500 31000
		limits 8 19 31 37
		expand 1.5
		lweight 3
		ticksize -1 0 -1 0
		box
		#
		rd grmonty.spec
		#
		histogram lw ll5 
		#
		xla "\nu [Hz]"
		yla "\nu L_\nu [erg s^{-1}]"
		#
		re speclab.m
		speclab
		#
		limits 8 19 31 37
		#
rd	1	#
		da $1
		#
		read {lw 1}
		read {l0 2  ta0 3  ts0 4  x10 5  x20 6  x30 7 }
		read {l1 8  ta1 9  ts1 10 x11 11 x21 12 x31 13}
		read {l2 14 ta2 15 ts2 16 x12 17 x22 18 x32 19}
		read {l3 20 ta3 21 ts3 22 x13 23 x23 24 x33 25}
		read {l4 26 ta4 27 ts4 28 x14 29 x24 30 x34 31}
		read {l5 32 ta5 33 ts5 34 x15 35 x25 36 x35 37}
		#
		set small = 1.e-12
		set ll0 = lg(l0+small) + lg(3.83e33)
		set ll1 = lg(l1+small) + lg(3.83e33)
		set ll2 = lg(l2+small) + lg(3.83e33)
		set ll3 = lg(l3+small) + lg(3.83e33)
		set ll4 = lg(l4+small) + lg(3.83e33)
		set ll5 = lg(l5+small) + lg(3.83e33)
		#
		set lw = lw + lg(9.1e-28*3.e10*3.e10/(6.626e-27))
		#
