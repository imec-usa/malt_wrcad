
.control
load netlist.nom
load netlist.envelope
setplot tran1
set group
plot v(phi.Xb0) v(phi.Xb0.XI1) v(phi.Xb1.XI12)
unset group
set single
set color2 = "black"
set color3 = "blue"
set color4 = "black"
plot tran2.hi0 v(phi.Xb0) tran2.lo0
plot tran2.hi1 v(phi.Xb0.XI1) tran2.lo1
plot tran2.hi2 v(phi.Xb1.XI12) tran2.lo2
.endc
