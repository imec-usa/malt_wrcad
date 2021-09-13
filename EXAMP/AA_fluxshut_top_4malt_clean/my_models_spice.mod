* MODEL Declarations
*unshunted JJ with scalable parameters
.subckt rql_junction_scale PLUS MINUS PHI area=0.25 jjmod=rql25
b1 PLUS MINUS phi jjmod area='area*Xa*Xj*Xjcomp/timingX'
.ends rql_junction_scale
* rsj with scalable parameters
.subckt rql_rsj_scale PLUS MINUS jjmod=rql25 ic=0.25 icrn=0.7
Xb0 PLUS MINUS phi rql_junction_scale area='ic' jjmod='jjmod'
Xr0 PLUS MINUS rql_resistor_scale r='icrn/ic'
.ends rql_rsj_scale
* inductor models with scalable parameters
.subckt rql_inductor_scale PLUS MINUS l=1p
L0 PLUS MINUS 'l*Xl*Xlcomp/timingX'
.ends rql_inductor_scale* resistor models with scalable parameters
.subckt rql_resistor_scale PLUS MINUS r=1
R0 PLUS MINUS 'r*Xr*timingX'
.ends rql_resistor_scale
*JJ models
.model rql25 jj(rtype=1,cct=1,icon=10m,vg=2.6m,delv=0.1m,icrit=1m,r0=40,rn=1.800,cap=0.70p)
* End MODEL Declarations
