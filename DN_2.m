(* ::Package:: *)

(* ::Title:: *)
(*DN2*)


(* ::Section:: *)
(*1.1 Ra\[CHacek]unanje vrednosti \[SHacek]tevila \[Pi] - Monte Carlo*)


(* ::Input:: *)
(*(*Tu pokli\[CHacek]em datote\[CHacek]no funkcijo, ki sem jo definiral v datoteki s kon\[CHacek]nico .m, katera vrne koordinate to\[CHacek]k, ki le\[ZHacek]ijo znotraj kvadrata in koordinate to\[CHacek]k, ki le\[ZHacek]ijo znotraj kroga*)*)
(*Get["I:\\FAKS\\3.LETNIK\\5. semester\\Napredna ra\[CHacek]unalni\[SHacek]ka orodja\\DN\\DN2\\mccpi.m"]*)


(* ::Input:: *)
(*(*Definiram glavno funkcijo, katera izra\[CHacek]una pribli\[ZHacek]no vrednost \[Pi], ter njeno napako*)*)
(*calcpi[st_]:=Module[{pipribl,napaka},*)
(*(*Prika\[ZHacek]em okno, v katerega vpi\[SHacek]em \[SHacek]tevilo to\[CHacek]k, ki jih \[ZHacek]elim uporabiti za izra\[CHacek]un*)*)
(*sttock=Input["Naklju\[CHacek]no \[SHacek]tevilo to\[CHacek]k:"];*)
(**)
(*(*Shranim koordinate s funkcijo mccpi*)*)
(*{ktkvad,ktkrog}=mccpi[sttock];*)
(*(*Pre\[SHacek]tejem to\[CHacek]ke, ki le\[ZHacek]ijo znotraj kroga*)*)
(*stkrog=Count[Norm/@ktkvad,_?(#<1&)];*)
(**)
(*(*Izra\[CHacek]unam pribli\[ZHacek]no vrednost \[Pi], ter njeno napako*)*)
(*pipribl=N[4 stkrog/sttock];*)
(*napaka=pipribl-Pi;*)
(**)
(*(*Izpi\[SHacek]em*)*)
(*Print["Priblji\[ZHacek]na vrednost \[Pi]: ",pipribl ];*)
(*Print["Napaka \[Pi]: ",napaka ];*)
(**)
(*(*Vse prika\[ZHacek]em na enem grafu s funkcijo Show*)*)
(*p1=ListPlot[ktkvad,PlotLegends->{"To\[CHacek]ke zunaj kroga"},PlotStyle->Thick];*)
(*p2=ListPlot[ktkrog,PlotLegends->{"To\[CHacek]ke znotraj kroga"},PlotStyle->Red,PlotStyle->Thick];*)
(*p3=ParametricPlot[{Sin[fi],Cos[fi]},{fi,0,2Pi},PlotStyle->Black];*)
(*Show[p1,p2,p3,PlotLabel->"Monte Carlo",AspectRatio->1,AxesLabel->{"x-koordinate","y-koordinate"},PlotRange->{{-1,1},{-1,1}}]*)
(*]*)


(* ::Input:: *)
(*(*Funkcijo calcpi pokli\[CHacek]em, za \[SHacek]tevilo to\[CHacek]k, ki jih vpi\[SHacek]em v vhodno okno.*)
(*Funkija izpi\[SHacek]e vrednosti in izri\[SHacek]e graf*)*)
(*calcpi[sttock]*)


(* ::Section:: *)
(*1.2 Razvoj v vrsto in funkcija Manipulate*)


(* ::Input:: *)
(*(*S funkcijo manipulate upravljamo red aproksimacijskih Taylorjevih vrst, s katerimi aproksimiram funkcijo f[t]. Na grafu izri\[SHacek]em pravo funkcijo ter njeno aproksimacijo, kjer spreminjam red Taylorjeve vrste.*)*)
(*Manipulate[Module[{f,f1},*)
(*(*Funkcija f[t]*)*)
(*f[t_]:=Sin[t]*t^2*Exp[-t];*)
(*(*Approksimacija funkcije f[t]-Taylorjeva vrsta razvita okoli to\[CHacek]ke t->2*)*)
(*f1=Normal[Series[f[t],{t,2,n}]];*)
(*(*Graf obeh funkcij*)*)
(*Plot[{f[t],f1},{t,0,4},PlotRange->All,PlotLegends->{"f[t]",Row[{"Taylorjeva aproksimacija reda: ",n}]}]],{{n,3,"Red aproksimirane Taylorjeve vrste: "},1,10,1}]*)
(**)


(* ::Text:: *)
(*Kar lahko opazimo  je, da se z vi\[SHacek]anjem reda vse bolj pribli\[ZHacek]ujemo pravi funkciji f[t]!*)
