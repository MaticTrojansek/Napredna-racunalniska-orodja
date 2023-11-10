(* ::Package:: *)

mccpi[sttock_]:=Module[{ktkvad,ktkrog},
ktkvad=RandomReal[{-1,1},{sttock,2}];
ktkrog=Select[ktkvad,Norm[#]<=1&];
{ktkvad,ktkrog}]
