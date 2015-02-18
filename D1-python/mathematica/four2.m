Domains[n_, opts___] :=
  Module[
    {tos, Ts, Ttot, Fs, Ftot,
      centered, param, timedomain, frequencydomain},

  Ts = DeltaT /. {opts} /. Options[Domains];
  Ttot = TotalT /. {opts} /. Options[Domains];
  Fs = DeltaF /. {opts} /. Options[Domains];
  Ftot = TotalF /. {opts} /. Options[Domains];
  tos = Toffset /. {opts} /. Options[Domains];
  centered = Centered /. {opts} /. Options[Domains];
  param = OnlyParameters /. {opts} /. Options[Domains];

  Switch[
      First[Flatten[Position[N[{Ts, Ttot, Fs, Ftot, tos}], _?NumberQ]]],
      1, Ttot = n Ts; Fs = 1/Ttot; Ftot = 1/Ts,
      2, Ts = Ttot/n; Fs = 1/Ttot; Ftot = 1/Ts,
      3, Ts = 1/(n Fs); Ttot = 1/Fs; Ftot = n Fs,
      4, Ts = 1/Ftot; Ttot = n Ts; Fs = Ftot/n,
      _, Ts = 1; Ttot = n; Fs = 1/n; Ftot = 1
  ];

  If[param,
      Return[{Ts, Ttot, Fs, Ftot}],
      timedomain = tos + (Range[n] - 1)Ts;
      frequencydomain =
          If[centered,
              (Range[n] - Ceiling[n/2])Fs,
              (Range[n] - 1)Fs
          ];
      Return[{timedomain, frequencydomain}]
  ]
]

Options[AddSignalRange] = {Toffset -> 0};

AddSignalRange[data_List, opts___] :=
  Transpose[{Domains[Length[data], opts][[1]], data}]

Options[AddSpectrumRange] = {
      Centered -> False,
      HalfSpectrum -> False
      };

AddSpectrumRange[data_List, opts___] :=
  Module[
      {n, lista, half, centered},

      half = HalfSpectrum /. {opts} /. Options[AddSpectrumRange];
      If[half,
          centered = False,
          centered = Centered /. {opts} /. Options[AddSpectrumRange];
      ];

      n = Length[data];
      lista = If[centered,

              Transpose[
                          {Domains[n, Centered -> True, opts][[2]],
                              RotateRight[data, Ceiling[n/2] - 1]}
              ],

              Transpose[
                          {Domains[n, Centered -> False, opts][[2]],
                              data}
              ]   
      ];
      If[half, (*then centered is forced to False*)
                      Take[lista, {1, Ceiling[(n + 1)/2]}],
                      lista (*can be centered or not *)
  ]
]

n = 150; Fs = 30;
data = Table[.7Sin[2Pi 4 t] - .4Cos[2Pi 12 t] + .2(Random[] - .5),
        {t, 0, (n - 1)/Fs, 1/Fs}];

plot1 = ListLinePlot[data, PlotRange -> Full]
signal = AddSignalRange[data, TotalF -> Fs];
plot2 = ListLinePlot[signal, PlotRange -> Full]
mag = Abs[Fourier[data]];
spectrum = AddSpectrumRange[mag, TotalF -> Fs, Centered -> True];
(* spectrum = AddSpectrumRange[mag, TotalF -> Fs, HalfSpectrum -> True];  *)
plot3 = ListLinePlot[spectrum, PlotRange -> Full]



Export["p1.png", plot1]
Export["p2.png", plot2]
Export["p3.png", plot3]

