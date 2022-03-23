func checklog ( ticklen , all=) {
  if (is_void(ticklen)) ticklen=0.006;
  if (is_void(logAdjMino_horiz)) logAdjMino_horiz=4;
  if (is_void(logAdjMino_vert)) logAdjMino_vert=4;
  
  if (!all)
    N = plsys();
  else N = [];
  isys = plsys();
  
  get_style,  landscape , sys, legends, clegends;

  sys(N).ticks.horiz.tickStyle.type = 1;
  sys(N).ticks.horiz.tickLen = [0.009, ticklen , ticklen ,0.0002,0.0002];
  //sys(N).ticks.horiz.tickLen = sys(N).ticks.horiz.tickLen *[1, 3, 1, 1, 0];
  sys(N).ticks.horiz.logAdjMinor= 100./sys(N).ticks.horiz.nMinor;
  //  sys(N).ticks.horiz.nMinor=1000;
  //sys(N).ticks.horiz.nDigits= 3;
  
  //sys(N).ticks.vert.tickLen = sys(N).ticks.vert.tickLen *  [1, 1, 0, 0, 0];
  sys(N).ticks.vert.tickLen = [0.009, ticklen , ticklen ,0.0002,0.0002];
  sys(N).ticks.vert.logAdjMinor= 100./sys(N).ticks.vert.nMinor;
  set_style,  landscape , sys , legends, clegends;

  plsys, isys;
}
  
