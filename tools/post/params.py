from sympy.parsing.sympy_parser import parse_expr

model_params = [('X Problem Size', 128), \
                ('Y Problem Size', 128), \
                ('Z Problem Size', 128), \
                ('X Block Size', 16), \
                ('Y Block Size', 16), \
                ('Z Block Size', 128), \
                ('Num Ghosts (NG)', 4), \
                ('Num Species', 9), \
                ('GB/s/thread', 1), \
                ('Gflop/s/thread', 10), \
                ('Int Regs/thread', 16), \
                ('FP Regs/thread', 16), \
                ('Threads', 1024), \
                ('$/thread group (kB)', 1024), \
                ('Word Size', 8), \
                ('Cache Line Size', 64), \
                ( 'R cost', 8), \
                ('RW cost', 16), \
                ( 'W cost', 16), \
                ('Division Cost', 39), \
                ('Special Cost', 125), \
                ('Streaming Writes', False), \
                ('NTA Hints', False), \
                ('NIC BW (GB/s)', 100), \
                ('NIC Latency (us)', 0.4), \
                ('CPU/DRAM Overlap', 1), \
                ('PWFlag', False)]

def defaultParams():
  params = {}
  for (name, value) in model_params:
    params[name] = value
  return params

def customizeParams(params, overrides):
  for (key, value) in overrides.iteritems():
    if key not in params:
      raise Exception('parameter: %s not valid' % key)
    if (value != params[key]):
      print 'Overriding parameter %s = %s' % (key, value)
      params[key] = value

# symbolic replacments
symrepl = [('dlo(1)', 'lo(1) - ng'), \
           ('dhi(1)', 'hi(1) + ng'), \
           ('dlo(2)', 'lo(2) - ng'), \
           ('dhi(2)', 'hi(2) + ng'), \
           ('dlo(3)', 'lo(3) - ng'), \
           ('dhi(3)', 'hi(3) + ng'), \
           ('hi(1) - lo(1) + 1', '__BoxX__'), \
           ('hi(2) - lo(2) + 1', '__BoxY__'), \
           ('hi(3) - lo(3) + 1', '__BoxZ__'), \
           ('ncons', 'nspecies + 5'), \
           ('nprim', '3 * nspecies + 7'), \
           ('iry1+n-1', 'iryn'), \
           ('qh1+n-1', 'qhn'), \
           ('qx1+n-1', 'qxn'), \
           ('qy1+n-1', 'qyn'), \
          ]

def doSymSubs(expr):
  if type(expr) == type(''):
    expr = parse_expr(expr)
  for r in symrepl:
    expr = expr.subs(r[0], r[1])
  return expr

# parameter replacements
paramrepl = [('__BoxX__', 'X Problem Size'), \
             ('__BoxY__', 'Y Problem Size'), \
             ('__BoxZ__', 'Z Problem Size'), \
             ('__BlockX__', 'X Block Size'), \
             ('__BlockY__', 'Y Block Size'), \
             ('__BlockZ__', 'Z Block Size'), \
             ('ng', 'Num Ghosts (NG)'), \
             ('nspecies', 'Num Species'), \
            ]

# non-parameter replacements
nonparamrepl = [('__LoopBlockX__', 'blockxcl'), \
                ('__LoopBlockY__', 'blocky'), \
                ('__LoopBlockZ__', 'blockz')
               ]

def doRefSubs_l(expr, lam=lambda x: x):
  if type(expr) != type(''):
    expr = str(expr)
  for r in paramrepl + nonparamrepl:
    expr = expr.replace(r[0], lam(r[1]))
  return '(%s)' % expr

def doParamSubs(expr, params, nonparams = {}):
  if type(expr) != type(0) and type(expr) != type(0.):
    for r in paramrepl:
      expr = expr.subs(r[0], params[r[1]])
    for r in nonparamrepl:
      if r[1] in nonparams:
        expr = expr.subs(r[0], nonparams[r[1]])
  return float(expr)

specComponents = ['n', 'qxn', 'qyn', 'qhn', 'iryn', 'idXn']

def isSpeciesArray(x):
  return x.split('.')[-1] in specComponents
