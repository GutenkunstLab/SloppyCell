#!/usr/bin/env python
"""
     Symbolic.

         Symbolic differentiation of python expressions.

Copyright 1999 Pearu Peterson, <pearu@ioc.ee>
March 11-12, 1999
Pearu Peterson

Modified by Ryan Gutenkunst, <rng7@cornell.edu>
April 28th, 2005
General clean-up of code.
Also changed saving mechanism to explicitly save the results, rather than doing it each time a derivative is taken.
"""

__version__ = "0.2"

import parser
import symbol,token
import os, cPickle

__Diff_saved = {}

def loadDiffs(filename):
    global __Diff_saved
    try:
        __Diff_saved = cPickle.load(file(filename, 'rb'))
    except IOError:
        print 'saved diffs not loaded'

def saveDiffs(filename):
    f = file(filename, 'wb')
    cPickle.dump(__Diff_saved, f, 2)
    f.close()

syms=token.tok_name
for s in symbol.sym_name.keys():
    syms[s]=symbol.sym_name[s]
    
def ast2shortlist(t):
    if type(t) is parser.ASTType: return ast2shortlist(t.tolist())
    if not isinstance(t, list): return t
    if t[1] == '': return None
    if not isinstance(t[1], list): return t
    if len(t) == 2 and isinstance(t[1], list):
        return ast2shortlist(t[1])
    o=[]
    for tt in map(ast2shortlist, t[1:]):
        if tt is not None:
            o.append(tt)
    if len(o)==1: return o[0]
    return [t[0]]+o

def sym2name(t):
    if type(t) is parser.ASTType: return sym2name(t.tolist())
    if not isinstance(t, list): return t
    return [syms[t[0]]]+map(sym2name,t[1:])

def string2ast(t):
    return sym2name(ast2shortlist(parser.expr(t)))

def ast2string(t):
    #if isinstance(t, str): return t
    if type(t) is parser.ASTType: return ast2string(t.tolist())
    if not isinstance(t, list): return None
    if not isinstance(t[1], list): return t[1]
    o=''
    for tt in map(ast2string,t):
        if isinstance(tt, str):
            o=o+tt
    return o

def getlftrt(t,a):
    lft=t[1]
    rt=t[3:]
    if len(rt)>1:
        rt=[t[0]]+rt
        drt=Diff(rt,a)
        if drt[0] not in ['NUMBER']:
            drt=['atom',['LPAR','('],drt,['RPAR',')']]
    else:
        rt=rt[0]
        drt=Diff(rt,a)
    dlft=Diff(lft,a)
    return lft,dlft,rt,drt

def trysimple(t):
    try: t = str(eval(t))
    except: pass
    return t

def ensureparen(t,flag=0):
    t=trysimple(t)
    if t[0]=='-':
        return '(%s)'%t
    tt=string2ast(t)
    if tt[0]=='arith_expr':
        return '(%s)'%t
    if flag>0 and tt[0] in ['term','power']:
        return '(%s)'%t
    return t

def dopower(l,r):
    if r=='0': return '1'
    if l=='0': return '0'
    if l=='1': return '1'
    if r=='1': return (l)
    return (trysimple('%s**%s'%(ensureparen(l),ensureparen(r))))

def domul(l,r):
    if l=='0' or r=='0': return '0'
    if l=='1': return (r)
    if r=='1': return (l)
    return (trysimple('%s*%s'%(ensureparen(l),ensureparen(r))))

def dodiv(l,r):
    if l=='0': return '0'
    if r=='1': return (l)
    if r=='-1': return (doneg(l))
    if r==l: return '1'
    return (trysimple('%s/%s'%(ensureparen(l),ensureparen(r,1))))

def doadd(l,r):
    if l=='0' and r=='0': return '0'
    if l=='0': return r
    if r=='0': return l
    if l==r: return trysimple(domul('2',l))
    if r[0]=='-': return trysimple('%s%s'%(l,r))
    return ensureparen(trysimple('%s+%s'%(l,r)))

def dosub(l,r):
    if l=='0' and r=='0': return '0'
    if l=='0': return doneg(r)
    if r=='0': return l
    if l==r: return '0'
    if r[0]=='-': return doadd(l,r[1:])
    return ensureparen(trysimple('%s-%s'%(l,r)))

def doneg(l):
    if l=='0': return '0'
    if l[0]=='-': return l[1:]
    return trysimple('-%s'%l)

def dofun(l,r):
    if l=='log':
        if r=='1': return '0'
        if r=='e': return '1'
    if l=='log10':
        if r=='1': return '0'
        if r=='10': return '1'
    if l=='log_0': return dodiv('1',r)
    if l=='log10_0': return dodiv(dofun('log','10'),r)
    if l=='ln_0': return dodiv('1',r)
    if l=='sin_0': return dofun('cos',r)
    if l=='sinh_0': return dofun('cosh',r)
    if l=='cosh_0': return dofun('sinh',r)
    if l=='tanh_0': return dodiv('1',dopower(dofun('cosh',r),'2'))
    if l=='coth_0': return dodiv('-1',dopower(dofun('sinh',r),'2'))
    if l=='asin_0': return dodiv('1',dofun('sqrt',dosub('1',dopower(r,'2'))))
    if l=='acos_0': return dodiv('-1',dofun('sqrt',dosub('1',dopower(r,'2'))))
    if l=='atan_0': return dodiv('1',doadd('1',dopower(r,'2')))
    if l=='acot_0': return dodiv('-1',doadd('1',dopower(r,'2')))
    if l=='tan_0': return dodiv('1',dopower(dofun('cos',r),'2'))
    if l=='cot_0': return doneg(dodiv('1',dopower(dofun('sin',r),'2')))
    if l=='cos_0': return doneg(dofun('sin',r))
    if l=='exp_0': return dofun('exp',r)
    if l=='sqrt_0': return dodiv('1',domul('2',dofun('sqrt',r)))
    return '%s(%s)'%(l,r)

def splitargs(da):
    if da.find('(')<0:
        return da.split(',')
    ll=[];o='';ii=0
    for i in da:
        if i==',' and ii==0:
            ll.append(o)
            o=''
        else:
            if i=='(': ii=ii+1
            if i==')': ii=ii-1
            o=o+i
    ll.append(o)
    return ll

def Diff(t, a='x'):
    if isinstance(a, list):
        s= ''.join([t, '__derivWRT__', str(a)])
        r=__Diff_saved.get(s, None)
        if r is not None: 
            return r

        o = string2ast(t)
        for aa in a:
            o = Diff(o,aa)
        res=ast2string(o)
        __Diff_saved[s] = res
        return res
    elif isinstance(t, str):
        s= ''.join([t, '__derivWRT__', str(a)])
        r=__Diff_saved.get(s, None)
        if r is not None: 
            return r

        res=ast2string(Diff(string2ast(t),a))
        __Diff_saved[s] = res
        return res

    if isinstance(t, list) and len(t) == 1:
        return Diff(t[0],a)
    if t[0]=='NUMBER': return [t[0],'0']
    if t[0]=='NAME':
        if t[1]==a: return ['NUMBER','1']
        return ['NUMBER','0']
    if t[0]=='factor':
        st=ast2string(t)
        return string2ast(doneg(ast2string(Diff(string2ast(st[1:]),a))))
    if t[0]=='arith_expr':
        o='0'
        for i in range(1,len(t),2):
            if t[i-1][1]=='-':
                o=dosub(o,ast2string(Diff(t[i],a)))
            else:
                o=doadd(o,ast2string(Diff(t[i],a)))
        return string2ast(o)
    if t[0]=='term':
        lft,dlft,rt,drt=map(ast2string,getlftrt(t,a))
        if t[2][0]=='STAR':
            return string2ast(doadd(domul(dlft,rt),domul(lft,drt)))
        if t[2][0]=='SLASH':
            return string2ast(dosub(dodiv(dlft,rt),domul(domul(lft,drt),dopower(rt,'-2'))))
    if t[0]=='power':
        if t[2][0]=='trailer':
            if len(t)>3:
                return string2ast(Diff(ast2string([t[0],string2ast('(%s)'%ast2string(t[:3])),t[3:]]),a))
        ts=[];o=[];dts=[]
        for i in t[1:]:
            if i[0] in ['DOUBLESTAR']:
                if len(o)==1: o=o[0]
                ts.append(o);
                dts.append(Diff(o,a))
                o=[]
            else: o.append(i)
        if len(o)==1: o=o[0]
        ts.append(o)
        dts.append(Diff(o,a))
        if t[2][0]=='DOUBLESTAR':
            st,lft,dlft,rt,drt=map(ast2string,[t,ts[0],dts[0],ts[1],dts[1]])
            rt1=trysimple('%s-1'%rt)
            if drt=='0':
                return string2ast(domul(domul(rt,dlft),dopower(lft,rt1)))
            return string2ast(domul(doadd(domul(dofun('log',lft),drt),dodiv(domul(dlft,rt),lft)),st))
        if t[2][0]=='trailer':
            return dts[0]
    if t[0] in ['arglist','testlist']:
        o=[]
        for i in ast2string(t).split(','): o.append(Diff(i,a))
        return string2ast(','.join(o))
    if t[0]=='atom':
        return string2ast(ensureparen(ast2string(Diff(t[2:-1],a))))
    if t[1][0]=='trailer': # t=[[NAME,f],[trailer,[(],[ll],[)]]]
        aa=t[1][2]
        ll=splitargs(ast2string(Diff(aa,a)))
        ii=0;o='0'
        for i in ll:
            o=doadd(o,domul(i,dofun(t[0][1]+'_'+`ii`,ast2string(aa))))
            ii=ii+1
        return string2ast(o)
    return t

def test1():
    print 'in test1:', __Diff_saved
    x='x'
    for str in ['x','x','a*x','a','a(x)','a(x)*x','a(x)*b(x)',
                'a+x','a(x)+x','a(x)+b(x)+c(x)','a(x)*b(x)*c(x)',
                'a(x)*(b(x)+c(x))','(a(x)+b(x))*c(x)',
                'a(x)/b(x)','(a(x))**b(x)',
                'x**n','x**5','x**-5','sin(x)','cos(x)','exp(x)',
                'ln(x)','log(x)','log10(x)','asin(x)','sinh(x)']:
        dstr=Diff(str,x)
        print 'Diff(%s,%s)=%s'%(str,x,dstr)

def test3():
    str='x**-1'
    for i in range(3):
        dstr=Diff(str,['x'])
        print 'Diff(%s,%s)=%s'%(str,'x',dstr)
        str=dstr

if __name__ == "__main__":
    loadDiffs('temp')
    test1()
    test3()
    print __Diff_saved
    saveDiffs('temp')
