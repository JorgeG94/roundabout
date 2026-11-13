"""Generates the fitted coefficient tables of src/shared_module_utilities/rdb_safe_math_consts.inc
(atan, erf, erfc and the cbrt seed) and prints them as Fortran on stdout.

MAINTAINER TOOL, not part of the build or the tests: needs mpmath, which the dev
python3.10 provides (run it with python3.10; do not pip install anything).  Run
`python3.10 tools/gen_safe_math_coeffs.py > tmp_local_artifacts/coeffs.inc` (about two
minutes, deterministic) and paste the result under the "fitted tables" line of the
constants file.  The fit report goes to stderr.

Method.  Polynomials: Chebyshev interpolation at Chebyshev nodes (near-minimax),
re-expanded in the natural variable.  erfc: R(x) = exp(x**2) erfc(x) on [0.5, 2] and
[2, 4], and F(w) = (x R sqrt(pi) - 1)/w, w = 1/x**2, on x >= 4, each as a (m,m) rational
function: Chebyshev-node rational interpolation, denominator rounded to double, numerator
refitted by least squares against the rounded denominator, then a seeded hill-climb of the
rounded coefficients against the error measured in TRUE double arithmetic (the
coefficients are only as good as Horner's rounding allows: about 4e-16 relative).
Method references: Cody & Waite (1980), Software Manual for the Elementary Functions;
Cody (1969), Math. Comp. 23, 631 (the erfc range split); Muller (2016), Elementary
Functions, ch. 3.
"""
import sys
from mpmath import mp, mpf, cos, pi, sqrt, erfc, erf, exp, atan, nstr, matrix, lu_solve
mp.prec = 200
def dbl(x):
    m,e=mp.frexp(x); return mp.nint(m*2**53)*mpf(2)**(e-53)
def ev(co,s):
    p=mpf(0)
    for v in reversed(co): p=p*s+v
    return p
def cheb_poly(f,a,b,n):
    """degree n-1 interpolant on [a,b]; returns monomial coefficients in x itself."""
    c=(a+b)/2;h=(b-a)/2
    nodes=[cos((2*j+1)*pi/(2*n)) for j in range(n)]
    vals=[f(c+h*s) for s in nodes]
    cc=[sum(vals[j]*cos(k*(2*j+1)*pi/(2*n)) for j in range(n))*2/n for k in range(n)]
    cc[0]/=2
    T=[[mpf(1)],[mpf(0),mpf(1)]]
    for k in range(2,n):
        new=[mpf(0)]+[2*v for v in T[k-1]]
        for i,v in enumerate(T[k-2]): new[i]-=v
        T.append(new)
    mono_s=[mpf(0)]*n
    for k in range(n):
        for i,v in enumerate(T[k]): mono_s[i]+=cc[k]*v
    # s = (x-c)/h  ->  monomial in x
    out=[mpf(0)]*n
    base=[mpf(1)]
    lin=[-c/h,1/h]
    for i in range(n):
        for j,v in enumerate(base): out[j]+=mono_s[i]*v
        nb=[mpf(0)]*(len(base)+1)
        for j,v in enumerate(base):
            nb[j]+=v*lin[0]; nb[j+1]+=v*lin[1]
        base=nb
    return out
def rat(f,a,b,m,n):
    c=(a+b)/2;h=(b-a)/2;N=m+n+1
    nodes=[cos((2*j+1)*pi/(2*N)) for j in range(N)]
    A=matrix(N,N);rhs=matrix(N,1)
    for j,s in enumerate(nodes):
        fx=f(c+h*s)
        for i in range(m+1): A[j,i]=s**i
        for i in range(1,n+1): A[j,m+i]=-fx*s**i
        rhs[j]=fx
    sol=lu_solve(A,rhs)
    return c,h,[sol[i] for i in range(m+1)],[mpf(1)]+[sol[m+i] for i in range(1,n+1)]
def fmt(name,vals,lo):
    items=[nstr(dbl(v),19)+"_wp" for v in vals]
    s=f"   real(wp), parameter :: {name}({lo}:{lo+len(vals)-1}) = [ &\n"
    lines=[]
    for i in range(0,len(items),3):
        lines.append("      "+", ".join(items[i:i+3]))
    return s+", &\n".join(lines)+" ]\n"
def polytable(name,g,a,b,deg,target=3e-17):
    co=[dbl(v) for v in cheb_poly(g,a,b,deg+1)]
    mx=max(abs(ev(co,a+(b-a)*i/600)/g(a+(b-a)*i/600)-1) for i in range(601))
    sys.stderr.write(f"{name}: deg {deg} max rel {nstr(mx,3)}\n")
    return fmt(name,co,0)
UMAX=(2-sqrt(3))**2
out=[]
out.append(polytable("SM_ATAN_C",lambda u: atan(sqrt(u))/sqrt(u) if u>0 else mpf(1),mpf(0),dbl(UMAX)*(1+mpf(2)**-50),11))
out.append(polytable("SM_ERF_C",lambda v: erf(sqrt(v))/sqrt(v) if v>0 else 2/sqrt(pi),mpf(0),mpf(1)/4,11))
from mpmath import cbrt
out.append(polytable("SM_CBRT_C",lambda x: cbrt(x),mpf("0.5"),mpf(4),4))
sp=sqrt(pi)
Rx=lambda x: exp(x*x)*erfc(x)
Tx=lambda x: x*exp(x*x)*erfc(x)*sp
def refitP(F,a,b,m,Q):
    N=2*(m+1)
    nodes=[(a+b)/2+(b-a)/2*cos((2*j+1)*pi/(2*N)) for j in range(N)]
    A=matrix(N,m+1);rhs=matrix(N,1)
    for j,t in enumerate(nodes):
        qv=ev(Q,t); w=1/F(t)
        for i in range(m+1): A[j,i]=t**i*w
        rhs[j]=qv*F(t)*w
    from mpmath import qr_solve
    sol=qr_solve(A,rhs)[0]
    return [sol[i] for i in range(m+1)]
def ratrun(name,F,a,b,m):
    N=2*m+1
    nodes=[(a+b)/2+(b-a)/2*cos((2*j+1)*pi/(2*N)) for j in range(N)]
    A=matrix(N,N);rhs=matrix(N,1)
    for j,t in enumerate(nodes):
        fx=F(t)
        for i in range(m+1): A[j,i]=t**i
        for i in range(1,m+1): A[j,m+i]=-fx*t**i
        rhs[j]=fx
    sol=lu_solve(A,rhs)
    Qd=[mpf(1)]+[dbl(sol[m+i]) for i in range(1,m+1)]
    Pd=[dbl(v) for v in refitP(F,a,b,m,Qd)]
    def dev(co,t):
        p=0.0
        for v in reversed(co): p=p*t+float(v)
        return p
    pts=[a+(b-a)*i/600 for i in range(601)]
    ts=[float(t) for t in pts]; ref=[F(t) for t in pts]
    def worst(P,Q):
        return max(abs(float(dev(P,t)/dev(Q,t)/r-1)) for t,r in zip(ts,ref))
    mx=worst(Pd,Qd)
    # polish: hill-climb the rounded coefficients (+-1..3 ulp) against the true double evaluation
    import random
    rnd=random.Random(12345)
    P=[float(v) for v in Pd]; Q=[float(v) for v in Qd]
    import math
    for it in range(6000):
        which=rnd.random()<0.5
        arr=P if which else Q
        j=rnd.randrange(len(arr)) if which else rnd.randrange(1,len(arr))
        old=arr[j]
        arr[j]=old+rnd.choice([-3,-2,-1,1,2,3])*math.ulp(old)
        w=worst(P,Q)
        if w<mx: mx=w
        else: arr[j]=old
    Pd=[mpf(v) for v in P]; Qd=[mpf(v) for v in Q]
    sys.stderr.write(f"{name}: deg ({m},{m}) on [{nstr(a,5)},{nstr(b,5)}] max rel (double eval) {mx:.2e}\n")
    return fmt(name+"_P",Pd,0)+fmt(name+"_Q",Qd,0)
out.append(ratrun("SM_ERFC_A",Rx,mpf('0.5'),mpf(2),6))
out.append(ratrun("SM_ERFC_B",Rx,mpf(2),mpf(4),6))
Fw=lambda w: (Tx(1/sqrt(w))-1)/w
out.append(ratrun("SM_ERFC_C",Fw,1/mpf('27.3')**2,mpf(1)/16,5))
print("".join(out))
