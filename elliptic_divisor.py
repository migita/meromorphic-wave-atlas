"""Exact torsion divisors on Y²=4X³-g2 X-g3.

A small standalone implementation of the chord-and-tangent group law and
Miller's line/vertical-divisor construction, sufficient for the atlas. It
checks the supplied point's exact order and returns the normalized function
with divisor order*(point)-order*(infinity).
"""
import sympy as s

X,Y=s.symbols("X Y")


def torsion_function(g2,g3,point,order,domain=s.QQ):
    if order<2:raise ValueError("The order must be at least two")
    if s.cancel(g2**3-27*g3**2)==0:raise ValueError("Use a nonsingular generic curve")
    P=tuple(map(s.sympify,point));px,py=P
    cubic=4*X**3-g2*X-g3
    if s.cancel(py**2-4*px**3+g2*px+g3)!=0:raise ValueError("Point is off the curve")

    def nf(poly):
        return s.Poly(s.expand(poly),Y).rem(s.Poly(Y**2-cubic,Y)).as_expr().expand()

    def step(A,B):
        ax,ay=A;bx,by=B
        if s.cancel(ax-bx)==0 and s.cancel(ay+by)==0:
            return None,X-ax
        if s.cancel(ax-bx)==0:
            slope=s.cancel((3*ax**2-g2/4)/ay)
        else:
            slope=s.cancel((by-ay)/(2*(bx-ax)))
        cx=s.cancel(slope**2-ax-bx)
        cy=s.cancel(-ay+2*slope*(ax-cx))
        line=Y/2-ay/2-slope*(X-ax)
        return (cx,cy),line/(X-cx)

    current=P;function=s.Integer(1);certificate=[]
    for multiple in range(2,order+1):
        following,factor=step(current,P)
        certificate.append(dict(multiple=multiple,point=following,factor=factor))
        function*=factor
        if following is None and multiple!=order:raise ValueError("The point has smaller order")
        current=following
    if current is not None:raise ValueError("The specified multiple is not infinity")
    numerator,denominator=s.fraction(s.together((-1)**(order-2)*function))
    numerator,denominator=nf(numerator),nf(denominator)
    A,B=numerator.coeff(Y,0),numerator.coeff(Y,1)
    C,D=denominator.coeff(Y,0),denominator.coeff(Y,1)
    norm=s.Poly(s.expand(C*C-D*D*cubic),X,domain=domain)
    real=s.Poly(s.expand(A*C-B*D*cubic),X,domain=domain).exquo(norm)
    imaginary=s.Poly(s.expand(B*C-A*D),X,domain=domain).exquo(norm)
    result=s.expand(real.as_expr()+Y*imaginary.as_expr())
    # The conjugate norm is an independent divisor/normalization check.
    expected=(-1)**order*(X-px)**order
    actual=nf(result*result.subs(Y,-Y))
    if s.cancel(actual-expected)!=0:raise AssertionError("Incorrect torsion norm")
    return result,certificate
