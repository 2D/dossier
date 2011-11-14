import rhinoscriptsyntax as rs, Rhino
import cmath, math
def zeta_ep(z,g2,g3):
    #Weierstrass epileptic function
    # http://mathworld.wolfram.com/WeierstrassEllipticFunction.html
    k=20
    c=generate_weier_coeffs(g2,g3,k)
    s=0
    for j in range(2,max(c)):
        s=c[j]*z**(2*j-2)+s
    return z**(-2)+s
def zeta_f(z,g2,g3):
    #Weierstrass zeta function
    # http://mathworld.wolfram.com/WeierstrassZetaFunction.html
    k=20
    c=generate_weier_coeffs(g2,g3,k)
    s=0
    for j in range(2,max(c)):
        s=c[j]*z**(2*j-1)/(2*j-1)+s
    return z**(-1)-s
def generate_weier_coeffs(g2,g3,k):
    #smotrim koefficienty zdes':
    # http://mathworld.wolfram.com/WeierstrassEllipticFunction.html
    c={}
    #c=zeros(7)
    c[2]=g2/20
    c[3]=g3/28
    for p in range(4,k+1):
        su=0
        for l in range(2,p):
            su=su+3/((2*p+1)*(p-1))*c[l]*c[p-2]
        c[p]=su
    return c
def x(u,v,R,g2,g3):
    # http://mathworld.wolfram.com/CostaMinimalSurface.html
    i=cmath.sqrt(-1)#mnimaya edinica
    pi=math.pi
    e1=6.87519
    return 0.5*R*(-zeta_f(u+i*v,g2,g3)+pi*u+pi**2/(4*e1)+pi/(2*e1)*
                  (zeta_f(u+i*v-0.5,g2,g3)-zeta_f(u+i*v-0.5*i,g2,g3)))
def y(u,v,R,g2,g3):
    i=cmath.sqrt(-1)#mnimaya edinica
    pi=math.pi#3.14159265 etc
    e1=6.87519
    return 0.5*R*(-i*zeta_f(u+i*v,g2,g3)+pi*u+pi**2/(4*e1)-pi/(2*e1)*
                  (zeta_f(u+i*v-0.5,g2,g3)-zeta_f(u+i*v-0.5*i,g2,g3)))
def z(u,v,R,g2,g3):
    i=cmath.sqrt(-1)#mnimaya edinica
    pi=math.pi#3.14159265 etc
    e1=6.87519
    return 1/4*math.sqrt(2*pi)*math.log(abs((zeta_ep(u+v*i,g2,g3)-e1)/(zeta_ep(u+v*i,g2,g3)+e1)))

if( __name__ == "__main__" ):
    g2=189.0727
    g3=0
    R=10#radius
    
    pts=[]
    
    umin = 0.001;umax = math.pi*2
    vmin = 0.001;vmax = math.pi*2
    v_N=u_N=10
    arrCrv = []
    for i in range(u_N):
        points=[]
        for j in range(v_N):
            u = (umin + i) * (umax - umin) / float(u_N)
            v = (vmin + j) * (vmax - vmin) / float(v_N)
            pt = [
                (x(u,v,R,g2,g3)),
                (y(u,v,R,g2,g3)),
                (z(u,v,R,g2,g3))
                ]
            pts.append(pt)
            print "pt:",pt
    #rs.AddSrfControlPtGrid( (len(u),len(v)), pts )
#print pts
    
    """
    Complex z1 = z1(a, b, n, k1);
         Complex z2 = z2(a, b, n, k2);


         float x_new = (float)(scal * z1.getRe());
         float y_new = (float)(scal * z2.getRe());
         float z_new = (float)(scal * (
         Math.cos(Alpha) * z1.getIm() +
           Math.sin(Alpha) * z2.getIm()));

         xPoints.add(x_new);
         yPoints.add(y_new);
         zPoints.add(z_new);

    Надеюсь, что называл функции понятно и тут все ясно написал.
    
    Alpha - угол камеры
    z1.getRe() - функция взятия действит. части (класс комплексного числа
    описан в другом фале (в процессинге должны появиться вкладки)).
    scal - просто масштаб.
    
    z1(a, b, n, k1) :
    a - меняется в диапазоне [-1, 1];
    b - [0, 1.7] (на самом деле до П/2);
    n - степень для вычисления числа (она не меняется, воспринимай, как параметр);
    k1 - k1 = [0; n - 1];
    возвращает комплексное число.
    
    """
