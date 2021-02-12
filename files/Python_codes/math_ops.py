# cos
# sin
# sqrt

E = 2.718281828459045
E10 = 22026.465794806703
E100 = 2.6881171418161212e+43
PI = 3.141592653589793

def exp(x):
    if x == 0:
        return 1
    if x > 0:
        sign = True
    else:
        sign = False
        x = -x
    y = 1
    while x > 100:
        y *= E100
        x -= 100
    while x > 10:
        y *= E10
        x -= 10
    while x > 1:
        y *= E
        x -= 1
    while x > 0.1:
        y *= 1.1051709180756477
        x -= 0.1
    y *= (1.+x+x*x/2.+x*x*x/6.+x*x*x*x/24.+x*x*x*x*x/120+x*x*x*x*x*x/720\
        +x*x*x*x*x*x*x/5040.+x*x*x*x*x*x*x*x/40320.)
    if not sign:
        y = 1/y
    return  y


def ln(x):
    if x < 0:
        return None
    if x == 0:
        return -float("inf")
    if x < 1:
        sign = False
        x = 1/x
    else:
        sign = True

    y_int = 0
    while x > E100:
        y_int += 100.
        x /= E100
    while x > E10:
        y_int += 10.
        x /= E10
    while x > E:
        y_int += 1.
        x /= E 
    y = 1
    y = y-1+x*exp(-y)
    y = y-1+x*exp(-y)
    y = y-1+x*exp(-y)
    y = y-1+x*exp(-y)
    y = y-1+x*exp(-y)

    y += y_int
    if not sign:
        y = -1*y
    return y

# def ln(x):
#     if x < 0:
#         return None
#     if x == 0:
#         return -float("inf")
#     if x == 1:
#         return 0
#     if x > 1:
#         if x > 20000:
#             l, r = 0.0, 709.0
#         else:
#             l, r = 0.0, 10.0
#     else:
#         if x < 1.0/20000:
#             l, r = -709.0, 0.0
#         else:
#             l, r = -10.0, 0.0
#     y = 0.0
#     while r-l > 1:
#         y = (r+l)/2.0
#         if exp(y) >= x:
#             r = y
#         else:
#             l = y
#     y = y-1+x*exp(-y)
#     y = y-1+x*exp(-y)
#     y = y-1+x*exp(-y)
#     y = y-1+x*exp(-y)
#     y = y-1+x*exp(-y)
#     return y


def power(x, y):
    return exp(ln(x)*y)


def ceil(x):
    x_int = int(x)
    if abs(x-x_int) < 1.e-15:
        return x_int
    if x > 0:
        return x_int+1
    else:
        return x_int


def floor(x):
    x_int = int(x)
    if abs(x-x_int) < 1.e-15:
        return x_int
    if x > 0:
        return x_int
    else:
        return x_int-1