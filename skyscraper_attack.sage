
p_BLS381 = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001

Fp = GF(p_BLS381)


def circular_shift(x, k, n):
    """
    Shift circularly the bits of an n-bit word of k positions towards the most significant bit.
    """
    return ((x  << k) | (x >> (n-k))) & (2**n-1)

def bitwise_not(x, n):
    """
    Compute the bitwise not of an n-bit word.
    """
    return x ^^ (2**n -1)

def chi8(x):
    """
    Compute the 8-bit Daemen chi function.
    """
    return x ^^ (circular_shift(bitwise_not(x, 8), 1, 8) & circular_shift(x, 2, 8) & circular_shift(x, 3, 8))

def L(x):
    """
    Compute the lookups in the SkyScraper hash function.
    """
    return circular_shift(chi8(x), 1, 8)


def Bar(x):
    """
    The Bar function defined in the article.
    """
    x_int = ZZ(x)
    
    x_list = []
    
    for i in range(32):
        x_list.append(x_int % 256)
        x_int //= 256
    y_list = [L(xi) for xi in x_list]
    y_list = y_list[16:] + y_list[:16]
    
    y = Fp(0)
    for i in range(32):
        y += 2**(8*i)*y_list[i]
    return y

def SkyScraper(xL, xR, rcons):
    """
    The full SkyScraper permutation.
    """
    yL, yR = xL, xR
    
    for i in range(2):
        yL, yR = yR + yL**2 + rcons[i], yL
    for i in range(2):
        yL, yR = yR + Bar(yL) + rcons[i+2], yL
    for i in range(2):
        yL, yR = yR + yL**2 + rcons[i+4], yL
    for i in range(2):
        yL, yR = yR + Bar(yL) + rcons[i+6], yL
    for i in range(2):
        yL, yR = yR + yL**2 + rcons[i+8], yL
    return yL, yR

def inbound_phase(delta_in_L, delta_in_R, delta_out_L, delta_out_R):
    """
    Find values in the middle to get the differential you need in input and output.
    """
    delta_R = delta_in_L
    delta_L = delta_out_R
    
    x_L =  (delta_out_L  - delta_L**2 - delta_R)    / (2*delta_L)
    x_R =  (delta_L      - delta_R**2 - delta_in_R) / (2*delta_R)
    
    return (x_L, x_R, delta_L, delta_R)

def select_delta_out(delta):
    """
    Select a delta_out which will propagate with high probability in outbound phase.
    """
    delta_out_L = delta*(2**129)
    delta_out_R = - 5*delta
    
    return delta_out_L, delta_out_R

def select_delta_in(delta):
    """
    Select a delta_in which will propagate with high probability in outbound phase.
    """
    delta_in_L = 5*delta
    delta_in_R = delta*(2**129)
    
    return delta_in_L, delta_in_R

def get_input_pair(delta_in, delta_out, rcons):
    """
    Given delta_in, delta_out, compute the input pair.
    """
    delta_in_L, delta_in_R = select_delta_in(delta_in)
    delta_out_L, delta_out_R = select_delta_out(delta_out)
    (x_L, x_R, delta_L, delta_R) = inbound_phase(delta_in_L, delta_in_R, delta_out_L, delta_out_R)
    
    x_L1 = x_L
    x_L2 = x_L + delta_L
    x_R1 = x_R
    x_R2 = x_R + delta_R
    
    x_L1, x_R1 = x_R1, x_L1
    x_L2, x_R2 = x_R2, x_L2

    x_L1, x_R1 = x_R1 - (x_L1**2 + rcons[4]), x_L1
    x_L2, x_R2 = x_R2 - (x_L2**2 + rcons[4]), x_L2
    
    for i in range(2):
        x_L1, x_R1 = x_R1 - (Bar(x_L1) + rcons[3-i]), x_L1
        x_L2, x_R2 = x_R2 - (Bar(x_L2) + rcons[3-i]), x_L2
    
    for i in range(2):
        x_L1, x_R1 = x_R1 - (x_L1**2 + rcons[1-i]), x_L1
        x_L2, x_R2 = x_R2 - (x_L2**2 + rcons[1-i]), x_L2
    
    x_L1, x_R1 = x_R1, x_L1
    x_L2, x_R2 = x_R2, x_L2

    return [[x_L1, x_R1], [x_L2, x_R2]]



def trunc_diff(rcons, good_delta):
    """
    Find input-output pairs mathing some truncated differential characteristic.
    """
    for delta_in in good_delta:
        for delta_out in good_delta:
            for i in range(15):
                for j in range(15):
                    [[x_L1, x_R1], [x_L2, x_R2]] = get_input_pair(delta_in*2**(8*i), delta_out*2**(8*j), rcons)
                    y_L1, y_R1 = SkyScraper(x_L1, x_R1, rcons)
                    y_L2, y_R2 = SkyScraper(x_L2, x_R2, rcons)

                    if (x_L2 - x_L1 == delta_in*2**(8*i)) and (y_R2 - y_R1 == -1*delta_out*2**(8*j)):
                        print("Input 1 is:")
                        print(hex(x_L1), hex(x_R1))
                        print("Input 2 is:")
                        print(hex(x_L2), hex(x_R2))
                        print("Output 1 is:")
                        print(hex(y_L1), hex(y_R1))
                        print("Output 2 is:")
                        print(hex(y_L2), hex(y_R2))
                        print("Differential: ({}, ?) -> (?, {})".format(hex(x_L2 - x_L1), hex(y_R2 - y_R1)))


def collision_9(rcons, good_delta):
    """
    Find input-output pairs mathing some truncated differential characteristic.
    """
    for delta in good_delta:
        for i in range(15):
            [[x_L1, x_R1], [x_L2, x_R2]] = get_input_pair(delta*2**(8*i), delta*2**(8*i), rcons)
            y_L1, y_R1 = SkyScraper(x_L1, x_R1, rcons)
            y_L2, y_R2 = SkyScraper(x_L2, x_R2, rcons)

            if (x_L2 - x_L1 + y_R2 - y_R1 == 0):
                # return 1 if a collision was found for this choice of constants
                return 1
    return 0


def test_trunc_diff():
    # Use round constans from the reference implementation
    rcons = [0, 17829420340877239108687448009732280677191990375576158938221412342251481978692, 
         27740342931201890067831390843279536630457710544396725670188095857896839417202, 17048088173265532689680903955395019356591870902241717143279822196003888806966, 
         4641041932484616674503150018277105937647322787229235415313079154504538247098, 23518768991468467328187394347260979305359711922005254253047385842741274989784, 
         42924498470449697215909973597747708755064028547754583139137172884417685452938, 4670171540012394890944659920922397025152994631853985835833059662854559397332, 
         16971509144034029782226530622087626979814683266929655790026304723118124142299, 0]

    print("Finding truncated differentials for one instance:")
    trunc_diff(rcons, good_delta)

def test_collision_9(n_iterations):
    score = 0

    print("Looking for collisions on the 9-round compression:")

    for k in range(n_iterations):
        rcons = [Fp.random_element() for i in range(10)]
        score += collision_9(rcons, good_delta)

    proba = float(score / n_iterations)

    print("A collision was found with probability {} over the choice of constants.".format(proba))
    print("Over {} instances of SkyScraper.".format(n_iterations))


if __name__ == "__main__":
    good_delta = [1, 2, 4, 8, 16, -1, -2, -4, -8, -16]
    good_delta = [Fp(delta) for delta in good_delta]

    n_iterations = 1000

    test_trunc_diff()
    test_collision_9(n_iterations)

