"""
This script generates the C++ code necessary to numerically evaluate the 
residual and Jacobian of the singular potential inversion.

First, it writes the code to calculate each of the degree 4 and degree 2
monomial integ
"""
import argparse

# the kth index of these vectors represents the row and column position
# respectively of the kth degree of freedom of traceless, symmetric tensors
i_idx = [0, 1, 0, 0, 1]
j_idx = [0, 1, 1, 2, 2]

def get_commandline_args():

    desc = '''Generates C++ code necessary to numerically evaluate residual 
    and Jacobian for singular potential inversion.'''

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--dim',
                        choices=['full_3D', 'quasi_2D'])
    args = parser.parse_args()

    return args.dim


def make_monomial_list(degree):
    """
    Create a list of all monomials of degree `degree`.
    Monomials are stored as a list of three integers, corresponding to the
    degrees of the x-, y-, and z-coordinates.

    Parameters
    ----------
    degree : int
        degree of each of the monomials in the list

    Returns
    -------
    monomial_list : list of list of ints
        Each item in the `monomial_list` list is a list of three ints, 
        representing the degrees of the x-, y-, and z-coordinates.
    """

    monomial_list = []

    for i in range(degree + 1):
        for j in range(i, degree + 1):
            monomial_list.append([degree - j, j - i, i])

    return monomial_list



def integral_expressions(idx, monomial_list):
    """
    Print expression for calculating integrals of monomials multiplied by
    exponential of Lambda contracted with coordinates.

    Parameters
    ----------
    idx : int
        Index of integral expression (based on the order of the monomials)
    monomial_list : list of lists of ints
        List of degree 2 or degree 4 monomials from `make_monomial_list`

    Returns
    -------
    int_expression : string
        C++ code which calculates the integral expression corresponding to the
        index `idx` and `monomial_list`
    """

    monomial = monomial_list[idx]
    degree = monomial[0] + monomial[1] + monomial[2]

    if degree == 2:
        return 'I2[{}] += {}{}{}exp_lambda * w[q];'.format(idx,
                                                           'x[q] * ' * monomial[0],
                                                           'y[q] * ' * monomial[1],
                                                           'z[q] * ' * monomial[2])
    elif degree == 4:
        return 'I4[{}] += {}{}{}exp_lambda * w[q];'.format(idx,
                                                           'x[q] * ' * monomial[0],
                                                           'y[q] * ' * monomial[1],
                                                           'z[q] * ' * monomial[2])
    else:
        raise ValueError('Integral expression degree must be 2 or 4')



def residual_expression(m, monomial_list_2):
    """
    Given an index `m` which denotes the degree of freedom, and an ordered
    list of second-degree monomials, prints the code to calculate the residual

    Parameters
    ----------
    m : int
        Indexes the degree of freedom of the residual to calculate
    monomial_list_2 : list of lists of ints
        Ordered list of monomials of degree 2 as calculated from 
        `make_monomial_list`

    Returns
    -------
    residual_code : string
        String which is the C++ code used to calculate the mth element of the
        residual
    """

    monomial = [0, 0, 0]
    i = i_idx[m]
    j = j_idx[m]

    monomial[i] += 1
    monomial[j] += 1

    integral_idx = monomial_list_2.index( monomial )

    return 'Res({}) = 1 / Z * I2[{}] - m({});'.format(m, integral_idx, m)



def jacobian_expression(m, n, monomial_list_4, monomial_list_2):
    """
    Given indices `m` and `n` which denotes (respectively) the degrees of 
    freedom of the residual and of the element of Lambda which the derivative 
    is being taken with respect to, this function prints the code to calculate 
    the Jacobian

    Parameters
    ----------
    m : int
        Indexes the degree of freedom of the residual
    n : int
        Indexes the degree of freedom of Lambda which the derivative is being
        taken with respect to
    monomial_list_4 : list of lists of ints
        Ordered list of monomials of degree 4 as calculated from 
        `make_monomial_list`
    monomial_list_2 : list of lists of ints
        Ordered list of monomials of degree 2 as calculated from 
        `make_monomial_list`

    Returns
    -------
    jacobian_code : string
        String which is the C++ code used to calculate the (m, n)th element of 
        the jacobian
    """

    monomial_4 = [0, 0, 0]
    monomial_2_1 = [0, 0, 0]
    monomial_2_2 = [0, 0, 0]

    i_m = i_idx[m]
    j_m = j_idx[m]

    i_n = i_idx[n]
    j_n = j_idx[n]

    monomial_4[i_m] += 1
    monomial_4[j_m] += 1
    monomial_4[i_n] += 1
    monomial_4[j_n] += 1

    monomial_2_1[i_m] += 1
    monomial_2_1[j_m] += 1
    monomial_2_2[i_n] += 1
    monomial_2_2[j_n] += 1

    I4_idx = monomial_list_4.index( monomial_4 )
    I21_idx = monomial_list_2.index( monomial_2_1 )
    I22_idx = monomial_list_2.index( monomial_2_2 )

    element_is_diagonal = n < 2
    if not element_is_diagonal:
        return 'Jac({}, {}) = 2 / Z * (I4[{}] - 1 / Z * I2[{}]*I2[{}]);'.format(m, n, I4_idx, I21_idx, I22_idx)

    monomial_4_diag = [0, 0, 2]
    monomial_4_diag[i_m] += 1
    monomial_4_diag[j_m] += 1
    monomial_2_diag = [0, 0, 2]

    I4_idx_diag = monomial_list_4.index(monomial_4_diag)
    I2_idx_diag = monomial_list_2.index(monomial_2_diag)

    diag_string = 'Jac({}, {}) = 1 / Z * (I4[{}] - I4[{}] - 1 / Z * I2[{}]*(I2[{}] - I2[{}]));'

    return diag_string.format(m, n, I4_idx, I4_idx_diag, I21_idx, I22_idx, I2_idx_diag)



def main():

    dim = get_commandline_args()
    vec_dim = 5 if dim == 'full_3D' else 3

    monomial_list_2 = make_monomial_list(2)
    monomial_list_4 = make_monomial_list(4)

    if dim == 'quasi_2D':
        allowed_z_powers = lambda coord_power: (coord_power[2] == 2) or (coord_power[2] == 0)
        monomial_list_2 = list( filter(allowed_z_powers, monomial_list_2) )
        monomial_list_4 = list( filter(allowed_z_powers, monomial_list_4) )

    # Print unique monomials
    print('Ordering of monomials is:', end='\n\n')
    for coord_power in monomial_list_2:
        print( '{}{}{}'.format('x'*coord_power[0],
                               'y'*coord_power[1],
                               'z'*coord_power[2]) )
    for coord_power in monomial_list_4:
        print( '{}{}{}'.format('x'*coord_power[0],
                               'y'*coord_power[1],
                               'z'*coord_power[2]) )
    print(end='\n\n')

    # Print code to calculate integral expressions
    print('Monomial integral code is:', end='\n\n')
    for i in range(len(monomial_list_2)):
        print(integral_expressions(i, monomial_list_2))
    for i in range(len(monomial_list_4)):
        print(integral_expressions(i, monomial_list_4))
    print(end='\n\n')

    # Print code to calculate residual and Jacobian
    print('Residual code is:', end='\n\n')
    for i in range(vec_dim):
        print( residual_expression(i, monomial_list_2) )
    for i in range(vec_dim):
        for j in range(vec_dim):
            print( jacobian_expression(i, j, monomial_list_4, monomial_list_2) )



if __name__ == '__main__':
    main()
