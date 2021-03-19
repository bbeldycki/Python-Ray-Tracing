from integrals import integrals
import math
import time


def carlson_intermediate_values_quartic_case() -> list:
    # Function which return a list of some intermediate values for given input
    # f1 = 2.7, g1 = -1.8, h1 = 0.9, f2= 2.0, g2 = 2.4, h2 = 0.8, a5 = 1.1 b5 = -0.4, x = 2.0, y = -3.0
    m2 = 0.36362947     # index 0
    lm2 = 0.53423014    # index 1
    lp2 = 24.673029     # index 2
    sig2 = 21.199185    # index 3
    sig02 = 6.1236295   # index 4
    a2 = 11237.193      # index 5
    b2 = 9741.4746      # index 6
    a02 = 844.71933     # index 7
    b02 = 247.52253     # index 8
    t2 = 10.288757      # index 9
    v2 = 1.1362990      # index 10
    t02 = 1.1328716     # index 11
    v02 = 1.1077327     # index 12
    ix = -1.1571677     # index 13
    ix0 = -0.093427949      # index 14
    rca2b2 = 0.00098889795  # index 15
    rca02b02 = 0.050085175  # index 16
    rct2v2 = 0.58372845     # index 17
    rct02v02 = 0.94657139   # index 18
    rfm2lm2lp2 = 0.54784092       # index 19
    rdm2lm2lp2 = 0.0422910488     # index 20
    rjm2lm2lp2sig2 = 0.048599080  # index 21
    rjm2lm2lp2sig02 = 0.12739513  # index 22
    gie = 10.495586         # index 23
    ha = 0.049905556        # index 24
    ha0 = -1.8557835        # index 25
    am111m1 = 1.5731367     # index 26
    a1111 = -0.49594737     # index 27
    am111m1m2 = 6.2622360   # index 28
    a1111m2 = 14.845682     # index 28
    result_list = [m2, lm2, lp2, sig2, sig02, a2, b2, a02, b02, t2, v2, t02, v02, ix, ix0, rca2b2, rca02b02, rct2v2,
                   rct02v02, rfm2lm2lp2, rdm2lm2lp2, rjm2lm2lp2sig2, rjm2lm2lp2sig02, gie, ha, ha0, am111m1, a1111,
                   am111m1m2, a1111m2]
    return result_list


def compute_delta_ij(f: list, g: list, h: list, k: int, l: int) -> float:
    # Function which computes coefficients delta_ij
    return 2 * f[k] * h[l] + 2 * f[l] * h[k] - g[k] * g[l]


def quartic_cases_test_data() -> tuple:
    a_list = [0.0, 0.0, 0.0, 0.0, 1.1]
    b_list = [0.0, 0.0, 0.0, 0.0, -0.4]
    # a_list = [0.0, 0.0, 0.0, 0.0, 1.0]
    # b_list = [0.0, 0.0, 0.0, 0.0, 0.0]
    fgh_1_list = [2.7, -1.8, 0.9]
    fgh_2_list = [2.0, 2.4, 0.8]
    y_val = -3.0
    x_val = 2.0
    p_list = [[-1, -1, -1, -1, 0], [-1, 1, 1, -1, 0], [1, 1, 1, 1, 0], [-1, 1, 1, -1, -2], [1, 1, 1, 1, -2]]
    # p_list is a list of lists
    # zero could be skipped but I will leave it to make sure every list has the same length
    # by default we will use p_list[0] = [ -1, -1, -1, 0, 0]
    # [-1, 1, 1, -1, 0]     index 1
    # [1, 1, 1, 1, 0]       index 2
    # [-1, 1, 1, -1, -2]    index 3
    # [1, 1, 1, 1, -2]      index 4
    return p_list, a_list, b_list, fgh_1_list, fgh_2_list, y_val, x_val


def comparison_result_lists(list1: list, list2: list, eps: float) -> bool:
    # Function to compare my computed values with values from Carlson computations
    # Remember to put my computed values in list1 and Carlson values in list2
    if len(list1) != len(list2):
        raise Exception('Error: Len of list1 is not equal to Len of list2.')
    for i in range(len(list1)):
        if abs(list1[i] - list2[i]) > eps:
            raise Exception('Error occurred during comparison of ours computed results with values from carlson '
                            'article. Reason: Absolute value of difference is larger than one part in a million. '
                            f'Error occurred for i value {i}, with my result = {list1[i]} and carlson = {list2[i]}')
    return True


def test_integrals_quartic_case_m2lm2lp2(fgh1: list, fgh2: list, y: float, x: float) -> list:
    if x < y:
        raise Exception('Error occurred in function test_integrals_quartic_case_m2lm2lp2. '
                        'We have situation where x is lower than y. The value of x must be greater than y')
    fi = [fgh1[0], fgh2[0]]
    gi = [fgh1[1], fgh2[1]]
    hi = [fgh1[2], fgh2[2]]
    ksi = []
    eta = []
    # coefficient ksi and eta from equation 2.1
    for i in range(len(fi)):
        ksi.append(math.sqrt(fi[i] + gi[i] * x + hi[i] * x * x))
        eta.append(math.sqrt(fi[i] + gi[i] * y + hi[i] * y * y))
    dzeta = []
    # coefficient dzeta from equation 2.5
    for i in range(len(fi)):
        dzeta.append(math.sqrt((ksi[i] + eta[i]) ** 2 - hi[i] * (x - y) ** 2))
    # coefficient M from equation 2.6
    m = dzeta[0] * dzeta[1] / (x - y)
    m2 = m * m
    # coefficients delta_ij and delta from equation 2.7
    d112 = compute_delta_ij(fi, gi, hi, 0, 0)
    d222 = compute_delta_ij(fi, gi, hi, 1, 1)
    d122 = compute_delta_ij(fi, gi, hi, 0, 1)
    delta = math.sqrt(d122 * d122 - d112 * d222)
    # coefficients delta_plus/minus and Lm2 and Lp2 from equation 2.8
    deltap = d122 + delta
    deltam = d122 - delta
    lm2 = m2 + deltam
    lp2 = m2 + deltap
    return [m2, lm2, lp2]


if __name__ == "__main__":
    start = time.time()
    p_test, a_test, b_test, fgh1_test, fgh2_test, y_test, x_test = quartic_cases_test_data()
    data_to_compare_with = carlson_intermediate_values_quartic_case()
    results_m2lm2lp2 = test_integrals_quartic_case_m2lm2lp2(fgh1_test, fgh2_test, y_test, x_test)
    eps = 1e-6
    if comparison_result_lists(results_m2lm2lp2, [data_to_compare_with[0],
                                                  data_to_compare_with[1], data_to_compare_with[2]], eps=eps):
        print(f'Correct values of: M2, Lm2, Lp2 with eps = {eps}')
    koniec = time.time()
    time = (koniec - start) * 1000
    print(f'Tests took: {time} ms')