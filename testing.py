from LossToleranceFunctions.Losstolerant_fusion import*
from LossToleranceFunctions.Erasure_threshold_codes import*
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline




def XZ(w, p_fail, transmission):
    # XY, ZX, ZY
    P_Z = p_fail * (1 - w) * (transmission ** (1 / p_fail))
    P_X = p_fail * w * (transmission ** (1 / p_fail))
    return P_Z, P_X




def adaptive_coeff_and_weigth(coeff_list):
    # P_S = 1, P_X = 2, P_Z = 3, P_L = 4
    # P_X_1, P_X_123,P_X_1231
    weight_dict = {}
    for cf in coeff_list:
        qbit = []
        for key in cf.keys():
            seq_numb = ''.join(str(x) for x in qbit)
            if seq_numb not in weight_dict.keys():
                weight_dict[seq_numb] = [0, 0]
            if cf[key] == "P_S":
                qbit.append(1)
            elif cf[key] == "P_X":
                weight_dict[seq_numb][1] += 1
                qbit.append(2)
            elif cf[key] == "P_Z":
                weight_dict[seq_numb][0] += 1
                qbit.append(3)
            elif cf[key] == "P_L":
                qbit.append(4)
        weight_dict.pop('', None)
    return weight_dict

def adaptive_coeff_with_Y(coeff_list):
    # P_S = 1, P_X = 2, P_Z = 3, P_L = 4
    # P_X_1, P_X_123,P_X_1231
    new_coeff_dict_full = []
    sequence_numbers = []
    for cf in coeff_list:
        new_coeff_dict = {}
        qbit = []
        for key in cf.keys():
            if cf[key] == "P_S":
                seq_numb = ''.join(str(x) for x in qbit)
                new_coeff_dict[key] = "P_S_" + seq_numb
                if seq_numb not in sequence_numbers:
                    sequence_numbers.append(seq_numb)
                qbit.append(1)
            elif cf[key] == "P_X":
                seq_numb = ''.join(str(x) for x in qbit)
                new_coeff_dict[key] = "P_X_" + seq_numb
                if seq_numb not in sequence_numbers:
                    sequence_numbers.append(seq_numb)
                qbit.append(2)
            elif cf[key] == "P_Z":
                seq_numb = ''.join(str(x) for x in qbit)
                new_coeff_dict[key] = "P_Z_" + seq_numb
                if seq_numb not in sequence_numbers:
                    sequence_numbers.append(seq_numb)
                qbit.append(3)
            elif cf[key] == "P_L":
                seq_numb = ''.join(str(x) for x in qbit)
                new_coeff_dict[key] = "P_L_" + seq_numb
                if seq_numb not in sequence_numbers:
                    sequence_numbers.append(seq_numb)
                qbit.append(4)
            elif cf[key] == "P_Y":
                seq_numb = ''.join(str(x) for x in qbit)
                new_coeff_dict[key] = "P_Y_" + seq_numb
                if seq_numb not in sequence_numbers:
                    sequence_numbers.append(seq_numb)
                qbit.append(5)
        new_coeff_dict_full.append(new_coeff_dict)
    return new_coeff_dict_full, sequence_numbers

def adaptive_coeff(coeff_list):
    # P_S = 1, P_X = 2, P_Z = 3, P_L = 4
    # P_X_1, P_X_123,P_X_1231
    new_coeff_dict_full = []
    sequence_numbers = []
    for cf in coeff_list:
        new_coeff_dict = {}
        qbit = []
        for key in cf.keys():
            if cf[key] == "P_S":
                seq_numb = ''.join(str(x) for x in qbit)
                new_coeff_dict[key] = "P_S_" + seq_numb
                if seq_numb not in sequence_numbers:
                    sequence_numbers.append(seq_numb)
                qbit.append(1)
            elif cf[key] == "P_X":
                seq_numb = ''.join(str(x) for x in qbit)
                new_coeff_dict[key] = "P_X_" + seq_numb
                if seq_numb not in sequence_numbers:
                    sequence_numbers.append(seq_numb)
                qbit.append(2)
            elif cf[key] == "P_Z":
                seq_numb = ''.join(str(x) for x in qbit)
                new_coeff_dict[key] = "P_Z_" + seq_numb
                if seq_numb not in sequence_numbers:
                    sequence_numbers.append(seq_numb)
                qbit.append(3)
            elif cf[key] == "P_L":
                seq_numb = ''.join(str(x) for x in qbit)
                new_coeff_dict[key] = "P_L_" + seq_numb
                if seq_numb not in sequence_numbers:
                    sequence_numbers.append(seq_numb)
                qbit.append(4)
        new_coeff_dict_full.append(new_coeff_dict)
    return new_coeff_dict_full, sequence_numbers

def turn_adaptive_coeff_to_weigth(adapt_coeff_in, numb_qubits = 4):
    tabs_list = {}
    for item in adapt_coeff_in:
        for i in range(1, numb_qubits):
            if len(item) > i:
                tabs_list[item[:i]] = [0, 0]
    # adapt_coeff = adapt_coeff_in
    for i in range(1, numb_qubits):
        # print("i {}".format(i))
        adapt_coeff = adapt_coeff_in
        for item in adapt_coeff:
            if len(item) == i:
                # adapt_coeff.remove(item)
                for item_1 in adapt_coeff:
                    if len(item_1) > i and item[:i] == item_1[:i]:
                        # print("adapt {} {} {} {}".format(item_1, i, item, adapt_coeff))
                        # if item[:i] == item_1[:i]:
                        if item_1[i] == "3":
                            tabs_list[str(item_1[:i])][0] += 1
                        elif item_1[i] == "2":
                            tabs_list[str(item_1[:i])][1] += 1
    #print(tabs_list)
    #print(adapt_coeff_in)
    return tabs_list


def turn_adaptive_coeff_to_weigth_with_Y(adapt_coeff_in, numb_qubits = 4):
    tabs_list = {}
    for item in adapt_coeff_in:
        for i in range(1, numb_qubits):
            if len(item) > i:
                tabs_list[item[:i]] = [0, 0, 0]
    # adapt_coeff = adapt_coeff_in
    for i in range(1, numb_qubits):
        # print("i {}".format(i))
        adapt_coeff = adapt_coeff_in
        for item in adapt_coeff:
            if len(item) == i:
                # adapt_coeff.remove(item)
                for item_1 in adapt_coeff:
                    if len(item_1) > i and item[:i] == item_1[:i]:
                        # print("adapt {} {} {} {}".format(item_1, i, item, adapt_coeff))
                        # if item[:i] == item_1[:i]:
                        if item_1[i] == "3":
                            tabs_list[str(item_1[:i])][0] += 1
                        elif item_1[i] == "2":
                            tabs_list[str(item_1[:i])][1] += 1
                        elif item_1[i] == "5":
                            tabs_list[str(item_1[:i])][2] += 1
    #print(tabs_list)
    #print(adapt_coeff_in)
    return tabs_list



def extract_wt_from_tab_list(tab_list):
    wt_dict = {}
    for key in tab_list.keys():
        Z,X = tab_list[key]
        # wt = X / (X + Z)
        if Z > X:
            wt = 0
        elif Z < X:
            wt = 1
        elif X == Z:
            wt = 1 / 2
        wt_dict[key] = wt
    return wt_dict


def extract_last_meas(coeff, numb_qubit):
    XZ_meas = []
    for pattern in coeff:
        initial = pattern["qubit_" + str(numb_qubit)][:3]
        if initial == "P_X":
            XZ_meas.append(pattern["qubit_" + str(numb_qubit)][4:] + str(2))
        elif initial == "P_Z":
            XZ_meas.append(pattern["qubit_" + str(numb_qubit)][4:] + str(3))
    return XZ_meas


def extract_last_meas_with_Y(coeff, numb_qubit):
    XZ_meas = []
    for pattern in coeff:
        initial = pattern["qubit_" + str(numb_qubit)][:3]
        if initial == "P_X":
            XZ_meas.append(pattern["qubit_" + str(numb_qubit)][4:] + str(2))
        elif initial == "P_Z":
            XZ_meas.append(pattern["qubit_" + str(numb_qubit)][4:] + str(3))
        elif initial == "P_Y":
            XZ_meas.append(pattern["qubit_" + str(numb_qubit)][4:] + str(5))
    return XZ_meas


def extract_first_weigth(coeff, numb_qubit=1):
    XZ_meas = {"X": 0,"Z": 0}
    for pattern in coeff:
        initial = pattern["qubit_" + str(numb_qubit)][:3]
        if initial == "P_X":
            XZ_meas["X"] += 1
        elif initial == "P_Z":
            XZ_meas["Z"] += 1
    if XZ_meas["X"] > XZ_meas["Z"]:
        wt = 1
    elif XZ_meas["X"] < XZ_meas["Z"]:
        wt = 0
    elif XZ_meas["X"] == XZ_meas["Z"]:
        wt = 1 / 2
    return wt


def extract_first_weigth_coeff_erasure(coeff, numb_qubit=1):
    XZ_meas = {"X": 0,"Z": 0}
    for pattern in coeff:
        initial = pattern["qubit_" + str(numb_qubit)][:3]
        if initial == "P_X":
            XZ_meas["X"] += 1
        elif initial == "P_Z":
            XZ_meas["Z"] += 1
    return XZ_meas


def extract_first_weigth_coeff_erasure_new(coeff, numb_qubit=1):
    XZ_meas = {"X": 0,"Z": 0}
    for i in range(len(coeff[''])):
        if i == 0:
            XZ_meas["Z"] += 1
        else:
            XZ_meas["X"] += 1
    return XZ_meas

def extract_first_weigth_coeff_erasure_with_Y(coeff, numb_qubit = 1):
    XZ_meas = {"X": 0, "Z": 0, "Y": 0}
    for pattern in coeff:
        initial = pattern["qubit_" + str(numb_qubit)][:3]
        if initial == "P_X":
            XZ_meas["X"] += 1
        elif initial == "P_Z":
            XZ_meas["Z"] += 1
        elif initial == "P_Y":
            XZ_meas["Y"] += 1
    return XZ_meas


def extract_first_weigth_erasure_with_Y(XZ_1, XZ_2):
    X = 0
    Z = 0
    Y = 0
    for key in XZ_1.keys():
        if key == "X":
            X = XZ_1[key] + XZ_2[key]
            # X_1, X_2 = XZ_1[key], XZ_2[key]
        elif key == "Z":
            Z = XZ_1[key] + XZ_2[key]
        elif key == "Y":
            Z = XZ_1[key] + XZ_2[key]
            # Z_1, Z_2 = XZ_1[key], XZ_2[key]
    if Z >= X:
        one = Z
        if Y >= X:
            two = Y
            meas_bases = "35"
        else:
            two = X
            meas_bases = "32"
    else:
        one = X
        if Y >= Z:
            two = Y
            meas_bases = "25"
        else:
            two = Z
            meas_bases = "32"
    if one > two:
        wt = 0
    elif one < two:
        wt = 1
    elif one == two:
        wt = 1 / 2
    return wt, meas_bases



def extract_first_weigth_erasure(XZ_1, XZ_2):
    X = 0
    Z = 0
    print(XZ_1, XZ_2)
    for key in XZ_1.keys():
        if key == "X":
            # X = XZ_1[key] + XZ_2[key]
            X_1, X_2 = XZ_1[key], XZ_2[key]
        else:
            # Z = XZ_1[key] + XZ_2[key]
            Z_1, Z_2 = XZ_1[key], XZ_2[key]
    if X_2 > Z_2:
        if X_1 > Z_1:
            wt = 1
        else:
            wt = 1 / 2
    else:
        if X_1 > Z_1:
            wt = 1 / 2
        else:
            wt = 0
    # if X > Z:
    #     wt = 1
    # elif X < Z:
    #     wt = 0
    # elif X == Z:
    #     wt = 1 / 2
    return wt


def convert_analystic_solution(coeff_list, wt_dict, wt_first,p_fail, transmission):
    P_S = (1 - p_fail) * (transmission ** (1 / p_fail))
    P_L = (1 - (transmission ** (1 / p_fail)))
    tot_term = 0
    for coeff_dict in coeff_list:
        term = 1
        for key in coeff_dict.keys():
            weigth_numb = coeff_dict[key][4:]
            if weigth_numb not in wt_dict.keys():
                weigth = wt_first
            else:
                weigth = wt_dict[weigth_numb]
            P_Z, P_X = XZ(weigth, p_fail, transmission)
            if coeff_dict[key][:3] == "P_S":
                term = term * (P_S)
            elif coeff_dict[key][:3] == "P_X":
                term = term * (P_X)
            elif coeff_dict[key][:3] == "P_Z":
                term = term * (P_Z)
            else:
                term = term * (P_L)
        tot_term += term
    return tot_term


def convert_analystic_solution_with_Y(coeff_list, wt_dict, p_fail, transmission, meas_bases, wt_first
                                      , first_meas_base):
    P_S = (1 - p_fail) * (transmission ** (1 / p_fail))
    P_L = (1 - (transmission ** (1 / p_fail)))
    tot_term = 0
    for coeff_dict in coeff_list:
        term = 1
        # print(coeff_dict)
        for key in coeff_dict.keys():
            weigth_numb = coeff_dict[key][4:]
            # print(weigth_numb)
            # print(meas_bases)
            if weigth_numb not in wt_dict.keys():
                weigth = wt_first
                meas_base = first_meas_base
            else:
                weigth = wt_dict[weigth_numb]
                meas_base = meas_bases[weigth_numb]
            # print(weigth_numb)
            # XY, ZX, ZY, Z=1, X=2 , Y = 3, 12, 13, 23
            # P_Z, P_X
            if meas_base == "32":
                P_1, P_2 = XZ(weigth, p_fail, transmission)
                meas_base_trans_1 = "Z"
                meas_base_trans_2 = "X"
            elif meas_base == "35":
                P_1, P_2 = XZ(weigth, p_fail, transmission)
                meas_base_trans_1 = "Z"
                meas_base_trans_2 = "Y"
            else:
                meas_base_trans_1 = "X"
                meas_base_trans_2 = "Y"
                P_1, P_2 = XZ(weigth, p_fail, transmission)
            # print("P_1, P_2 {} {}".format(P_1, P_2))
            # print(coeff_dict[key][:3])
            if coeff_dict[key][:3] == "P_S":
                term = term * (P_S)
            elif coeff_dict[key][:3] == "P_" + meas_base_trans_1:
                term = term * (P_1)
            elif coeff_dict[key][:3] == "P_" + meas_base_trans_2:
                term = term * (P_2)
            elif coeff_dict[key][:3] == "P_L":
                term = term * (P_L)
            else:
                # print("HERE")
                term = 0
        # if term != 0:
        #     print(coeff_dict)
        #     print(wt_dict[weigth_numb])
        #     print(weigth_numb)
        #     print(meas_bases[weigth_numb])
        tot_term += term
        # print("P_1 " + str(P_1) + " P_2 " + str(P_2))
        # print(meas_base_trans_1, meas_base_trans_2)
        # print("term" + str(term))
    return tot_term



def add_X_Z_wt_list_new(tab_list_X, tab_list_Z):
    final_dict = {}
    for key in tab_list_X.keys():
        final_dict[key] = tab_list_X[key]
    for key in tab_list_Z.keys():
        if key in final_dict.keys():
            first_X = final_dict[key][0]
            first_Z = final_dict[key][1]
            first_Y = final_dict[key][2]
            second_X = tab_list_Z[key][0]
            second_Z = tab_list_Z[key][1]
            second_Y = tab_list_Z[key][2]
            if min(first_X, second_X) > min(first_Z, second_Z) and min(first_X, second_X) > min(first_Y, second_Y):
                # wt = 1
                final_dict[key][0] = 1
                final_dict[key][1] = 0
                final_dict[key][2] = 0
            elif min(first_Z, second_Z) > min(first_X, second_X) and min(first_Z, second_Z) > min(first_Y, second_Y):
            # elif min(XX, ZX) < min(XZ, ZZ):
                # wt = 0
                final_dict[key][0] = 0
                final_dict[key][1] = 1
                final_dict[key][2] = 0
            elif min(first_Y, second_Y) > min(first_X, second_X) and min(first_Y, second_Y) > min(first_Z, second_Z):
            # elif min(XX, ZX) == min(XZ, ZZ):
                # wt = 1
                final_dict[key][0] = 0
                final_dict[key][1] = 0
                final_dict[key][1] = 1
            else:
                final_dict[key][0] = 1
                final_dict[key][1] = 0
                final_dict[key][1] = 0
        else:
            final_dict[key] = tab_list_Z[key]
    return final_dict



def add_X_Z_wt_list(tab_list_X, tab_list_Z):
    final_dict = {}
    for key in tab_list_X.keys():
        final_dict[key] = tab_list_X[key]
    for key in tab_list_Z.keys():
        if key in final_dict.keys():
            final_dict[key][0] += tab_list_Z[key][0]
            final_dict[key][1] += tab_list_Z[key][1]
        else:
            final_dict[key] = tab_list_Z[key]
    return final_dict





def convert_coeff(coeff_list, p_fail, w, transmission):
    P_X = p_fail * w * (transmission ** (1 / p_fail))
    P_Z = p_fail * (1 - w) * (transmission ** (1 / p_fail))
    P_s = (1 - p_fail) * (transmission ** (1 / p_fail))
    P_l = (1 - (transmission ** (1 / p_fail)))
    tot_term = 0
    for coeff_dict in coeff_list:
        term = 1
        for key in coeff_dict.keys():
            if key == "P_s":
                term = term * (P_s ** coeff_dict[key])
            elif key == "P_X":
                term = term * (P_X ** coeff_dict[key])
            elif key == "P_Z":
                term = term * (P_Z ** coeff_dict[key])
            else:
                term = term * (P_l ** coeff_dict[key])
        tot_term += term
    return tot_term




def extract_wt_from_tab_list_with_Y(tab_list):
    wt_dict = {}
    list_of_bases = {}
    for key in tab_list.keys():
        Z, X, Y = tab_list[key]
        one = 0
        two = 0
        if Z >= X:
            one = Z
            if Y >= X:
                two = Y
                meas_bases = "35"
            else:
                two = X
                meas_bases = "32"
        else:
            # one = X
            if Y >= Z:
                one = X
                two = Y
                meas_bases = "25"
            else:
                one = Z
                two = X
                # two = Z
                meas_bases = "32"
        list_of_bases[key] = meas_bases
        if one > two:
            wt = 0
        elif one < two:
            wt = 1
        elif one == two:
            wt = 1 / 2
        wt_dict[key] = wt
    return wt_dict, list_of_bases



def get_weights(coef_1, coef_2):
    adapt_dict_1, adpat_cof_1 = adaptive_coeff_with_Y(coef_1)
    adapt_dict_2, adpat_cof_2 = adaptive_coeff_with_Y(coef_2)
    last_qubits_1 = extract_last_meas_with_Y(adapt_dict_1, numb_qubits - 1)
    adpat_cof_final_1 = adpat_cof_1 + last_qubits_1
    last_qubits_2 = extract_last_meas_with_Y(adapt_dict_2, numb_qubits - 1)
    adpat_cof_final_2 = adpat_cof_2 + last_qubits_2
    tab_list_1 = turn_adaptive_coeff_to_weigth_with_Y(adpat_cof_final_1, numb_qubits)
    tab_list_2 = turn_adaptive_coeff_to_weigth_with_Y(adpat_cof_final_2, numb_qubits)
    # print(tab_list_1)
    # print(tab_list_2)
    summed_wt_list_12 = add_X_Z_wt_list(tab_list_1, tab_list_2)
    # print(summed_wt_list_12)
    wt_dict_12, meas_base_12 = extract_wt_from_tab_list_with_Y(summed_wt_list_12)
    first_wt_1 = extract_first_weigth_coeff_erasure_with_Y(adapt_dict_1)
    first_wt_2 = extract_first_weigth_coeff_erasure_with_Y(adapt_dict_2)
    first_wt, first_meas_base = extract_first_weigth_erasure_with_Y(first_wt_1, first_wt_2)
    return adapt_dict_1, adapt_dict_2, wt_dict_12, meas_base_12, first_wt, first_meas_base


if __name__ == '__main__':
    # TODO : LOOKING ONLY AT GRAPHS: 2 * N < numb edges < 2.5 * N

    #################### VARIABLES #######################
    in_qubit = 0
    p_fails = [1 / 2]  # [1 / 2, 1 / 4, 1 / 8, 1 / 16, 1 / 32]
    transvals = np.linspace(0.955, 1, 100)
    weigths = np.linspace(0.001, 0.999, 100)
    ring_threshold = 1 - 0.1198
    # ghz_threshold = 1 - 0.069
    numb_qubits = 10
    p_fail = 1 / 2
    trans_val = np.linspace(0.915, 0.905, 10)
    best_loss = 1
    best_graph = 0
    for i in range(20):
        print(i)
        #################### INITIALIZE #######################
        # graph = nx.Graph()
        # edges = [(0, 1), (1, 2), (2, 3), (0, 3), (1, 3)]
        # graph.add_edges_from(edges)
        # graph = gen_fullyconnected_graph(numb_qubits)
        # graph = gen_star_graph(numb_qubits)
        # graph = gen_random_connected_graph(numb_qubits)
        graph = gen_random_connected_fixed_numb_edges_graph(numb_qubits)
        gstate = GraphState(graph)
        decoder_X = AdaptiveErasureDecoderNew(gstate, "X", numb_qubits)
        decoder_Z = AdaptiveErasureDecoderNew(gstate, "Z", numb_qubits)
        coef_X, tab_list_X_new, adapt_dict_X_new = decoder_X.run_dec_updated()
        coef_Z, tab_list_Z_new, adapt_dict_Z_new = decoder_Z.run_dec_updated()
        print("decoding done")
        summed_wt_list_XZ = add_X_Z_wt_list_new(tab_list_X_new, tab_list_Z_new)
        wt_dict_XZ, meas_base_XZ = extract_wt_from_tab_list_with_Y(summed_wt_list_XZ)
        first_wt = 0
        first_meas_base = "32"
        for trans in trans_val:
            analytic_sol_Z_XZ = convert_analystic_solution_with_Y(adapt_dict_Z_new, wt_dict_XZ, 1 / 2, trans,
                                                                  meas_base_XZ, first_wt, first_meas_base)
            analytic_sol_X_XZ = convert_analystic_solution_with_Y(adapt_dict_X_new, wt_dict_XZ, 1 / 2, trans,
                                                                  meas_base_XZ, first_wt, first_meas_base)
            print(analytic_sol_Z_XZ, analytic_sol_X_XZ)
            if min(analytic_sol_X_XZ, analytic_sol_Z_XZ) > ring_threshold:
                if trans <= best_loss:
                    best_loss = trans
                    trans_val = np.linspace(trans, 0.905, 10)
                    print(trans)
                    best_graph = gstate
                    print("XZ")
                    print(trans_val)
                else:
                    break
            else:
                break
    print("Best loss threshold value {}".format(best_loss))
    best_graph.image(input_qubits=[in_qubit])
    plt.show()
    '''
        summed_wt_list_new_XZ = add_X_Z_wt_list(tab_list_X_new, tab_list_Z_new)
        summed_wt_list_new_YZ = add_X_Z_wt_list(tab_list_Y_new, tab_list_Z_new)
        summed_wt_list_new_XY = add_X_Z_wt_list(tab_list_X_new, tab_list_Y_new)
        wt_dict_new_XY, meas_base_XY = extract_wt_from_tab_list_with_Y(summed_wt_list_new_XY)
        wt_dict_new_YZ, meas_base_YZ = extract_wt_from_tab_list_with_Y(summed_wt_list_new_YZ)
        wt_dict_new_XZ, meas_base_XZ = extract_wt_from_tab_list_with_Y(summed_wt_list_new_XZ)
        # print("adapt dict " + str(tab_list_X_new))
        # print("adapt dict " + str(tab_list_Z_new))
        # print(wt_dict_new_XZ)
        # print(meas_base_XZ)
        analytic_sol_X_XY = convert_analystic_solution_with_Y(adapt_dict_X_new, wt_dict_new, 1 / 2, trans, meas_base)
        analytic_sol_Y_XY = convert_analystic_solution_with_Y(adapt_dict_Y_new, wt_dict_new, 1 / 2, trans, meas_base)
        analytic_sol_Z_YZ = convert_analystic_solution_with_Y(adapt_dict_Z_new, wt_dict_new_YZ, 1 / 2, trans, meas_base_YZ)
        analytic_sol_Y_YZ = convert_analystic_solution_with_Y(adapt_dict_Y_new, wt_dict_new_YZ, 1 / 2, trans, meas_base_YZ)
        analytic_sol_Z_XZ = convert_analystic_solution_with_Y(adapt_dict_Z_new, wt_dict_new_XZ, 1 / 2, trans,
                                                           meas_base_XZ)
        analytic_sol_X_XZ = convert_analystic_solution_with_Y(adapt_dict_X_new, wt_dict_new_XZ, 1 / 2, trans,
                                                           meas_base_XZ)
        if min(analytic_sol_X_XZ, analytic_sol_Z_XZ) > 0.88:
            print("XZ")
            print(analytic_sol_X_XZ, analytic_sol_Z_XZ)
            gstate.image(input_qubits=[in_qubit])
            plt.show()
        elif min(analytic_sol_Y_YZ, analytic_sol_Z_YZ) > 0.88:
            print("YZ")
            print(analytic_sol_Y_YZ, analytic_sol_Z_YZ)
            gstate.image(input_qubits=[in_qubit])
            plt.show()
        elif min(analytic_sol_Y_XY, analytic_sol_X_XY) > 0.88:
            print("XY")
            print(analytic_sol_Y_XY, analytic_sol_X_XY)
            gstate.image(input_qubits=[in_qubit])
            plt.show()
        '''
