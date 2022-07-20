from LossToleranceFunctions.Losstolerant_fusion import*
from LossToleranceFunctions.Erasure_threshold_codes import*
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
from tqdm import tqdm




## NOTE 13 COULD BE 134....
def XZ(w, p_fail, transmission):
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
        adapt_coeff = adapt_coeff_in
        for item in adapt_coeff:
            if len(item) == i:
                # adapt_coeff.remove(item)
                for item_1 in adapt_coeff:
                    if len(item_1) > i and item[:i] == item_1[:i]:
                        # if item[:i] == item_1[:i]:
                        if item_1[i] == "3":
                            tabs_list[str(item_1[:i])][0] += 1
                        elif item_1[i] == "2":
                            tabs_list[str(item_1[:i])][1] += 1
    #print(tabs_list)
    #print(adapt_coeff_in)
    return tabs_list

def extract_wt_from_tab_list_with_Y(tab_list):
    wt_dict = {}
    for key in tab_list.keys():
        Z, X, Y = tab_list[key]
        one = 0
        two = 0
        if Z >= X:
            one = Z
            if Y >= X:
                two = Y
                meas_bases = "13"
            else:
                two = X
                meas_bases = "12"
        else:
            one = X
            if Y >= Z:
                two = Y
                meas_bases = "23"
            else:
                two = Z
                meas_bases = "21"

        if one > two:
            wt = 0
        elif one < two:
            wt = 1
        elif one == two:
            wt = 1 / 2
        wt_dict[key] = wt
    return wt_dict, meas_bases


def extract_wt_from_tab_list(tab_list):
    wt_dict = {}
    for key in tab_list.keys():
        Z,X = tab_list[key]
        if Z > X:
            wt = 0
        elif Z < X:
            wt = 1
        elif X == Z:
            wt = 1 / 2
        # wt = X / (X + Z)
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

def extract_first_weigth_erasure(XZ_1, XZ_2):
    X = 0
    Z = 0
    for key in XZ_1.keys():
        if key == "X":
            X = XZ_1[key] + XZ_2[key]
        else:
            Z = XZ_1[key] + XZ_2[key]
    if X > Z:
        wt = 1
    elif X < Z:
        wt = 0
    elif X == Z:
        wt = 1 / 2
    return wt


def extract_first_weigth_erasure_new(XZ_1, XZ_2):
    X = 0
    Z = 0
    for key in XZ_1.keys():
        if key == "X":
            X = XZ_1[key] + XZ_2[key]
        else:
            Z = XZ_1[key] + XZ_2[key]
    if X > Z:
        wt = 1
    elif X < Z:
        wt = 0
    elif X == Z:
        wt = 1 / 2
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


def add_X_Z_wt_list_new(tab_list_X, tab_list_Z):
    final_dict = {}
    for key in tab_list_X.keys():
        final_dict[key] = tab_list_X[key]
    for key in tab_list_Z.keys():
        if key in final_dict.keys():
            XX = final_dict[key][0]
            XZ = final_dict[key][1]
            ZX = tab_list_Z[key][0]
            ZZ = tab_list_Z[key][1]
            if min(XX, ZX) > min(XZ, ZZ):
                # wt = 1
                final_dict[key][0] = 1
                final_dict[key][1] = 0
            elif min(XX, ZX) < min(XZ, ZZ):
                # wt = 0
                final_dict[key][0] = 0
                final_dict[key][1] = 1
            elif min(XX, ZX) == min(XZ, ZZ):
                # wt = 1
                final_dict[key][0] = 1
                final_dict[key][1] = 0
            # final_dict[key][0] += tab_list_Z[key][0]
            # final_dict[key][1] += tab_list_Z[key][1]
        else:
            final_dict[key] = tab_list_Z[key]
    return final_dict


def adaptive_erasure_search(gstate):
    decoder_X = AdaptiveErasureDecoder(gstate, "X", numb_qubits)
    decoder_Z = AdaptiveErasureDecoder(gstate, "Z", numb_qubits)
    coef_X, tab_list_X_new, adapt_dict_X_new = decoder_X.run_dec_updated()
    coef_Z, tab_list_Z_new, adapt_dict_Z_new = decoder_Z.run_dec_updated()

    summed_wt_list_new = add_X_Z_wt_list_new(tab_list_X_new, tab_list_Z_new)
    # summed_wt_list_new.pop('', None)
    wt_dict = extract_wt_from_tab_list(summed_wt_list_new)
    # first_wt_X_new = extract_first_weigth_coeff_erasure(adapt_dict_X_new)
    # first_wt_Z_new = extract_first_weigth_coeff_erasure(adapt_dict_Z_new)
    # first_wt_new = extract_first_weigth_erasure(first_wt_X_new, first_wt_Z_new)
    first_wt = 0
    return adapt_dict_X_new, adapt_dict_Z_new, wt_dict, first_wt
    # return adapt_dict_X_new, adapt_dict_Z_new, wt_dict_new, first_wt_new


def calculate_adaptive_solution(adapt_dict_X, adapt_dict_Z, wt_dict, first_wt, p_fail, trans):
    analytic_sol_X = convert_analystic_solution(adapt_dict_X, wt_dict, first_wt, p_fail, trans)
    analytic_sol_Z = convert_analystic_solution(adapt_dict_Z, wt_dict, first_wt, p_fail, trans)
    return analytic_sol_X, analytic_sol_Z


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



def get_p_erase(p_fails, coeff_X, coeff_Z, weigths, transvals):
    loss_threshold = []
    beat = False
    for p_f in p_fails:
        for t in transvals:
            average_best = 0
            placeholder_weigth = 1
            wt_before = 0.001
            #for wt in weigths:
            for i in range(len(weigths)):
                X = convert_coeff(coeff_X, p_f, wt_before, t)
                Z = convert_coeff(coeff_Z, p_f, wt_before, t)
                average_before = (X + Z) / 2
                wt = weigths[i]
                wt_before = wt
                X = convert_coeff(coeff_X, p_f, wt, t)
                Z = convert_coeff(coeff_Z, p_f, wt, t)
                average = (X + Z) / 2
                # print(wt)
                if average < average_before:
                    break
                if average > average_best:
                    average_best = average
                    placeholder_weigth = wt
            if average_best > ring_threshold:
                print(average_best, placeholder_weigth, t)
                loss_threshold.append(1 - t)
                print(p_f)
                beat = True
                break
        if average_best < ring_threshold:
            loss_threshold.append(0)
    return loss_threshold, beat



def smooth_plt(p_fails, p_er, col, numb):
    model = make_interp_spline(p_fails, p_er)
    xs = np.linspace(0, 1 / 2, 500)
    ys = model(xs)
    plt.plot(xs, ys, '--', color=col, label="N=" + str(numb))


def smooth_shor(p_fails, p_er, col, numb):
    model = make_interp_spline(p_fails, p_er)
    xs = np.linspace(0, 1 / 2, 500)
    ys = model(xs)
    plt.plot(xs, ys, '--', color=col, label="shor-code")




if __name__ == '__main__':
    # TODO : Find Better way of deciding the weight
    # IN OLD COUNTING WE ONLY USED ONE OF THE LAST SEQUENCE FOR LOST OR SUCCESFULL FUSION
    #################### VARIABLES #######################
    in_qubit = 0
    p_fails = [1 / 2]  # [1 / 2, 1 / 4, 1 / 8, 1 / 16, 1 / 32]
    transvals = np.linspace(0.955, 1, 100)
    weigths = np.linspace(0.001, 0.999, 100)
    # weigths = [1/2]
    ring_threshold = 1 - 0.1198
    # ghz_threshold = 1 - 0.069
    numb_qubits = 6
    p_fail = 1 / 2
    trans = 0.945
    '''
    for i in range(1):
        print(i)
        #################### INITIALIZE #######################
        # edges = [(1, 2), (0, 3), (3, 2)]
        # graph = nx.Graph()
        # graph = gen_star_graph(numb_qubits)
        graph = gen_random_connected_graph(numb_qubits)
        gstate = GraphState(graph)
        decoder_X = AdaptiveErasureDecoder(gstate, "X", numb_qubits)
        decoder_Z = AdaptiveErasureDecoder(gstate, "Z", numb_qubits)
        coef_X, tab_list_X_new, adapt_dict_X_new = decoder_X.run_dec_updated()
        coef_Z, tab_list_Z_new, adapt_dict_Z_new = decoder_Z.run_dec_updated()
        print("decode done")
        
        adapt_dict_X, adpat_cof_X = adaptive_coeff(coef_X)
        # print(adapt_dict_X)
        adapt_dict_Z, adpat_cof_Z = adaptive_coeff(coef_Z)
        last_qubits_X = extract_last_meas(adapt_dict_X, numb_qubits - 1)
        adpat_cof_final_X = adpat_cof_X + last_qubits_X
        last_qubits_Z = extract_last_meas(adapt_dict_Z, numb_qubits - 1)
        adpat_cof_final_Z = adpat_cof_Z + last_qubits_Z
        tab_list_X = turn_adaptive_coeff_to_weigth(adpat_cof_final_X, numb_qubits)
        tab_list_Z = turn_adaptive_coeff_to_weigth(adpat_cof_final_Z, numb_qubits)

        # wt_dict_X = extract_wt_from_tab_list(tab_list_X)
        # wt_dict_Z = extract_wt_from_tab_list(tab_list_Z)
        # print(tab_list_Z)
        # print(tab_list_Z_new)
        # summed_wt_list = add_X_Z_wt_list(tab_list_X, tab_list_Z)
        summed_wt_list = add_X_Z_wt_list_new(tab_list_X_new, tab_list_Z_new)
        # print(summed_wt_list)
        # print(summed_wt_list)
        wt_dict = extract_wt_from_tab_list(summed_wt_list)
        first_wt_X = extract_first_weigth_coeff_erasure(adapt_dict_X)
        first_wt_Z = extract_first_weigth_coeff_erasure(adapt_dict_Z)
        first_wt = extract_first_weigth_erasure(first_wt_X, first_wt_Z)
        print(wt_dict)
        
        # wt_dict_X = extract_wt_from_tab_list(tab_list_X_new)
        # wt_dict_Z = extract_wt_from_tab_list(tab_list_Z_new)
        # summed_wt_list_new = add_X_Z_wt_list(tab_list_X_new, tab_list_Z_new)
        # print(summed_wt_list_new)
        # summed_wt_list_new.pop('', None)
        # print("Z list {}".format(tab_list_Z))
        # print("X list {}".format(tab_list_X))
        # print("new Z {}".format(tab_list_Z_new))
        # print("new X {}".format(tab_list_X_new))
        print("extracted weigths")
        # wt_dict_new = extract_wt_from_tab_list(summed_wt_list_new)
        # first_wt_X_new = extract_first_weigth_coeff_erasure(adapt_dict_X_new)
        # first_wt_Z_new = extract_first_weigth_coeff_erasure(adapt_dict_Z_new)
        # first_wt_new = extract_first_weigth_erasure(first_wt_X_new, first_wt_Z_new)

        #################### FIND SOLUTION #######################
        # analytic_sol_X = convert_analystic_solution(adapt_dict_X_new, wt_dict_new, first_wt_new, p_fail, trans)
        # analytic_sol_Z = convert_analystic_solution(adapt_dict_Z_new, wt_dict_new, first_wt_new, p_fail, trans)
        # print("NEW {} {}".format(analytic_sol_Z, analytic_sol_X))
        analytic_sol_X = convert_analystic_solution(adapt_dict_X, wt_dict, first_wt, p_fail, trans)
        analytic_sol_Z = convert_analystic_solution(adapt_dict_Z, wt_dict, first_wt, p_fail, trans)
        print("OLD {} {}".format(analytic_sol_Z, analytic_sol_X))
        if min(analytic_sol_Z, analytic_sol_X) > ring_threshold:
            print(analytic_sol_X, analytic_sol_Z)
            gstate.image(input_qubits=[in_qubit])
            plt.show()
    '''

    path = r'C:\Users\Admin\data_9_qubits.json'

    f = open(path)
    the_dict = json.load(f)
    loss_threshold_final = []
    best_val = 0
    loop_dict = the_dict["6"]
    trans_val = np.linspace(0.99, 0.98, 10)
    best_loss = 1
    print(len(loop_dict))
    i = 0
    for edges in loop_dict:
        print(i)
        i += 1
        graph = nx.Graph()
        graph.add_edges_from(edges)
        # graph = gen_random_connected_graph(numb_qubits)
        gstate = GraphState(graph)
        adapt_dict_X, adapt_dict_Z, wt_dict, first_wt = adaptive_erasure_search(gstate)
        for t in trans_val:
            X_sol, Z_sol = calculate_adaptive_solution(adapt_dict_X, adapt_dict_Z, wt_dict, first_wt
                                                       , p_fail, t)
            if (Z_sol + X_sol)/2 > ring_threshold:
            # if min(Z_sol, X_sol) > ring_threshold:
                if t < best_loss:
                    best_loss = t
                    trans_val = np.linspace(t, 0.95, 10)
                    print(t)
                    print(Z_sol, X_sol)
                    best_graph = gstate
            else:
                # print(Z_sol, X_sol)
                break
    print("best loss {}".format(best_loss))
    best_graph.image(input_qubits=[in_qubit])
    plt.show()

    '''
    for i in range(100):
        #################### INITIALIZE #######################
        graph = gen_random_connected_graph(numb_qubits)
        gstate = GraphState(graph)
        decoder_fusion = AdaptiveFusionDecoder(gstate)
        coeff_fusion = decoder_fusion.run_dec_updated()
        adapt_dict, adpat_cof = adaptive_coeff(coeff_fusion)
        last_qubits = extract_last_meas(adapt_dict, numb_qubits-1)
        adpat_cof_final = adpat_cof + last_qubits
        tab_list = turn_adaptive_coeff_to_weigth(adpat_cof_final, numb_qubits)
        wt_dict = extract_wt_from_tab_list(tab_list)
        first_wt = extract_first_weigth(adapt_dict)
        
        #################### FIND SOLUTION #######################
        analytic_sol = convert_analystic_solution(adapt_dict, wt_dict, first_wt, 1 / 4, 0.95)
        print(i)
        if analytic_sol > ring_threshold:
            print(analytic_sol)
            gstate.image(input_qubits=[in_qubit])
            plt.show()
    '''
    '''
    # paths_X = ["X_coeff_4_qubits_QD_p_fail_0.5.json", "X_coeff_5_qubits_QD_p_fail_0.5.json",
    #            "X_coeff_6_qubits_QD_p_fail_0.5.json","X_coeff_7_qubits_QD_p_fail_0.5.json",
    #            "X_coeff_9_qubits_QD_p_fail_0.5.json", "X_coeff_5_qubits_arb_p_fail_0.5.json", "X_coeff_5_qubits_QD_p_fail_0.5.json", "X_coeff_7_qubits.json",
    #            "X_coeff_7_qubits_arb_p_fail_0.5.json", "X_coeff_9_qubits_arb_p_fail_0.5.json"]
    
    paths_X = [ "X_coeff_5_qubits_arb_p_fail_0.5.json", "X_coeff_5_qubits_QD_p_fail_0.5.json",
                "X_coeff_7_qubits.json",
               "X_coeff_7_qubits_arb_p_fail_0.5.json", "X_coeff_9_qubits_arb_p_fail_0.5.json"]
    # paths_Z = ["Z_coeff_4_qubits_QD_p_fail_0.5.json", "Z_coeff_5_qubits_QD_p_fail_0.5.json",
    #            "Z_coeff_6_qubits_QD_p_fail_0.5.json","Z_coeff_7_qubits_QD_p_fail_0.5.json",
    #            "Z_coeff_9_qubits_QD_p_fail_0.5.json", "Z_coeff_5_qubits_arb_p_fail_0.5.json", "Z_coeff_5_qubits_QD_p_fail_0.5.json", "Z_coeff_7_qubits.json",
    #            "Z_coeff_7_qubits_arb_p_fail_0.5.json", "Z_coeff_9_qubits_arb_p_fail_0.5.json"
    #            ]
    
    paths_Z = ["Z_coeff_5_qubits_arb_p_fail_0.5.json", "Z_coeff_5_qubits_QD_p_fail_0.5.json",
               "Z_coeff_7_qubits.json",
               "Z_coeff_7_qubits_arb_p_fail_0.5.json", "Z_coeff_9_qubits_arb_p_fail_0.5.json"]
    loss_4 = 0
    loss_5 = 0
    loss_6 = 0
    loss_7 = 0
    loss_8 = 0
    loss_9 = 0
    loss_10 = 0
    loss_11 = 0
    loss_12 = 0
    for i in range(5):
        dict_x = open(paths_X[i])
        dict_z = open(paths_Z[i])
        dict_x = json.load(dict_x)
        dict_z = json.load(dict_z)
        print(paths_Z[i])
        loss_threshold, beat = get_p_erase(p_fails, dict_x, dict_z, weigths, transvals)
        val_plot = []
        for val in loss_threshold:
            if val == 0:
                val_plot.append(None)
            else:
                val_plot.append(val)
        if i == 0:
            loss_4 = val_plot[0]
        elif i == 1:
            loss_5 = val_plot[0]
        elif i == 2:
            loss_6 = val_plot[0]
        elif i == 3:
            loss_7 = val_plot[0]
        elif i == 4:
            loss_8 = val_plot[0]
        elif i == 5:
            loss_9 = val_plot[0]
        elif i == 6:
            loss_10 = val_plot[0]
        elif i == 7:
            loss_11 = val_plot[0]
        else:
            loss_12 = val_plot[0]
        #plt.plot(p_fails, val_plot, "-o", label="N="+str(i+5))
    
    # paths_X = ["X_coeff_4_qubits_QD_p_fail_0.25.json", "X_coeff_5_qubits_QD_p_fail_0.25.json",
    #            "X_coeff_6_qubits_QD_p_fail_0.25.json","X_coeff_7_qubits_QD_p_fail_0.25.json",
    #            "X_coeff_9_qubits_QD_p_fail_0.25.json", "X_coeff_5_qubits_2.json", "X_coeff_6_qubits.json", "X_coeff_7_qubits.json",
    #            "X_coeff_8_qubits.json", "X_coeff_9_qubits_arb_p_fail_0.25.json"
    #            ]
    paths_X = [ "X_coeff_5_qubits_2.json", "X_coeff_6_qubits.json", "X_coeff_7_qubits.json",
               "X_coeff_8_qubits.json", "X_coeff_9_qubits_arb_p_fail_0.25.json"]
    # paths_Z = ["Z_coeff_4_qubits_QD_p_fail_0.25.json", "Z_coeff_5_qubits_QD_p_fail_0.25.json",
    #            "Z_coeff_6_qubits_QD_p_fail_0.25.json","Z_coeff_7_qubits_QD_p_fail_0.25.json",
    #            "Z_coeff_9_qubits_QD_p_fail_0.25.json", "Z_coeff_5_qubits_2.json", "Z_coeff_6_qubits.json", "Z_coeff_7_qubits.json",
    #            "Z_coeff_8_qubits.json", "Z_coeff_9_qubits_arb_p_fail_0.25.json"
    #            ]
    
    paths_Z = ["Z_coeff_5_qubits_2.json", "Z_coeff_6_qubits.json", "Z_coeff_7_qubits.json",
               "Z_coeff_8_qubits.json", "Z_coeff_9_qubits_arb_p_fail_0.25.json"]
    
    saved_dict = {}
    for i in range(5):
        p_fails = [1 / 4, 1 / 8, 1 / 16, 1 / 32]
        dict_x = open(paths_X[i])
        dict_z = open(paths_Z[i])
        dict_x = json.load(dict_x)
        dict_z = json.load(dict_z)
        loss_threshold, beat = get_p_erase(p_fails, dict_x, dict_z, weigths, transvals)
        val_plot = []
        if i == 0:
            val_plot.append(loss_4)
        elif i == 1:
            val_plot.append(loss_5)
        elif i == 2:
            val_plot.append(loss_6)
        elif i == 3:
            val_plot.append(loss_7)
        elif i == 4:
            val_plot.append(loss_8)
        elif i == 5:
            val_plot.append(loss_9)
        elif i == 6:
            val_plot.append(loss_10)
        elif i == 7:
            val_plot.append(loss_11)
        else:
            val_plot.append(loss_12)
        for val in loss_threshold:
            if val == 0:
                val_plot.append(None)
            else:
                val_plot.append(val)
        p_fails = [1 / 2, 1 / 8, 1 / 16, 1 / 32]
        p_fails_inv = np.array([1/32, 1/16, 1/8, 1/4,1/2])
        val_plot_inv = val_plot[::-1]
        saved_dict[str(i)]=val_plot
        print(val_plot_inv)
        val_25 = val_plot[1]
        del val_plot[1]
        val_plot_inv = np.array(val_plot_inv)
        print(i)
        if i <= 4:
            if i == 0:
                plt.plot(p_fails, val_plot, "o", color="red", markersize=7)
                plt.plot([1/4], val_25, "*", color="red", markersize=10)
                smooth_plt(p_fails_inv, val_plot_inv, "red", i+5)
                # plt.plot(p_fails, val_plot, color="red", label="N=" + str(i + 5))
            elif i == 1:
                plt.plot(p_fails, val_plot, "o", color="blue", markersize=7)
                plt.plot([1 / 4], val_25, "*", color="blue", markersize=10)
                smooth_plt(p_fails_inv, val_plot_inv, "blue", i + 5)
                # plt.plot(p_fails, val_plot, color="blue", label="N=" + str(i + 5))
            elif i == 2:
                plt.plot(p_fails, val_plot, "o", color="green", markersize=7)
                plt.plot([1 / 4], val_25, "*", color="green", markersize=10)
                smooth_plt(p_fails_inv, val_plot_inv, "green", i + 5)
                # plt.plot(p_fails, val_plot, color="green", label="N=" + str(i + 5))
            elif i == 3:
                plt.plot(p_fails, val_plot, "o", color="black", markersize=7)
                plt.plot([1 / 4], val_25, "*", color="black", markersize=10)
                smooth_plt(p_fails_inv, val_plot_inv, "black", i + 5)
                # plt.plot(p_fails, val_plot, color="black", label="N=" + str(i + 5))
            else:
                plt.plot(p_fails, val_plot, "o", color="orange", markersize=7)
                plt.plot([1 / 4], val_25, "*", color="orange", markersize=10)
                smooth_plt(p_fails_inv, val_plot_inv, "orange", i + 5)
                # plt.plot(p_fails, val_plot, color="orange", label="N=" + str(i + 5))
        else:
            if i == 4:
                plt.plot(p_fails, val_plot, "-o", color="red")
            elif i == 5:
                plt.plot(p_fails, val_plot, "-o", color="blue")
            elif i == 6:
                plt.plot(p_fails, val_plot, "-o", color="green")
            elif i == 7:
                plt.plot(p_fails, val_plot, "-o", color="black")
            else:
                plt.plot(p_fails, val_plot, "-o", color="orange")
            # plt.plot(p_fails, val_plot, "-o", label="N=" + str(i + 5))
    with open('GHZ-network.json', 'w') as outfile:
        json.dump(saved_dict, outfile)
    # shors = [0.027070707070707134, 0.021818181818181848, 0.012929292929292902, 0.006868686868686913]
    shors = np.array([0.006868686868686913, 0.012929292929292902,0.021818181818181848,0.027070707070707134])
    p_fails = np.array([1/32, 1/16, 1/8, 1/4])
    smooth_shor(p_fails, shors, "purple", "shor-code")
    p_fails = np.array([1/32, 1/16, 1/8])
    shors = np.array([0.006868686868686913, 0.012929292929292902, 0.021818181818181848])
    plt.plot(p_fails, shors, "o", color="purple", markersize=7)
    p_fails = np.array([1/4])
    shors = np.array([0.027070707070707134])
    plt.plot(p_fails, shors, "*", color="purple", markersize=10)
    
    # plt.plot(p_fails, shors, "-*", color="purple", label="shor-code")
    plt.ylim(0, None)
    plt.xlim(0, None)
    plt.title("Ring-network")
    plt.ylabel("Photon loss threshold ($p_{loss}$)")
    plt.xlabel("Fusion failure probability ($p_{fail}$)")
    plt.legend()
    plt.show()
    
    '''
    path = r'C:\Users\Admin\data_9_qubits.json'
    path_save = r'C:\Users\Admin\Desktop\project_outside_course\erasure_beating_shors'
    '''
    f = open(path)
    the_dict = json.load(f)
    loss_threshold_final = []
    best_val = 0
    # graph_list = [i for i in range(300)]
    # for i in tqdm(graph_list):
    p_fail = [[1 / 2]]
    keys = ["4", "5", "6", "7", "8"]
    loop_dict = the_dict["4"]
    print(len(loop_dict))
    for p_fails in p_fail:
        current_loss_threshold = [1]
        i = 0
        for edges in loop_dict:
            print(i)
            i += 1
            if i > 1:
                break
            graph = nx.Graph()
            graph.add_edges_from(edges)
            # graph = gen_random_connected_graph(9)
            gstate = GraphState(graph)
            decoder_Z = ErasureDecoder(gstate, "Z", 1 / 2, 1, 0, in_qubit)
            succ_Z, coeff_Z = decoder_Z.run_dec_updated()
            decoder_X = ErasureDecoder(gstate, "X", 1 / 2, 1, 0, in_qubit)
            succ_X, coeff_X = decoder_X.run_dec_updated()
            loss_threshold, beat = get_p_erase(p_fails, coeff_X, coeff_Z, weigths, current_loss_threshold)
            print("first beat {}".format(beat))
            print("current threshold {}".format(current_loss_threshold))
            if beat == True:
                loss_threshold, beat = get_p_erase(p_fails, coeff_X, coeff_Z, weigths, transvals)
            if (1-loss_threshold[0]) < current_loss_threshold[0]:
                current_loss_threshold = [(1-loss_threshold[0])]
                transvals = np.linspace(0.94, current_loss_threshold[0], 50)
            if max(loss_threshold) > best_val:
                best_val = max(loss_threshold)
                loss_threshold_final = (loss_threshold)
                best_graph = gstate
                best_X_coeff = coeff_X
                best_Z_coeff = coeff_Z
                print("We here")
                print(loss_threshold_final)

        # fig1 = plt.figure()
        # best_graph.image(input_qubits=[in_qubit])
        # plt.savefig(path_save + "\graph-QD-9-qubits-" + str(p_fails[0]) + '.png')
        # plt.close(fig1)
    
        # with open('X_coeff_9_qubits_QD_p_fail_' + str(p_fails[0]) + '.json', 'w') as outfile:
        #     json.dump(best_X_coeff, outfile)
        # with open('Z_coeff_9_qubits_QD_p_fail_' + str(p_fails[0]) + '.json', 'w') as outfile:
        #     json.dump(best_Z_coeff, outfile)
    '''