from LossToleranceFunctions.Losstolerant_fusion import*
import numpy as np
import qecc as q

#TODO: CALCULATE THE COEFFICIENTS AND THEN REUSE FOR DIFFERENT TRANSMISSION
class ErasureDecoder(object):
    def __init__(self, gstate, indmeas_pauli, p_fail, transmission, w, in_qubit=0):
        self.gstate = gstate
        self.indmeas_pauli = indmeas_pauli
        self.in_qubit = in_qubit
        self.p_fail = p_fail
        self.transmission = transmission
        self.w = w
        self.all_stabs = self.remove_paulis()
        self.configurations = self.get_configuration()
        self.succ_config = self.get_succ_conf()
        self.pattern = self.gen_meas_pattern()
        self.p_x = self.er_x()
        self.p_z = self.er_z()
        self.pattern_update = self.gen_meas_pattern_new()
        self.p_x_succ = self.p_fail*self.w*(self.transmission**(1/self.p_fail))
        self.p_z_succ = self.p_fail*(1-self.w)*(self.transmission**(1/self.p_fail))
        self.p_fusion_succ = (1 - p_fail)*(self.transmission**(1/self.p_fail))
        self.p_fusion_lost = (1 - (self.transmission**(1/self.p_fail)))

    def get_indirect_meas_stabs(self):
        num_qubits = len(self.gstate)
        # get all possible 2^N stabilizers.
        # TODO: to be improved checking only "smart" stabilizers.
        all_stabs = q.from_generators(self.gstate.stab_gens)
        poss_stabs_list = []
        for this_stab0 in all_stabs:
            this_stab = this_stab0.op
            if this_stab[self.in_qubit] == self.indmeas_pauli:
                meas_weight = num_qubits - this_stab.count('I')
                Z_weight = this_stab.count('Z') / (meas_weight + 1)
                meas_qubits = [qbt for qbt in range(num_qubits)
                               if ((qbt != self.in_qubit) and (this_stab[qbt] != 'I'))]
                ### order them such that we always prefer mstrategies with smaller weight, and with more Zs in the non-trivial paulis.
                poss_stabs_list.append([this_stab, meas_qubits, meas_weight, Z_weight])
        poss_stabs_list.sort(key=lambda x: x[2] - x[3])
        return poss_stabs_list

    def remove_paulis(self):
        new_stabs = []
        for stab in self.get_indirect_meas_stabs():
            Y_count = 0
            for idx in stab[0]:
                if idx == "Y":
                    Y_count += 1
            if Y_count == 0:
                new_stabs.append(stab)
        return new_stabs

    def remove_qubits(self, indx, stabs, pauli):
        # Pauli op is what we keep
        new_stabs = []
        for stab in stabs:
            if stab[0][indx] == pauli or stab[0][indx] == "I":
                new_stabs.append(stab)
        return new_stabs

    def get_configuration(self):
        numb = len(self.all_stabs[0][0])-1
        configurations = []
        #iter = product('+-', repeat=numb)
        iter = product("abcd", repeat=numb)
        for ch in iter:
            configurations.append(ch)
        return configurations

    def get_succ_conf(self):
        config = []
        for stab in self.all_stabs:
            config.append(stab[1])
        return config

    # def get_succ_conf_fixed(self):
    #     config = []
    #     for stab in self.all_stabs:
    #         operators = stab[0]
    #         tracker = 0
    #         for idnx in stab[1]:
    #             if operators[idnx] == self.pattern[idnx] or operators[idnx] == "I":
    #                 tracker += 0
    #             else:
    #                 tracker += 1
    #         if tracker == 0:
    #             config.append(stab[1])
    #     return config

    def gen_meas_pattern(self):
        pattern = {}
        for stabs in self.all_stabs:
            for qbits in stabs[1]:
                if qbits in pattern.keys():
                    continue
                else:
                    pattern[qbits] = (stabs[0][qbits])
        return pattern

    def gen_meas_pattern_new(self):
        patterns = []
        for stabs in self.all_stabs:
            pattern = {}
            for qbits in stabs[1]:
                pattern[qbits] = (stabs[0][qbits])
            patterns.append(pattern)
        return patterns

    def er_z(self):
        return (1 - (1 - self.w * self.p_fail) * (self.transmission ** (1 / self.p_fail)))

    def er_x(self):
        return (1 - (1 - (1 - self.w) * self.p_fail) * (self.transmission ** (1 / self.p_fail)))


    def run_dec_updated(self):
        coeff_list = []
        tot_term = 0
        for config in self.configurations:
            term = 1
            current_run = self.all_stabs
            current_dict = {"P_s":0, "P_X":0, "P_Z":0, "P_l":0}
            for i in range(len(config)):
                indx = i + 1
                if config[i] == "a":
                    term = term * self.p_fusion_succ
                    current_dict["P_s"] += 1
                elif config[i] == "b":
                    # REMOVE Z
                    term = term * self.p_x_succ
                    current_run = self.remove_qubits(indx, current_run, "X")
                    current_dict["P_X"] += 1
                    #current_run = self.remove_qubits(indx, current_run, "I")
                elif config[i] == "c":
                    # REMOVE X
                    term = term * self.p_z_succ
                    current_run = self.remove_qubits(indx, current_run, "Z")
                    current_dict["P_Z"] += 1

                    #current_run = self.remove_qubits(indx, current_run, "I")
                else:
                    # REMOVE BOTH Z AND X
                    term = term * self.p_fusion_lost
                    current_run = self.remove_qubits(indx, current_run, "I")
                    current_dict["P_l"] += 1
            if len(current_run) > 0:
                #print(term)
                #print(config)
                tot_term += term
                coeff_list.append(current_dict)
            else:
                tot_term = tot_term
        return tot_term, coeff_list

    def run_dec(self):
        tot_term = 0
        for config in self.configurations:
            succ_meas = []
            failed_meas = []
            for i in range(len(config)):
                if config[i] == "+":
                    succ_meas.append(i + 1)
                else:
                    failed_meas.append(i + 1)
            for set_in_list in self.succ_config:
                set_in_list = set(set_in_list)
                succ_meas = set(succ_meas)
                if set_in_list.issubset(succ_meas):
                    term = 1
                    for qbit in succ_meas:
                        if self.pattern[qbit] == "X":
                            term = term * (1 - self.p_x)
                        else:
                            term = term * (1 - self.p_z)

                    for qbit in failed_meas:

                        if self.pattern[qbit] == "X":
                            term = term * (self.p_x)
                        else:
                            term = term * (self.p_z)
                    tot_term += term
                    break
        return tot_term

class AdaptiveErasureDecoder(object):
    def __init__(self, gstate, indmeas_pauli, numb_qubits,in_qubit=0):
        self.gstate = gstate
        self.numb_qubits = numb_qubits
        self.indmeas_pauli = indmeas_pauli
        self.in_qubit = in_qubit
        self.all_stabs = self.remove_paulis()
        self.configurations = self.get_configuration()
        self.succ_config = self.get_succ_conf()
        self.pattern = self.gen_meas_pattern()
        self.pattern_update = self.gen_meas_pattern_new()

    def get_indirect_meas_stabs(self):
        num_qubits = len(self.gstate)
        # get all possible 2^N stabilizers.
        # TODO: to be improved checking only "smart" stabilizers.
        all_stabs = q.from_generators(self.gstate.stab_gens)
        poss_stabs_list = []
        for this_stab0 in all_stabs:
            this_stab = this_stab0.op
            if this_stab[self.in_qubit] == self.indmeas_pauli:
                meas_weight = num_qubits - this_stab.count('I')
                Z_weight = this_stab.count('Z') / (meas_weight + 1)
                meas_qubits = [qbt for qbt in range(num_qubits)
                               if ((qbt != self.in_qubit) and (this_stab[qbt] != 'I'))]
                ### order them such that we always prefer mstrategies with smaller weight, and with more Zs in the non-trivial paulis.
                poss_stabs_list.append([this_stab, meas_qubits, meas_weight, Z_weight])
        poss_stabs_list.sort(key=lambda x: x[2] - x[3])
        return poss_stabs_list

    def remove_paulis(self):
        new_stabs = []
        for stab in self.get_indirect_meas_stabs():
            Y_count = 0
            for idx in stab[0]:
                if idx == "Y":
                    Y_count += 1
            if Y_count == 0:
                new_stabs.append(stab)
        return new_stabs

    def remove_qubits(self, indx, stabs, pauli):
        # Pauli op is what we keep
        new_stabs = []
        for stab in stabs:
            if stab[0][indx] == pauli or stab[0][indx] == "I":
                new_stabs.append(stab)
        return new_stabs

    def get_configuration(self):
        numb = len(self.all_stabs[0][0])-1
        configurations = []
        #iter = product('+-', repeat=numb)
        iter = product("abcd", repeat=numb)
        for ch in iter:
            configurations.append(ch)
        return configurations

    def get_succ_conf(self):
        config = []
        for stab in self.all_stabs:
            config.append(stab[1])
        return config


    def gen_meas_pattern(self):
        pattern = {}
        for stabs in self.all_stabs:
            for qbits in stabs[1]:
                if qbits in pattern.keys():
                    continue
                else:
                    pattern[qbits] = (stabs[0][qbits])
        return pattern

    def gen_meas_pattern_new(self):
        patterns = []
        for stabs in self.all_stabs:
            pattern = {}
            for qbits in stabs[1]:
                pattern[qbits] = (stabs[0][qbits])
            patterns.append(pattern)
        return patterns

    def turn_adaptive_coeff_to_weigth(self, adapt_coeff_in, numb_qubits=4):
        tabs_list = {}
        for item in adapt_coeff_in:
            for i in range(1, numb_qubits):
                if len(item) > i:
                    tabs_list[item[:i]] = [0, 0]
        for i in range(1, numb_qubits):
            adapt_coeff = adapt_coeff_in
            # start = time.time()
            for item in adapt_coeff:
                if len(item) == i:
                    for item_1 in adapt_coeff:
                        if len(item_1) > i:
                            if item[:i] == item_1[:i]:
                                if item_1[i] == "3":
                                    tabs_list[str(item_1[:i])][0] += 1
                                elif item_1[i] == "2":
                                    tabs_list[str(item_1[:i])][1] += 1
        return tabs_list

    def extract_last_meas(self, coeff, numb_qubit):
        XZ_meas = []
        for pattern in coeff:
            initial = pattern["qubit_" + str(numb_qubit)][:3]
            if initial == "P_X":
                XZ_meas.append(pattern["qubit_" + str(numb_qubit)][4:] + str(2))
            elif initial == "P_Z":
                XZ_meas.append(pattern["qubit_" + str(numb_qubit)][4:] + str(3))
        return XZ_meas

    def run_dec_updated(self):
        coeff_list = []
        weight_dict = {}
        new_coeff_dict_full = []
        for config in self.configurations:
            qbit = []
            current_run = self.all_stabs
            current_dict = {}
            for i in range(len(config)):
                seq_numb = ''.join(str(x) for x in qbit)
                if seq_numb not in weight_dict.keys():
                    weight_dict[seq_numb] = [0, 0]
                indx = i + 1
                if config[i] == "a":
                    current_dict["qubit_"+str(indx)] = "P_S"
                elif config[i] == "b":
                    # REMOVE Z
                    current_run = self.remove_qubits(indx, current_run, "X")
                    current_dict["qubit_" + str(indx)] = "P_X"
                elif config[i] == "c":
                    # REMOVE X
                    current_run = self.remove_qubits(indx, current_run, "Z")
                    current_dict["qubit_" + str(indx)] = "P_Z"
                    #current_run = self.remove_qubits(indx, current_run, "I")
                else:
                    # REMOVE BOTH Z AND X
                    current_run = self.remove_qubits(indx, current_run, "I")
                    current_dict["qubit_" + str(indx)] = "P_L"
            if len(current_run) > 0:
                new_coeff_dict = {}
                coeff_list.append(current_dict)
                qbit = []
                for i in range(len(config)):
                    seq_numb = ''.join(str(x) for x in qbit)
                    if seq_numb not in weight_dict.keys():
                        weight_dict[seq_numb] = [0, 0]
                    indx = i + 1
                    if config[i] == "a":
                        seq_numb = ''.join(str(x) for x in qbit)
                        new_coeff_dict["qubit_" + str(indx)] = "P_S_" + seq_numb
                    elif config[i] == "b":
                        current_run = self.remove_qubits(indx, current_run, "X")
                        seq_numb = ''.join(str(x) for x in qbit)
                        new_coeff_dict["qubit_" + str(indx)] = "P_X_" + seq_numb
                        weight_dict[seq_numb][1] += 1
                    elif config[i] == "c":
                        current_run = self.remove_qubits(indx, current_run, "Z")
                        seq_numb = ''.join(str(x) for x in qbit)
                        new_coeff_dict["qubit_" + str(indx)] = "P_Z_" + seq_numb
                        weight_dict[seq_numb][0] += 1
                    else:
                        current_run = self.remove_qubits(indx, current_run, "I")
                        seq_numb = ''.join(str(x) for x in qbit)
                        new_coeff_dict["qubit_" + str(indx)] = "P_L_" + seq_numb
                    if config[i] == "a":
                        qbit.append(1)
                    elif config[i] == "b":
                        qbit.append(2)
                    elif config[i] == "c":
                        qbit.append(3)
                    elif config[i] == "d":
                        qbit.append(4)
                new_coeff_dict_full.append(new_coeff_dict)
        return coeff_list, weight_dict, new_coeff_dict_full



class AdaptiveErasureDecoderNew(object):
    def __init__(self, gstate, indmeas_pauli, numb_qubits,in_qubit=0):
        self.gstate = gstate
        self.numb_qubits = numb_qubits
        self.indmeas_pauli = indmeas_pauli
        self.in_qubit = in_qubit
        self.all_stabs = self.get_indirect_meas_stabs()
        self.configurations = self.get_configuration()
        self.succ_config = self.get_succ_conf()
        self.pattern = self.gen_meas_pattern()
        self.pattern_update = self.gen_meas_pattern_new()

    def get_indirect_meas_stabs(self):
        num_qubits = len(self.gstate)
        # get all possible 2^N stabilizers.
        # TODO: to be improved checking only "smart" stabilizers.
        all_stabs = q.from_generators(self.gstate.stab_gens)
        poss_stabs_list = []
        for this_stab0 in all_stabs:
            this_stab = this_stab0.op
            if this_stab[self.in_qubit] == self.indmeas_pauli:
                meas_weight = num_qubits - this_stab.count('I')
                Z_weight = this_stab.count('Z') / (meas_weight + 1)
                meas_qubits = [qbt for qbt in range(num_qubits)
                               if ((qbt != self.in_qubit) and (this_stab[qbt] != 'I'))]
                ### order them such that we always prefer mstrategies with smaller weight, and with more Zs in the non-trivial paulis.
                poss_stabs_list.append([this_stab, meas_qubits, meas_weight, Z_weight])
        poss_stabs_list.sort(key=lambda x: x[2] - x[3])
        return poss_stabs_list



    def remove_qubits(self, indx, stabs, pauli):
        # Pauli op is what we keep
        new_stabs = []
        for stab in stabs:
            if stab[0][indx] == pauli or stab[0][indx] == "I":
                new_stabs.append(stab)
        return new_stabs

    def get_configuration(self):
        numb = len(self.all_stabs[0][0])-1
        configurations = []
        #iter = product('+-', repeat=numb)
        iter = product("abcgd", repeat=numb)
        for ch in iter:
            configurations.append(ch)
        return configurations

    def get_succ_conf(self):
        config = []
        for stab in self.all_stabs:
            config.append(stab[1])
        return config


    def gen_meas_pattern(self):
        pattern = {}
        for stabs in self.all_stabs:
            for qbits in stabs[1]:
                if qbits in pattern.keys():
                    continue
                else:
                    pattern[qbits] = (stabs[0][qbits])
        return pattern

    def gen_meas_pattern_new(self):
        patterns = []
        for stabs in self.all_stabs:
            pattern = {}
            for qbits in stabs[1]:
                pattern[qbits] = (stabs[0][qbits])
            patterns.append(pattern)
        return patterns


    def run_dec_updated(self):
        coeff_list = []
        weight_dict = {}
        new_coeff_dict_full = []
        for config in self.configurations:
            qbit = []
            current_run = self.all_stabs
            current_dict = {}
            for i in range(len(config)):
                seq_numb = ''.join(str(x) for x in qbit)
                if seq_numb not in weight_dict.keys():
                    weight_dict[seq_numb] = [0, 0, 0]
                indx = i + 1
                if config[i] == "a":
                    current_dict["qubit_"+str(indx)] = "P_S"
                elif config[i] == "b":
                    # REMOVE Z
                    current_run = self.remove_qubits(indx, current_run, "X")
                    current_dict["qubit_" + str(indx)] = "P_X"
                elif config[i] == "c":
                    # REMOVE X
                    current_run = self.remove_qubits(indx, current_run, "Z")
                    current_dict["qubit_" + str(indx)] = "P_Z"
                    #current_run = self.remove_qubits(indx, current_run, "I")
                elif config[i] == "g":
                    current_run = self.remove_qubits(indx, current_run, "Y")
                    current_dict["qubit_" + str(indx)] = "P_Y"
                else:
                    # REMOVE BOTH Z AND X
                    current_run = self.remove_qubits(indx, current_run, "I")
                    current_dict["qubit_" + str(indx)] = "P_L"
            if len(current_run) > 0:
                new_coeff_dict = {}
                coeff_list.append(current_dict)
                qbit = []
                for i in range(len(config)):
                    seq_numb = ''.join(str(x) for x in qbit)
                    if seq_numb not in weight_dict.keys():
                        weight_dict[seq_numb] = [0, 0, 0]
                    indx = i + 1
                    if config[i] == "a":
                        seq_numb = ''.join(str(x) for x in qbit)
                        new_coeff_dict["qubit_" + str(indx)] = "P_S_" + seq_numb
                    elif config[i] == "b":
                        # current_run = self.remove_qubits(indx, current_run, "X")
                        seq_numb = ''.join(str(x) for x in qbit)
                        new_coeff_dict["qubit_" + str(indx)] = "P_X_" + seq_numb
                        weight_dict[seq_numb][1] += 1
                    elif config[i] == "c":
                        # current_run = self.remove_qubits(indx, current_run, "Z")
                        seq_numb = ''.join(str(x) for x in qbit)
                        new_coeff_dict["qubit_" + str(indx)] = "P_Z_" + seq_numb
                        weight_dict[seq_numb][0] += 1
                    elif config[i] == "g":
                        # current_run = self.remove_qubits(indx, current_run, "Z")
                        seq_numb = ''.join(str(x) for x in qbit)
                        new_coeff_dict["qubit_" + str(indx)] = "P_Y_" + seq_numb
                        weight_dict[seq_numb][2] += 1
                    else:
                        # current_run = self.remove_qubits(indx, current_run, "I")
                        seq_numb = ''.join(str(x) for x in qbit)
                        new_coeff_dict["qubit_" + str(indx)] = "P_L_" + seq_numb
                    if config[i] == "a":
                        qbit.append(1)
                    elif config[i] == "b":
                        qbit.append(2)
                    elif config[i] == "c":
                        qbit.append(3)
                    elif config[i] == "g":
                        qbit.append(5)
                    elif config[i] == "d":
                        qbit.append(4)
                new_coeff_dict_full.append(new_coeff_dict)
        return coeff_list, weight_dict, new_coeff_dict_full


class AdaptiveFusionDecoder(object):
    def __init__(self, gstate, in_qubit=0):
        self.gstate = gstate
        self.in_qubit = in_qubit
        self.all_stabs_X = self.remove_paulis("X")
        self.all_stabs_Z = self.remove_paulis("Z")
        self.configurations = self.get_configuration()

    def get_indirect_meas_stabs(self, meas_qubit):
        num_qubits = len(self.gstate)
        # get all possible 2^N stabilizers.
        # TODO: to be improved checking only "smart" stabilizers.
        all_stabs = q.from_generators(self.gstate.stab_gens)
        poss_stabs_list = []
        for this_stab0 in all_stabs:
            this_stab = this_stab0.op
            if this_stab[self.in_qubit] == meas_qubit:
                meas_weight = num_qubits - this_stab.count('I')
                Z_weight = this_stab.count('Z') / (meas_weight + 1)
                meas_qubits = [qbt for qbt in range(num_qubits)
                               if ((qbt != self.in_qubit) and (this_stab[qbt] != 'I'))]
                ### order them such that we always prefer mstrategies with smaller weight, and with more Zs in the non-trivial paulis.
                poss_stabs_list.append([this_stab, meas_qubits, meas_weight, Z_weight])
        poss_stabs_list.sort(key=lambda x: x[2] - x[3])
        return poss_stabs_list

    def remove_paulis(self, meas_qubit):
        new_stabs = []
        for stab in self.get_indirect_meas_stabs(meas_qubit):
            Y_count = 0
            for idx in stab[0]:
                if idx == "Y":
                    Y_count += 1
            if Y_count == 0:
                new_stabs.append(stab)
        return new_stabs

    def remove_qubits(self, indx, stabs, pauli):
        # Pauli op is what we keep
        new_stabs = []
        for stab in stabs:
            if stab[0][indx] == pauli or stab[0][indx] == "I":
                new_stabs.append(stab)
        return new_stabs

    def get_configuration(self):
        numb = len(self.all_stabs_X[0][0])-1
        configurations = []
        #iter = product('+-', repeat=numb)
        iter = product("abcd", repeat=numb)
        for ch in iter:
            configurations.append(ch)
        return configurations


    def run_dec_updated(self):
        coeff_list = []
        tot_term = 0
        for config in self.configurations:
            term = 1
            current_run_X = self.all_stabs_X
            current_run_Z = self.all_stabs_Z
            current_dict = {}
            for i in range(len(config)):
                indx = i + 1
                if config[i] == "a":
                    current_dict["qubit_"+str(indx)] = "P_S"
                elif config[i] == "b":
                    # REMOVE Z
                    current_run_X = self.remove_qubits(indx, current_run_X, "X")
                    current_run_Z = self.remove_qubits(indx, current_run_Z, "X")
                    current_dict["qubit_" + str(indx)] = "P_X"
                elif config[i] == "c":
                    # REMOVE X
                    current_run_X = self.remove_qubits(indx, current_run_X, "Z")
                    current_run_Z = self.remove_qubits(indx, current_run_Z, "Z")
                    current_dict["qubit_" + str(indx)] = "P_Z"

                    #current_run = self.remove_qubits(indx, current_run, "I")
                else:
                    # REMOVE BOTH Z AND X
                    current_run_X = self.remove_qubits(indx, current_run_X, "I")
                    current_run_Z = self.remove_qubits(indx, current_run_Z, "I")
                    current_dict["qubit_" + str(indx)] = "P_L"
            if len(current_run_X) > 0 and len(current_run_Z) > 0:
                tot_term += term
                coeff_list.append(current_dict)
            else:
                tot_term = tot_term
        return coeff_list


'''
def erasure_z(p_fail, t, w):
    return (1 - (1-w*p_fail)*(t**(1/p_fail)))


def erasure_x(p_fail, t, w):
    return (1 - (1-(1-w)*p_fail)*(t**(1/p_fail)))

def gen_erase_prob_variable(p_z, p_x, configurations, succ_conf, pattern):
    tot_term = 0
    for config in configurations:
        succ_meas = []
        failed_meas = []
        for i in range(len(config)):
            if config[i] == "+":
                succ_meas.append(i+1)
            else:
                failed_meas.append(i+1)
        for set_in_list in succ_conf:
            set_in_list = set(set_in_list)
            succ_meas = set(succ_meas)
            if set_in_list.issubset(succ_meas):
                term = 1
                for qbit in succ_meas:
                    if pattern[qbit] == "X":
                        term = term * (1-p_x)
                    else:
                        term = term* (1- p_z)
                for qbit in failed_meas:
                    if pattern[qbit] == "X":
                        term = term * (p_x)
                    else:
                        term = term* (p_z)
                tot_term += term
                break
    return tot_term


def run_new_erase_decoder_zx(gstate, p_fail, loss, s):
    values = []
    decod0 = LT_IndMeasDecoder(gstate, 'X', in_qubit)
    stab_X = decod0.all_stabs
    stab_X = remove_Y(stab_X)
    decod1 = LT_IndMeasDecoder(gstate, 'Z', in_qubit)
    stab_Z = decod1.all_stabs
    stab_Z = remove_Y(stab_Z)
    configurations = get_configuration(len(stab_X[0][0])-1)
    succ_config_Z = get_succ_conf(stab_Z)
    succ_config_X = get_succ_conf(stab_X)
    pattern_Z = gen_meas_pattern(stab_Z)
    pattern_X  = gen_meas_pattern(stab_X)
    for t in loss:
        p_x = erasure_x(p_fail, t, s)
        p_z = erasure_z(p_fail, t, s)
        succ_X = gen_erase_prob_variable(p_z, p_x, configurations, succ_config_X, pattern_X)
        succ_Z = gen_erase_prob_variable(p_z, p_x, configurations, succ_config_Z, pattern_Z)
        fail_X = 1 - succ_X
        fail_Z = 1 - succ_Z
        print(succ_X)
        average = (fail_X + fail_Z) / 2
        values.append(average)
    return values
'''
if __name__ == '__main__':
    in_qubit = 0
    # graph = gen_star_graph(5)
    # graph = nx.Graph()
    # graph_edges = [(1, 2), (2, 0), (0, 3), (3, 4)]
    # graph.add_edges_from(graph_edges)
    # gstate = GraphState(graph)
    path = r'C:\Users\Admin\data_9_qubits.json'
    f = open(path)
    the_dict = json.load(f)
    p_fail = 1 / 2
    transmission = 1
    weigth = np.linspace(0.001, 0.99, 100)
    transmissions = [0.98]#np.linspace(0.97, 1, 30)
    average_best = 0.15
    path_save = r'C:\Users\Admin\Desktop\project_outside_course\graph_beating_GHZ\QD-graphs'
    graph = gen_star_graph(4)
    # graph = nx.Graph()
    # graph_edges = [(1, 2), (2, 3), (3, 4), (3, 0), (0, 2)]
    # graph.add_edges_from(graph_edges)
    gstate = GraphState(graph)
    # decoder_Z = AdaptiveErasureDecoder(gstate, "Z", in_qubit)
    # z_succ = (decoder_Z.run_dec_updated())
    # print(decoder_Z.all_stabs)
    # decoder_X = ErasureDecoder(gstate, "X", 1 / 2, 0.96, 0.0001, in_qubit)
    # print(decoder_X.all_stabs)
    # x_succ = (decoder_X.run_dec_updated())
    # print((x_succ + z_succ) / 2)
    '''
    keys = ["4", "5", "6", "7", "8"]
    i = 0
    print(len(the_dict["5"]))
    # graph = gen_star_graph(4)
    # gstateGHZ = GraphState(graph)
    for key in keys:
        loop_dict = the_dict[key]
        graph = gen_star_graph(int(key)+1)
        gstateGHZ = GraphState(graph)
        for edges in loop_dict:
            i += 1
    # for i in range(5, 200):
        #print(i)
        # graph = gen_random_connected_graph(4)
            graph = nx.Graph()
            graph.add_edges_from(edges)
            # graph = nx.Graph()
            # graph_edges = [(1, 2), (2, 3), (3, 4), (3, 0), (0, 2)]
            # graph.add_edges_from(graph_edges)
            # graph = gen_star_graph(6)
            gstate = GraphState(graph)
            for transmission in transmissions:
                placeholder_weigth = 0
                for wt in weigth:
                    decoder_Z = ErasureDecoder(gstate, "Z", p_fail, transmission, wt, in_qubit)
                    decoder_X = ErasureDecoder(gstate, "X", p_fail, transmission, wt, in_qubit)
                    succ_Z = decoder_Z.run_dec_updated()
                    succ_X = decoder_X.run_dec_updated()
                    average = (1-succ_X + 1-succ_Z) / 2
                    if average <= average_best:
                        average_best = average
                        placeholder_weigth = wt
                        print(average)
                if placeholder_weigth > 0:
                    wt = placeholder_weigth
                    t = np.linspace(0.8, 1, 100)
                    best_erasures = []
                    best_erasures_GHZ = []
                    for trans in t:
                        decoder_Z = ErasureDecoder(gstate, "Z", p_fail, trans, wt, in_qubit)
                        decoder_X = ErasureDecoder(gstate, "X", p_fail, trans, wt, in_qubit)
                        succ_Z = decoder_Z.run_dec_updated()
                        succ_X = decoder_X.run_dec_updated()
                        average = (1-succ_X + 1-succ_Z) / 2
                        best_erasures.append(average)
                        decoder_Z = ErasureDecoder(gstateGHZ, "Z", p_fail, trans, 0.001, in_qubit)
                        decoder_X = ErasureDecoder(gstateGHZ, "X", p_fail, trans, 0.001, in_qubit)
                        succ_Z = decoder_Z.run_dec_updated()
                        succ_X = decoder_X.run_dec_updated()
                        average = (1-succ_X + 1-succ_Z) / 2
                        best_erasures_GHZ.append(average)

                    plot_line_val = 0
                    for j in range(len(best_erasures_GHZ)):
                        if best_erasures_GHZ[j] >= best_erasures[j]:
                            plot_line_val = t[j]
                    fig1 = plt.figure()
                    gstate.image(input_qubits=[in_qubit])
                    plt.savefig(path_save + "\graph-" + str(i) + '.png')
                    plt.close(fig1)
                    fig1 = plt.figure()
                    plt.plot(t, best_erasures, label="graph-encoded")
                    plt.plot(t, best_erasures_GHZ, label="GHZ-encoded")
                    plt.axvline(x=plot_line_val, linestyle='--', color='r', label=str(plot_line_val)[:6]+"%")
                    best_erasures = np.array(best_erasures)
                    best_erasures_GHZ = np.array(best_erasures_GHZ)
                    plt.fill_between(t, best_erasures_GHZ, best_erasures, where=best_erasures_GHZ >= best_erasures, facecolor='lightgrey', interpolate=True)
                    plt.ylabel("$p_{erasure}$")
                    plt.xlabel("transmission")
                    plt.title("$p_{fail}$ = 1/2")
                    plt.legend()
                    plt.savefig(path_save + "\graph-plot-" + str(i) + '.png')
                    plt.close(fig1)
    
    '''
    '''
    weigths = []
    best_erasures = []
    loss_value = []
    for transmission in transmissions:
        #average_best = 1
        for wt in weigth:
            decoder_Z = ErasureDecoder(gstate, "Z", p_fail, transmission, wt, in_qubit)
            decoder_X = ErasureDecoder(gstate, "X", p_fail, transmission, wt, in_qubit)
            succ_Z = decoder_Z.run_dec()
            succ_X = decoder_X.run_dec()
            average = (1-succ_X + 1-succ_Z) / 2
            if average < average_best:
                average_best = average
                print("average {}".format(average))
                print("weigth {}".format(wt))
                print(succ_X, succ_Z)
                weigths.append(wt)
                best_erasures.append(average)
                loss_value.append(transmission)
    plt.plot(loss_value, weigths, label="weigths")
    plt.plot(loss_value, best_erasures, label="erasure_prob")
    plt.xlabel("Transmission")
    plt.legend()
    plt.show()
    '''

    '''
    in_qubit = 0
    ring_network_threshold = 0.1198
    p_fail = 1 / 4
    transmission = [0.98]
    weigth = 0.001
    arms = [3, 4, 5, 6, 7, 8, 9, 10]
    actual_arms = [i - 1 for i in arms]
    eras_plot_vals = []
    thresold_ring = []
    for i in arms:
        graph = gen_star_graph(i, 0)
        gstate = GraphState(graph)
        decoder = ErasureDecoder(gstate, "X", 1/4, 0.98, 0.001)
        print("my decoder {}".format(decoder.run_dec()))
        encoded_erasure = run_new_erase_decoder_zx(gstate, p_fail, transmission, weigth)[0]
        thresold_ring.append(ring_network_threshold)
        eras_plot_vals.append(encoded_erasure)
        print(encoded_erasure)
    plt.plot(actual_arms, eras_plot_vals, "-o", label="encoded-GHZ")
    plt.plot(actual_arms, thresold_ring, "k:", label="ring-network-threshold-value")
    plt.title("$p_{fail} =$" + " {}, t = {} % and weigth = {}".format(p_fail, transmission[0], weigth))
    plt.xlabel("Number of arms")
    plt.ylabel("Encoded erasure prob.")
    plt.legend()
    plt.show()
    '''
