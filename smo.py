import sys, getopt
from glpk import *
import islpy as isl

#from enum import Enum
#import impipe
#import matplotlib.pyplot as plt
#from matplotlib.ticker import MultipleLocator
#from mpl_toolkits.mplot3d import Axes3D
#plt.switch_backend("QT4Agg")

C = 8

global DEBUG
DEBUG = False

global GLPKSOLVE
GLPKSOLVE = False

def dprint(obj):
    if DEBUG == True:
        print obj
    else:
        return

    

########################################################
# Name management
########################################################
class NameType(object):
    param = 0
    input = 1
    output = 2
    farkas_multiplier = 3
    decision = 4
    hyperplane_coeff = 5
    misc = 6

def generate_names(name, start, n):
    names = []
    for i in range(n):
        names.append(name + str(start + i))
    return names

def get_names(n, name_type):
    if name_type == NameType.param:
        names = generate_names('N', get_names.counter[NameType.param], n)
    elif name_type == NameType.input: 
        names = generate_names('x', get_names.counter[NameType.input], n)
    elif name_type == NameType.output: 
        names = generate_names('y', get_names.counter[NameType.output], n)
    elif name_type == NameType.farkas_multiplier: 
        names = generate_names('L', get_names.counter[NameType.farkas_multiplier], n)
    elif name_type == NameType.decision: 
        names = generate_names('d', get_names.counter[NameType.decision], n)
    elif name_type == NameType.hyperplane_coeff: 
        names = generate_names('c', get_names.counter[NameType.hyperplane_coeff], n)
    elif name_type == NameType.misc: 
        names = generate_names('u', get_names.counter[NameType.misc], n)
    get_names.counter[name_type] += n
    return names
get_names.counter = [0, 0, 0, 0, 0, 0, 0]

def analyze_constraint(constraint):
    dprint(constraint)

# skips comments and empty lines
def read_next_line(f):
    line = f.readline()
    while line[0] == '#' or line[0] == '\n':
        line = f.readline()
    return line

def read_smo_input(ctx, inputfile):
    f = open(inputfile, 'r')

    parameters = get_names(1, NameType.param)
    # read the constraints on the array space
    domain = isl.BasicSet(read_next_line(f))
    domain = domain.reset_space(isl.Space.create_from_names(ctx, set=get_names(len(domain.get_space().get_var_ids(isl._isl.dim_type.set)), NameType.input), params=parameters)).get_basic_sets()[0]

    dprint("DOMAIN is")
    dprint(domain.lexmin())

    image = isl.BasicSet(read_next_line(f))
    image = image.reset_space(isl.Space.create_from_names(ctx, set=get_names(len(image.get_space().get_var_ids(isl._isl.dim_type.set)), NameType.output), params=parameters)).get_basic_sets()[0]
    
    # create an order pair given the domain and range
    ordered_pair = isl.BasicMap.from_domain_and_range(domain, image)

    # read the number of conflict polyhedra in the conflict set
    n_conflict_polyhedra = int(read_next_line(f))
    
    # !!! Ideally, the order pair should be constructed using a range which has thse output names
    output_names = read_next_line(f).split()
    # image = domain.reset_space(isl.Space.create_from_names(ctx, set=output_names))

    cs = []
    for i in range(n_conflict_polyhedra):
        polyhedron = isl.BasicMap(read_next_line(f), ctx)
	polyhedron = polyhedron.add_constraints(ordered_pair.get_constraints())
	cs.append(polyhedron)
    
    f.close()
    return cs

class Smo:
    def apply_farkas_lemma(self, polyhedron, cst):

        def gather_coefficients(cst_spec, cst_set, dim_type, polyhedron):
            constraints = []
            faces = polyhedron.get_constraints()
            for i in range(len(faces)):
                # get the coeffcient of a particular in_ dim in a given face
                for var_index in range(len(polyhedron.get_space().get_var_ids(dim_type))):
                    name = polyhedron.get_space().get_var_ids(dim_type)[var_index]
                    if cst_spec.get(name) != None:
                        cst_spec[name].update({multipliers[i] : -faces[i].get_coefficient_val(dim_type, var_index).get_num_si()})
                    else:
                        cst_spec[name] = {multipliers[i] : -faces[i].get_coefficient_val(dim_type, var_index).get_num_si()}

            eq_faces = filter(lambda x : x.is_equality(), faces)
            for i in range(len(eq_faces)):
                # get the coeffcient of a particular in_ dim in a given face
                for var_index in range(len(polyhedron.get_space().get_var_ids(dim_type))):
                    name = polyhedron.get_space().get_var_ids(dim_type)[var_index]
                    if cst_spec.get(name) != None:
                        cst_spec[name].update({multipliers[i + len(faces)] : eq_faces[i].get_coefficient_val(dim_type, var_index).get_num_si()})
                    else:
                        cst_spec[name] = {multipliers[i + len(faces)] : eq_faces[i].get_coefficient_val(dim_type, var_index).get_num_si()}

            for name in polyhedron.get_space().get_var_ids(dim_type):
                if cst_spec.get(name) != None:
                    constraint = isl.Constraint.eq_from_names(cst_set.get_space(), cst_spec[name])
                    cst_set = cst_set.add_constraint(constraint)
                    dprint("Adding constraint due to " + name)
                    dprint(constraint)
                    dprint("\n")
                    dprint(cst_set)
            return cst_spec, cst_set

        faces = polyhedron.get_constraints()
        dprint("=== Faces ===")
        dprint(faces)
        dprint("\n")
        self.nvariables = len(self.space.get_var_ids(isl._isl.dim_type.set))

        nmultipliers = len(faces) + len(filter(lambda x : x.is_equality(), faces)) + 1
        cst_set = isl.BasicSet.universe(self.space)
        cst_set = cst_set.add_dims(isl._isl.dim_type.set, nmultipliers)
        multipliers = get_names(nmultipliers, NameType.farkas_multiplier)

        # get a farkas multiplier for each face
        for i in range(len(faces)):
            cst_set = cst_set.set_dim_name(isl._isl.dim_type.set, self.nvariables + i, multipliers[i])
            cst_set = cst_set.add_constraint(isl.Constraint.ineq_from_names(cst_set.get_space(), {multipliers[i] : 1}))

        # get another farkas multiplier for each 'equality' face
        eq_faces = len(filter(lambda x : x.is_equality(), faces))
        for i in range(len(filter(lambda x : x.is_equality(), faces))):
            cst_set = cst_set.set_dim_name(isl._isl.dim_type.set, self.nvariables + i + len(faces), multipliers[len(faces) + i])
            cst_set = cst_set.add_constraint(isl.Constraint.ineq_from_names(cst_set.get_space(), {multipliers[len(faces) + i] : 1}))

        # get another farkas multiplier for the constant part in the affine combination of the faces
        cst_set = cst_set.set_dim_name(isl._isl.dim_type.set, self.nvariables + nmultipliers - 1, multipliers[nmultipliers - 1])
        cst_set = cst_set.add_constraint(isl.Constraint.ineq_from_names(cst_set.get_space(), {multipliers[nmultipliers - 1] : 1}))
        
        cst, cst_set = gather_coefficients(cst, cst_set, isl._isl.dim_type.param, polyhedron)
        cst, cst_set = gather_coefficients(cst, cst_set, isl._isl.dim_type.set, polyhedron)
        cst, cst_set = gather_coefficients(cst, cst_set, isl._isl.dim_type.in_, polyhedron)

        # constant values due the faces
        for i in range(len(faces)):
            if cst.get(1) != None:
                cst[1].update({multipliers[i] : -faces[i].get_constant_val().get_num_si()})
            else:
                cst[1] = {multipliers[i] : -faces[i].get_constant_val().get_num_si()}


        # constant values due the eq_faces
        eq_faces = filter(lambda x : x.is_equality(), faces)
        for i in range(len(eq_faces)):
            if cst.get(1) != None:
                cst[1].update({multipliers[i + len(faces)] : eq_faces[i].get_constant_val().get_num_si()})
            else:
                cst[1] = {multipliers[i + len(faces)] : eq_faces[i].get_constant_val().get_num_si()}

        # constant values due the additional farkas multiplier
        cst[1].update({multipliers[nmultipliers - 1] : -1})

        if cst.get(1) != None:
            constraint = isl.Constraint.eq_from_names(cst_set.get_space(), cst[1])
            dprint("Adding constraint due to constant")
            dprint(constraint)
            dprint("\n")
            cst_set = cst_set.add_constraint(constraint)

        dprint("====== Constraints =======")
        dprint(cst_set)
        dprint("\n")

        return cst_set

    def eliminate(self, polyhedron, cst_set):
        faces = polyhedron.get_constraints()
        eq_faces = filter(lambda x : x.is_equality(), faces)

        #dprint("======= elimination ===========")

        for i in range(len(faces) + len(eq_faces) + 1):
            #dprint("\n")
            cst_set = cst_set.project_out(isl._isl.dim_type.set, self.nvariables, 1)
            #dprint(cst_set)
            #dprint("\n")

        """
        cst_set = cst_set.project_out(isl._isl.dim_type.set, self.nvariables, len(faces) + len(eq_faces) + 1)
        #cst_set = cst_set.remove_dims(isl._isl.dim_type.set, self.nvariables, len(faces) + len(eq_faces) + 1)
        """

        dprint("======= After elimination ===========")
        dprint(cst_set)
        dprint("\n")
        return cst_set

    # (u.P + w) <= cP + c
    def approximate_upper_bound_constraint(self, polyhedron):
        dprint("APPROXIMATE UPPER BOUND")
        cst = {1 : {1 : C, self.ubound[-1] : -1}}
        for i in range(len(self.ubound) - 1):
            cst[polyhedron.get_space().get_dim_name(isl._isl.dim_type.param, i)] = { self.ubound[i] : -1, 1 : C}
        cst_set = self.apply_farkas_lemma(polyhedron, cst)
        cst_set = self.eliminate(polyhedron, cst_set)
        return cst_set


    # flag ? (h.s - h.t) <= u.P + w :  -(h.s - h.t) <= u.P + w
    def get_upper_bound_constraints(self, polyhedron, flag):
        dprint("UPPER BOUND CONSTRAINTS")
        cst = {1 : {self.ubound[-1] : 1}}
        for i in range(len(self.ubound) - 1):
            cst[polyhedron.get_space().get_dim_name(isl._isl.dim_type.param, i)] = { self.ubound[i] : 1}

        for i in range(len(polyhedron.get_space().get_var_ids(isl._isl.dim_type.set))):
            cst[polyhedron.get_space().get_dim_name(isl._isl.dim_type.set, i)] = { self.coefficients[i] : -1 if flag else 1}
            cst[polyhedron.get_space().get_dim_name(isl._isl.dim_type.in_, i)] = { self.coefficients[i] : 1 if flag else -1}

        cst_set = self.apply_farkas_lemma(polyhedron, cst)
        cst_set = self.eliminate(polyhedron, cst_set)
        return cst_set

    # (decision_index % 2 == 0) : (h.s - h.t) >= 1 - (1 - d1)(cP + c + 1) : (h.s - h.t) <= -1 + (1 - d2)(cP + c + 1)
    def get_decision_constraints(self, polyhedron, decision_index):
        dprint("DECISION CONSTRAINTS")
        cst = {1 : {1 : -1 + C + 1, self.bdvs[decision_index] : - C - 1}}
        for i in range(len(self.ubound) - 1):
            name = polyhedron.get_space().get_dim_name(isl._isl.dim_type.param, i)
            cst[name] = { self.bdvs[decision_index] : -C, 1 : C}
            #cst[1].update({ name : C})

        for i in range(len(polyhedron.get_space().get_var_ids(isl._isl.dim_type.set))):
            cst[polyhedron.get_space().get_dim_name(isl._isl.dim_type.set, i)] = { self.coefficients[i] : 1 if (decision_index % 2 == 0) else -1}
            cst[polyhedron.get_space().get_dim_name(isl._isl.dim_type.in_, i)] = { self.coefficients[i] : -1 if (decision_index % 2 == 0) else 1}

        cst_set = self.apply_farkas_lemma(polyhedron, cst)
        cst_set = self.eliminate(polyhedron, cst_set)
        return cst_set

    def add_constraints(self, cst_set, new_constraints):
        #dprint("Before adding constraint")
        #dprint(cst_set)
        #dprint(new_constraints)
        cst_set = cst_set.intersect(new_constraints)
        #dprint("After adding constraint")
        #dprint(cst_set)
        #if cst_set.is_empty():
        #    dprint("EMPTY SET RESULTED!")
        return cst_set

    def formulate_bounding_constraints(self, polyhedron, index):
        #dprint("====== Conflict Polyhedron ========")
        #dprint(polyhedron)
        #dprint("\n")
        cst_set = isl.BasicSet.universe(self.space)
        cst_set = self.add_constraints(cst_set, self.approximate_upper_bound_constraint(polyhedron))
        cst_set = self.add_constraints(cst_set, self.get_upper_bound_constraints(polyhedron, True))
        cst_set = self.add_constraints(cst_set, self.get_upper_bound_constraints(polyhedron, False))
        cst_set = self.add_constraints(cst_set, self.get_decision_constraints(polyhedron, 2*index))
        cst_set = self.add_constraints(cst_set, self.get_decision_constraints(polyhedron, 2*index+1))

        cst_set = cst_set.add_constraint(isl.Constraint.ineq_from_names(cst_set.get_space(), {self.bdvs[2*index] : 1}))
        cst_set = cst_set.add_constraint(isl.Constraint.ineq_from_names(cst_set.get_space(), {self.bdvs[2*index] : -1, 1 : 1}))
        cst_set = cst_set.add_constraint(isl.Constraint.ineq_from_names(cst_set.get_space(), {self.bdvs[2*index+1] : 1}))
        cst_set = cst_set.add_constraint(isl.Constraint.ineq_from_names(cst_set.get_space(), {self.bdvs[2*index+1] : -1, 1 : 1}))

        return cst_set

    def get_eta_bound(self, cs):
        cst = {1 : len(cs), self.eta[0] : -1}
        for i in range(len(cs)):
            cst[self.bdvs[2*i]] = -1
            cst[self.bdvs[2*i+1]] = -1
        return cst

    def revise_conflict_set_deprecated(self, cs, lexmin):
        basic_set = lexmin.get_basic_sets()[0]
        hyperplane = {}
        cst = {}
        input_names = cs[0].get_space().get_var_names(isl._isl.dim_type.set)
        output_names = cs[0].get_space().get_var_names(isl._isl.dim_type.in_)
        for i in range(len(self.coefficients)):
            hyperplane[self.coefficients[-1-i]] = -basic_set.get_constraints()[i].get_constant_val().get_num_si()
            cst[input_names[-1-i]] = hyperplane[self.coefficients[-1-i]]
            cst[output_names[-1-i]] = -hyperplane[self.coefficients[-1-i]]

        constraint = isl.Constraint.eq_from_names(cs[0].get_space(), cst)
        print constraint

        cs_revised = []
        for i in range(len(cs)):
            cs[i] = cs[i].add_constraint(constraint)
            if cs[i].is_empty() == False:
                cs_revised.append(cs[i])
        return cs_revised

    def revise_conflict_set(self, cs, coeff_values):
        cst = {}
        input_names = cs[0].get_space().get_var_names(isl._isl.dim_type.set)
        output_names = cs[0].get_space().get_var_names(isl._isl.dim_type.in_)
        for i in range(len(self.coefficients)):
            cst[input_names[i]] = coeff_values[i]
            cst[output_names[i]] = -coeff_values[i]

        constraint = isl.Constraint.eq_from_names(cs[0].get_space(), cst)
        print constraint

        cs_revised = []
        for i in range(len(cs)):
            cs[i] = cs[i].add_constraint(constraint)
            if cs[i].is_empty() == False:
                cs_revised.append(cs[i])
        return cs_revised

    def optimize(self, inputfile):
        def is_empty(cs):
            for i in range(len(cs)):
                if cs[i].is_empty() == False:
                    return False
            return True

        def get_storage_hyperplane(cs, lexmin):
            def get_value(basic_set, var_index):
                for constraint in basic_set.get_constraints():
                    if constraint.get_coefficient_val(isl._isl.dim_type.set, var_index) == 1:
                        return constraint.get_constant_val().get_num_si()
                print "Error! Value unassigned!"
                self.exit(2)

            lexmin_set = lexmin.get_basic_sets()[0]
            coeff_values = []
            for i in range(len(self.coefficients)):
                coeff_values.append(-get_value(lexmin_set, len(lexmin_set.get_var_ids(isl._isl.dim_type.set))-i-1))
            coeff_values.reverse()

            max_diff_values = []
            for i in range(len(self.params) + 1):
                max_diff_values.append(-get_value(lexmin_set, 1+i))

            return coeff_values, max_diff_values

        def print_access(names, params, coeff_max_diff_tuple):
            coeff_values, max_diff_values = coeff_max_diff_tuple
            print("("),
            for i in range(len(coeff_values)):
                if abs(coeff_values[i]) == 1:
                    print("+ " + names[i] if coeff_values[i] == 1 else "- " + names[i]),
                elif coeff_values[i] != 0:
                    print(("+ " if coeff_values[i] > 0 else "") + str(coeff_values[i]) + names[i]),
            print(") mod ("),

            print("" if max_diff_values[-1] == -1 else (max_diff_values[-1] + 1)),
            for i in range(len(max_diff_values) - 1):
                if abs(max_diff_values[i]) == 1:
                    print(("+ " + params[i]) if max_diff_values[i] == 1 else ("- " + params[i])),
                elif max_diff_values[i] != 0:
                    print(("+ " if max_diff_values[i] > 0 else "") + str(max_diff_values[i]) + params[i]),
            print(")"),

        def glpk_solve(cst_set, names, params):
            size = 1000+1
            ia = intArray(size)
            ja = intArray(size)
            ar = doubleArray(size)
            prob = glp_create_prob()

            glp_set_prob_name(prob, "problem")
            glp_set_obj_dir(prob, GLP_MIN)

            glp_add_rows(prob, len(cst_set.get_constraints()))
            for i in range(len(cst_set.get_constraints())):
                print cst_set.get_constraints()[i]
                glp_set_row_name(prob, i+1, "s"+str(i+1))
                if cst_set.get_constraints()[i].is_equality():
                    glp_set_row_bnds(prob, i+1, GLP_FX, -cst_set.get_constraints()[i].get_constant_val().get_num_si(), 0.0)
                else:
                    glp_set_row_bnds(prob, i+1, GLP_LO, -cst_set.get_constraints()[i].get_constant_val().get_num_si(), 0.0)

            ndims = len(cst_set.get_space().get_var_ids(isl._isl.dim_type.set))
            glp_add_cols(prob, ndims)
            for i in range(ndims):
                glp_set_col_name(prob, i+1, cst_set.get_space().get_var_ids(isl._isl.dim_type.set)[i])
                glp_set_col_kind(prob, i+1, GLP_IV)
                glp_set_col_bnds(prob, i+1, GLP_FR, 0.0, 0.0)

            glp_set_obj_coef(prob, 1, 100.0)
            glp_set_obj_coef(prob, 2, 10.0)
            glp_set_obj_coef(prob, 3, .1)

            count = 1
            for i in range(len(cst_set.get_constraints())):
                for j in range(ndims):
                    value = cst_set.get_constraints()[i].get_coefficient_val(isl._isl.dim_type.set, j).get_num_si()
                    if value != 0:
                        ia[count] = i+1
                        ja[count] = j+1
                        ar[count] = value
                        count = count+1

            glp_load_matrix(prob, count-1, ia, ja, ar)
            glp_write_lp(prob, None, "./lp.lp")
            problem = glpk("./lp.lp")
            problem.solve()
            solution = problem.solution()

            coeff_values, max_diff_values = [], []
            for i in range(ndims):
                name = cst_set.get_space().get_var_ids(isl._isl.dim_type.set)[i]
                print name + " = " + str(solution[name])
                if i > 0 and i <= len(params) + 1:
                    max_diff_values.append(solution[name])
                if i >= (ndims - len(names)):
                    coeff_values.append(solution[name])
            del prob
            return coeff_values, max_diff_values 

        self.ctx = isl.Context()
        cs = read_smo_input(self.ctx, inputfile)
        dim_names = cs[0].get_space().get_var_ids(isl._isl.dim_type.set)
        params = cs[0].get_space().get_var_ids(isl._isl.dim_type.param)

        hyperplanes = []
        while is_empty(cs) == False:
            dprint(cs)
            self.bdvs = get_names(2*len(cs), NameType.decision)
            self.ubound = get_names(len(cs[0].get_space().get_var_ids(isl._isl.dim_type.param)) + 1, NameType.misc)
            self.eta = get_names(1, NameType.misc)
            self.coefficients = get_names(len(cs[0].domain().get_space().get_var_ids(isl._isl.dim_type.set)), NameType.hyperplane_coeff)
            self.params = cs[0].get_space().get_var_names(isl._isl.dim_type.param)
            self.space = isl.Space.create_from_names(self.ctx, set = self.eta + self.ubound + self.bdvs + self.coefficients, params=self.params)

            cst_set = isl.BasicSet.universe(self.space)
            for i in range(len(cs)):
                dprint("============ TREATING A NEW CONFLICT POLYHEDRA ===========")
                bounding_constraints = self.formulate_bounding_constraints(cs[i], i)
                dprint(bounding_constraints)
                cst_set = cst_set.intersect(bounding_constraints)

            cst_set = cst_set.add_constraint(isl.Constraint.eq_from_names(cst_set.get_space(), self.get_eta_bound(cs)))

            # setting a few bounds on the hyperplane coefficients

            for i in range(len(self.coefficients)):
                cst_set = cst_set.add_constraint(isl.Constraint.ineq_from_names(cst_set.get_space(), {self.coefficients[i]: 1, 1: 4}))
                cst_set = cst_set.add_constraint(isl.Constraint.ineq_from_names(cst_set.get_space(), {self.coefficients[i]: -1, 1 : 4}))

            dprint(cst_set)

            #f = open('./smo.output', 'a')
	    #cst_set.print_(f, 1, '', '', 1)
	    #f.close()
            dprint("============= <LEXMIN> ==============")
            if GLPKSOLVE == True:
                cst_set = cst_set.remove_divs()
                coeff_values, max_diff_values = glpk_solve(cst_set, dim_names, params)
            else:
                lexmin = cst_set.lexmin()
                dprint(lexmin)
                coeff_values, max_diff_values = get_storage_hyperplane(cs, lexmin)
            print("Hyperplane, Modulo ="),
            print coeff_values, max_diff_values
            hyperplanes.append((coeff_values, max_diff_values))
            
            dprint("============= <LEXMIN\> ==============")
            cs = self.revise_conflict_set(cs, coeff_values)

        print "\nStorage mapping found"
        print(dim_names),
        print("---> ["),
        for i in range(len(hyperplanes)):
            print_access(dim_names, params, hyperplanes[i])
            if i < len(hyperplanes) - 1:
                print(","),
        print("]")

        
def main(argv):
    inputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hdgi:",["ifile=", "debug", "glpk-solve"])
        for opt, arg in opts:
            if opt == '-h':
               print'smo.py -i <inputfile>\n'
               print '--------------- OPTIONS ----------------'
               print '-g/--glpk-solve\t\t\tUse glpk solver'
               print '-d/--debug\t\t\tdebug print'
               print '-i/--ifile\t\t\tconflict specification input file'
               sys.exit()
            elif opt in ("-i", "--ifile"):
               inputfile = arg
               Smo().optimize(inputfile)
               sys.exit()
            elif opt in ("-d", "--debug"):
	       global DEBUG
               DEBUG = True
            elif opt in ("-g", "--glpk-solve"):
	       global  GLPKSOLVE
               GLPKSOLVE = True
    except getopt.GetoptError:
        print 'smo.py -i <inputfile>'
        sys.exit(2)
    print 'smo.py -i <inputfile>'
    sys.exit(2)


#islApiExamples()
main(sys.argv[1:])
