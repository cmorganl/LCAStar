#!/usr/bin/python

from __future__ import division

try:
    import sys
    import traceback
    import re
    import itertools
    import operator
    import logging
    from math import log, pow, sqrt
    from python_resources.ncbiTaxonomyTree import *
except ImportError:
    sys.stderr.write(""" Could not load some user defined  module functions""" + "\n")
    sys.stderr.write(""" Make sure your typed \"source MetaPathwaysrc\"""" + "\n")
    sys.stderr.write(""" """ + "\n")
    sys.stderr.write(traceback.print_exc(10) + "\n")
    sys.exit(3)

try:
    from math import erf  # only available in Python 2.7
except ImportError:
    logging.warning("Could not import math.erf, will use interpolation for p-values." + "\n")


def copyList(a, b): 
    [ b.append(x) for x in a ] 


# Class to estimate and cumulative distribution functions
class CDF:

    def __init__(self, obs):
        self.obs = obs

    def __call__(self, x):
        counter = 0.0
        for obs in self.obs:
            if obs <= x:
                counter += 1
        return counter / len(self.obs)


class LCAStar(object):

    def __init__(self):

        # initialize class variables
        self.begin_pattern = re.compile("#")

        # a map from taxid to value, which has the S = sum n,  value for each taxid
        self.id_to_R = {}
        # a map from taxid to value, which has the S = sum n,  value for each taxid
        self.id_to_S = {}
        # a map from taxid to value, which has the L = sum n log n,  value for each taxid
        self.id_to_L = {}
        # a map from taxid to value, which has the entropy H value for each taxid
        self.id_to_H = {}
        # a map to keep track of visited nodes
        self.id_to_V = {}
        self.lca_min_score = 50   # an LCA parameter for min score for a hit to be considered
        self.lca_top_percent = 10  # an LCA param to confine the hits to within the top hits score upto the top_percent%
        self.lca_min_support = 5   # a minimum number of reads in the sample to consider a taxon to be present
        self.results_dictionary = None
        
        # LCAStar parameters
        self.lca_star_min_reads = 10
        self.lca_star_min_depth = 3
        self.lca_star_alpha = 0.51
        self.chi_squared_cdf = None

        # Standard NCBI Tree parameters
        self.ncbi_tree = None  # type: NcbiTaxonomyTree

    def setParameters(self, min_score, top_percent, min_support):
        self.lca_min_score = min_score
        self.lca_top_percent = top_percent
        self.lca_min_support = min_support
         
    def sizeTaxnames(self):
        return len(self.ncbi_tree.name_to_taxid)

    def sizeTaxids(self):
        return len(self.ncbi_tree.dic)
          
    def get_a_valid_taxid(self, name_group):
        for name in name_group:
            try:
                return self.ncbi_tree.name_to_taxid[name]
            except KeyError:
                return -1

    # given a taxon name it returns the corresponding unique ncbi tax taxid
    def translateNameToID(self, name):
        if name not in self.ncbi_tree.name_to_taxid:
            return None
        return self.ncbi_tree.name_to_taxid[name]

    # given a taxon taxid to taxon name map
    def translateIdToName(self, taxid: int):
        if taxid not in self.ncbi_tree.dic:
            return None
        return self.ncbi_tree.dic[taxid].name

    # given a name it returns the parents name
    def getParentName(self, name):
        if name not in self.ncbi_tree.name_to_taxid:
            return None
        taxid = self.ncbi_tree.name_to_taxid[name]
        pid = self.getParentTaxId(taxid)
        return self.translateIdToName(pid)

    # given a ncbi tax taxid returns the parents tax taxid
    def getParentTaxId(self, ID):
        if ID not in self.ncbi_tree.taxid_to_ptaxid:
            return None
        return self.ncbi_tree.taxid_to_ptaxid[ID][0]

    def get_lca(self, taxon_ids, return_id=False):
        """
        Given a set of ids it returns the lowest common ancestor, without caring about min support.
        Here LCA for a set of ids are computed as follows: 
        first we consider one ID at a time for each taxid we traverse up the NCBI tree using the taxid to parent taxid map
        at the same time increasing the count on the second value of the 3-tuple 
        note that at the node where all the of the individual ids (limit in number)
        converges the counter matches the limit for the first time, while climbing up. 
        This also this enables us to make the selection of taxid arbitrary 

        :param taxon_ids:
        :param return_id: 
        :return: 
        """
        limit = len(taxon_ids)
        for taxid in taxon_ids:
            tmp_id = int(taxid)
            while tmp_id in self.ncbi_tree.dic and tmp_id != 1:
                self.ncbi_tree.taxid_to_ptaxid[tmp_id][1] += 1
                if self.ncbi_tree.taxid_to_ptaxid[tmp_id][1] == limit:
                    if return_id:
                        return tmp_id
                    else:
                        return self.ncbi_tree.dic[tmp_id].name
                tmp_id = self.ncbi_tree.taxid_to_ptaxid[tmp_id][0]

        if return_id:
            return 1
        else:
            return "root"

    def update_taxon_support_count(self, taxonomy):
        taxid = self.get_a_valid_taxid([taxonomy])
        tid = taxid
        while tid in self.ncbi_tree.taxid_to_ptaxid and tid != 1:
            self.ncbi_tree.taxid_to_ptaxid[tid][2] += 1
            tid = self.ncbi_tree.taxid_to_ptaxid[tid][0]

    def get_supported_taxon(self, taxonomy):
        taxid = self.get_a_valid_taxid([taxonomy])
        tid = taxid
        while tid in self.ncbi_tree.taxid_to_ptaxid and tid !=1:
            if self.lca_min_support > self.ncbi_tree.taxid_to_ptaxid[tid][2]:
                tid = self.ncbi_tree.taxid_to_ptaxid[tid][0]
            else:
                return self.translateIdToName(tid)

        return self.translateIdToName(tid)

    def clear_cells(self, taxon_ids):
        """
        Need to call this to clear the counts of reads at every node

        :param taxon_ids:
        :return:
        """
        limit = len(taxon_ids)
        for taxid in taxon_ids:
            tid = taxid
            while tid in self.ncbi_tree.taxid_to_ptaxid and tid != 1:
                self.ncbi_tree.taxid_to_ptaxid[tid][1] = 0
                tid = self.ncbi_tree.taxid_to_ptaxid[tid][0]
        return ""

    def getTaxonomy(self, name_groups, return_id=False):
        """
        Given a set of sets of names it computes an lca in the format [ [name1, name2], [name3, name4,....namex] ...]
        here name1 and name2 are synonyms and so are name3 through namex

        :param name_groups:
        :param return_id:
        :return:
        """
        taxon_ids = []
        for name_group in name_groups:
            taxid = self.get_a_valid_taxid(name_group)
            if taxid != -1:
                taxon_ids.append(taxid)
        consensus = self.get_lca(taxon_ids, return_id)
        self.clear_cells(taxon_ids)
        return consensus

    # given an ID gets the lineage
    def get_lineage(self, taxid):
        tid = taxid
        lineage = []
        lineage.append(taxid)
        while tid in self.ncbi_tree.taxid_to_ptaxid and tid !=1:
            lineage.append(self.ncbi_tree.taxid_to_ptaxid[tid][0])
            tid = self.ncbi_tree.taxid_to_ptaxid[tid][0]
        return lineage
    
    # given two names calculates the distance on the tree
    def get_distance(self, taxa1, taxa2, debug=False):
        id1 = self.get_a_valid_taxid([taxa1])
        id2 = self.get_a_valid_taxid([taxa2])  # real
        lin1 = self.get_lineage(id1)
        lin2 = self.get_lineage(id2)
        if debug:
            logging.debug("id1:", str(id1) + "\n")
            logging.debug("id2:", str(id2) + "\n")
            logging.debug("lin1:", str(lin1) + "\n")
            logging.debug("lin2:", str(lin2) + "\n")
        large = None
        if len(lin1) <= len(lin2):
           large = lin2
           small = lin1
        else:
           large = lin1
           small = lin2
        for i in range(len(large)):
           for j in range(len(small)):
               if large[i] == small[j]:
                  return i + j
        return None

    # extracts taxon names for a refseq annotation
    def get_species(self, hit):
       if not 'product' in hit: 
           return None
       species = []
       try:
           m = re.findall(r'\[([^\[]+)\]', hit['product'])
           if m is not None:
             copyList(m, species)
       except:
             return None
   
       if species:
          return species
       else:
          return None
 
    # used for optimization
    def set_results_dictionary(self, results_dictionary):
        self.results_dictionary = results_dictionary

    def taxon_depth(self, taxon):
        taxid = self.translateNameToID(taxon)
        if taxid is None:
           return 0

        tid = taxid 
        depth = 0
        #climb up the tree from the taxon to the root 
        # the number of climbing steps is the depth
        while( tid in self.ncbi_tree.taxid_to_ptaxid and tid !=1 ):
            tid = self.ncbi_tree.taxid_to_ptaxid[tid][0]
            depth += 1

        return depth

    def filter_taxa_list(self, taxalist):
        # filter based on depth 
        newlist = []
        for taxon in taxalist:
            depth = self.taxon_depth(taxon)
            if depth < self.lca_star_min_depth:
                continue
            newlist.append(taxon)

        # filter based on min_reads /decide if we should return 'root'
        # or compute the taxon using lca_star
        if len(newlist) < self.lca_star_min_reads:
            return None 
        else:
            return newlist

    def setLCAStarParameters(self,  min_depth=3, alpha=0.51,  min_reads=5):
        self.lca_star_min_reads = min_reads
        self.lca_star_min_depth = min_depth
        self.lca_star_alpha = alpha
        
    def __read_counts(self, taxalist):
        read_counts = {}
        total = 0
        for taxon in taxalist:
            if taxon not in read_counts:
                read_counts[taxon] = 0
            read_counts[taxon] += 1
            total += 1
        return read_counts, total

    # find the taxon with the highest count but also has count higher than the 
    # majority threshold
    def lca_majority(self, taxalist):
        taxalist = self.filter_taxa_list(taxalist)

        majority = 'all'
        if taxalist is None:
            return majority

        majority = self.__lca_majority(taxalist)

        if majority is None:
            return 'all'

        return majority

    def __lca_majority(self, taxalist):
        # create the read counts
        read_counts, total = self.__read_counts(taxalist)
         
        # find the taxon with the highest count but also has count higher than the 
        # majority threshold
        majority = None
        maxcount = 0
        for taxon in taxalist: 
            if maxcount < read_counts[taxon] and total*self.lca_star_alpha < read_counts[taxon]:
                maxcount = read_counts[taxon]
                majority = taxon

        # majority exists 
        if majority is not None:
            return majority

        return None

    def __color_tree(self, read_counts):
        for taxon in read_counts:
            taxid = self.translateNameToID(taxon)
            tid = taxid
            # climb up the tree from the taxon to the root
            #  and mark the parent to child structure with True
            while tid in self.ncbi_tree.taxid_to_ptaxid and tid != 1:
                pid = self.ncbi_tree.taxid_to_ptaxid[tid][0]
                self.ncbi_tree.ptaxid_to_taxid[pid][tid] = True
                tid = pid

    def __annotate_tree_counts(self, read_counts):
        for taxon in read_counts:
            value = read_counts[taxon]
            taxid = self.translateNameToID(taxon)
            if taxid is None:
                continue
            self.id_to_R[taxid] = value

    def __decolor_tree(self):
        S = [1]
        while len(S) > 0:
            taxid = S.pop()

            C = []
            if taxid in self.ncbi_tree.ptaxid_to_taxid:
              C = self.ncbi_tree.ptaxid_to_taxid[taxid].keys()
           
            for child in C:
               if self.ncbi_tree.ptaxid_to_taxid[taxid][child]:
                  self.ncbi_tree.ptaxid_to_taxid[taxid][child] = False
                  S.append(child)

    def __create_majority(self, root, read_name_counts):
        read_counts = {}
        total = 0
        for taxon in read_name_counts:
            count = read_name_counts[taxon]
            taxid = self.translateNameToID(taxon)
            read_counts[taxid] = count
            total += count

        candidate = [1, 10000000.00]
        Stack = [root]
        while len(Stack) > 0:
            taxid = Stack.pop()
            # calculate here
            if taxid in self.id_to_V:
                C = []
                if taxid in self.ncbi_tree.ptaxid_to_taxid:
                    C = self.ncbi_tree.ptaxid_to_taxid[taxid].keys()
                # I am coming up
                self.id_to_H[taxid] = 0

                if taxid in read_counts:
                    self.id_to_S[taxid] = float(read_counts[taxid])
                    self.id_to_L[taxid] = float(read_counts[taxid])*log(float(read_counts[taxid]))
                else:
                    self.id_to_S[taxid] = 0
                    self.id_to_L[taxid] = 0

                for child in C:
                    if self.ncbi_tree.ptaxid_to_taxid[taxid][child]:
                        self.id_to_S[taxid] += self.id_to_S[child]
                        self.id_to_L[taxid] += self.id_to_L[child]
                try:
                    self.id_to_H[taxid] = -(self.id_to_L[taxid]/self.id_to_S[taxid] - log(self.id_to_S[taxid]))
                except:
                    logging.debug("ID: " + str(taxid) + "\n")
                    exit(-1)
                if self.id_to_S[taxid] > total*self.lca_star_alpha:
                    if candidate[1] > self.id_to_H[taxid]:
                        candidate[0] = taxid
                        candidate[1] = self.id_to_H[taxid]
            else:  # going down
                self.id_to_V[taxid] = True
                C = []
                if taxid in self.ncbi_tree.ptaxid_to_taxid:
                    C = self.ncbi_tree.ptaxid_to_taxid[taxid].keys()

                Stack.append(taxid)
                for child in C:
                    if self.ncbi_tree.ptaxid_to_taxid[taxid][child]:
                        Stack.append(child)
        
        return candidate[0]

    def __clear_lca_star_data_structure(self):
        self.id_to_R = {}
        self.id_to_S = {}
        self.id_to_L = {}
        self.id_to_H = {}
        self.id_to_V = {}

    # monotonically decreasing function of depth of divergence d
    def step_cost(self, d):
        return 1 / pow(2, d)
    
    def wtd_distance(self, obs, exp, debug=False):
        """ weighted taxonomic distance calculates the distance between 
            observed and expected positions on the NCBI Taxonomy Database
            weighting each step by a cost function step_cost base in the 
            depth in the tree.
        """
        exp_id = self.get_a_valid_taxid([exp])
        obs_id = self.get_a_valid_taxid([obs])
        exp_lin = self.get_lineage(exp_id)
        obs_lin = self.get_lineage(obs_id)
        if debug:
            logging.debug("Expected: ", exp + "\n")
            logging.debug("Observed: ", obs + "\n")
            logging.debug("Expected ID: ", str(exp_id) + "\n")
            logging.debug("Observed ID: ", str(obs_id) + "\n")
            logging.debug("Expected Lineage:", str(exp_lin) + "\n")
            logging.debug("Observed Lineage :", str(obs_lin) + "\n")
        sign = -1

        # check to see if expected in observed lineage
        # if so distance sign is positive
        if exp_id in obs_lin:
            sign = 1
        large = None
        if len(obs_lin) <= len(exp_lin):
            # expected longer than observed
            large = exp_lin
            small = obs_lin
        else:
            large = obs_lin
            small = exp_lin

        # calculate cost
        a_cost = 0
        b_cost = 0
        for i in range(len(large)):
            if i > 0:
                a_cost += self.step_cost(len(large)-i-1)
            b_cost = 0
            for j in range(len(small)):
                if j > 0:
                    b_cost += self.step_cost(len(small)-j-1)
                if large[i] == small[j]:
                    return (a_cost + b_cost) * sign
        return None  # did not find lineages

    # Function takes an LCAStar candidate and read_counts object and returns a revised
    # taxonomy distribution to calculate a p-value
    def __collapse_distribution(self, result_id, read_counts):
        new_taxa_dist_ids = []
        for r in read_counts:
            r_id = self.get_a_valid_taxid([r])
            lin = self.get_lineage(r_id)
            if result_id in lin:
                # found in lineage so replace with candidate
                r_id = result_id
            # create a list of taxa
            for i in range(read_counts[r]):
                new_taxa_dist_ids.append(r_id)
        new_data_dist = list(map(self.translateIdToName, map(int, new_taxa_dist_ids)))
        return new_data_dist

    def lca_star(self, taxa_list, return_id=False):
        # filter taxa dist by depth
        taxa_list = self.filter_taxa_list(taxa_list)
        
        if taxa_list is None:
            return 'all', None
        
        majority = self.__lca_majority(taxa_list)
        
        if majority is not None:
            p_val = self.calculate_pvalue(taxa_list, majority)
            if return_id:
                majority = self.translateNameToID(majority)
            return majority, str(p_val)
        
        read_counts, total = self.__read_counts(taxa_list)
        
        # Calculate LCA Star
        self.__annotate_tree_counts(read_counts)
        self.__color_tree(read_counts)
        result_id = self.__create_majority(1, read_counts)
        collapsed_taxa_list = self.__collapse_distribution(result_id, read_counts)
        self.__clear_lca_star_data_structure()
        result_taxon = self.translateIdToName(result_id)
        p_val = self.calculate_pvalue(collapsed_taxa_list, result_taxon)
        self.__decolor_tree()
        if return_id:
            return result_id, str(p_val)
        return result_taxon, str(p_val)

    # Chi-squared with one degree of freedom
    def chi_squared(self, x):
        return erf(sqrt(x/2))

    # Calculate pvalue based on Nettleton result
    def calculate_pvalue(self, taxa_list, taxa):
        if len(list(taxa_list)) <= 1:
            logging.warning("p-value not defined for taxa lists of <2\n")
            return None
        elif taxa not in taxa_list:
            logging.warning("Tried to calculate a p-value for a taxa not in taxa_list:\n" + str(taxa_list) + "\n")
            return None

        X = {}  # hash of taxa counts
        M = 0  # maximum count

        for t in taxa_list:
            if t not in X:
                X[t] = 0
            X[t] += 1
        try:
            X_k = X[taxa]  # taxa count for test statistic
        except KeyError:
            logging.error("Unable to properly parse taxa_list into a dictionary." + "\n")
            sys.exit()
        for t in X:
            if t != taxa:
                if X[t] > M:
                    M = X[t]

        T = 0  # test statistic

        if X_k <= M:
            # trivial case
            return round(1 - self.chi_squared(T), 3)
        else:
            first = 0
            if M > 0:
                first = M * log((2 * M) / (M + X_k))
            second = X_k * log((2 * X_k) / (M + X_k))
            T = 2 * (first + second)
            return round(1 - self.chi_squared(T), 3)

    # Returns the simple majority taxa: calculate the most common taxa in a given list
    def simple_majority(self, L, return_id=False):
        if len(L) == 0:
            return "root"
        # get an iterable of (item, iterable) pairs
        SL = sorted((x, i) for i, x in enumerate(L))
        # print 'SL:', SL
        groups = itertools.groupby(SL, key=operator.itemgetter(0))
        # auxiliary function to get "quality" for an item

        def _auxfun(g):
            item, iterable = g
            count = 0
            min_index = len(L)
            for _, where in iterable:
                count += 1
                min_index = min(min_index, where)
            # print 'item %r, count %r, minind %r' % (item, count, min_index)
            return count, -min_index
        # pick the highest-count/earliest item
        majority = max(groups, key=_auxfun)[0]
        if return_id:
            majority = self.get_a_valid_taxid([majority])
        return majority
