

from suffix_tree import SuffixTree
from suffix_tree import GeneralisedSuffixTree



class IgorSuffixTree(object):

    def __init__(self, string):
        self.__stree = SuffixTree(string)


    def search_stree(self, string):
        current_index = -1
        terminal = False
        string = list(string)
        str_index = 0
        current = self.__stree.root
        while True:
            current_child = current.firstChild
            found = False
            while True:
                if current_child == None:
                    break
                edgeLabel = current_child.edgeLabel
                #print edgeLabel
                if string[str_index] == edgeLabel[0]:
                    found = True
                    break
                else:
                    #elif current_child.edgeLabel == current_child.lastChild.edgeLabel:
                    #    break
                    current_child = current_child.next
            if found:
                #print edgeLabel
                len_matched = 0
                #print edgeLabel
                current_index = current_child.index
                while len_matched < len(edgeLabel) and str_index < len(string) and string[str_index] == edgeLabel[len_matched]:
                    str_index += 1
                    len_matched += 1
                    #print str_index, len_matched
                if len_matched < len(edgeLabel):
                    if edgeLabel[len_matched] in ["|", "$"]:
                        terminal = True
                    break
                else:
                    if str_index == len(string):
                        terminal = True
                        break
                    current = current_child
            else:
                break
        return str_index, current_index, terminal




