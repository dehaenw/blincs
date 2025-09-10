from rdkit import Chem
from collections import Counter
import re
import math
from rdkit import RDLogger  
RDLogger.DisableLog('rdApp.*')       

pt = Chem.GetPeriodicTable()
bonddict = {"-":Chem.BondType.SINGLE,"=":Chem.BondType.DOUBLE,"#":Chem.BondType.TRIPLE,":":Chem.BondType.AROMATIC}
tokenizer = re.compile(r"[=#:-]?(?:[EZ])?(?:%\d+|\[.*?\]\d*|(?:Cl|Br|[A-Z])\d*|[cnbosp]\d*)")

VERBOSE = False


def seq_is_ok(seq):
    """
    check if ring openings are closed and if edges and degrees match.
    """
    ok = True
    rings = [num for num in seq if num < 0]
    degrees = [num if num > 0 else 0 for num in seq ]
    parity = sum([0 if c%2==0 else 1 for c in dict(Counter(rings)).values()])
    if parity > 0:
        ok = False
    if sum(degrees)/2 != len(seq)-len(rings)/2-1:
        ok = False
    return ok
    

def seq_to_adj(seq):
    """
    this function takes an sequence of vertex degrees and turns it into
    connectivity information (an array of index pairs). For example:
    `[3,2,1,1,2,1]`
    turns into
    `[(0,1),(0,2),(0,3),(1,4),(4,5)]`
    the integers correspond to the degree of vertex, if the degree is n,
    the n-1 (or n for the first, parentless element) next elements are the 
    children of that vertex. This can describe any tree.
    To broaden this approach for graphs, ring opening and closing operations
    can be added, analogously to those used in SMILES notation. These numbers
    are notated using a negative sign. For example:
    `1,3,-1,2,2,-1`
    corresponds to
    `[(0,1),(1,3),(1,2),(2,3)]`
    """
    adj = []
    if not seq_is_ok(seq):
        return adj
    ring_dict = {}
    n_atoms = sum([num>0 for num in seq])
    curr_idx_parent = 0
    curr_idx_child = 0
    curr_atom = 0
    curr_free_atom = 0
    while curr_idx_parent < len(seq) and curr_atom<n_atoms:
        while seq[curr_idx_parent] < 0:
            curr_idx_parent += 1
        deg = seq[curr_idx_parent]
        parentless = 1 if curr_atom == 0 else 0
        for i in range(deg-1+parentless):
            curr_idx_child += 1
            if seq[curr_idx_child] > 0:
                curr_free_atom += 1
                adj.append((curr_atom,curr_free_atom))
            else:
                if seq[curr_idx_child] in ring_dict:
                    ring_dict[seq[curr_idx_child]].append(curr_atom)
                else:
                    ring_dict[seq[curr_idx_child]] = [curr_atom]
        curr_atom += 1
        curr_idx_parent += 1
    for k in ring_dict:
        ring_atoms = ring_dict[k]
        for i in range(len(ring_atoms)//2):
            adj.append((ring_atoms[2*i],ring_atoms[2*i+1]))
    return adj

def blincs_parser(bli):
    """
    takes a blincs string and separate it into bond, atom and degree tokens
    """
    tokens = re.findall(tokenizer, bli)
    if bli != "".join(tokens):
        if VERBOSE:
            print("tokenization failure")
            print(tokens)
        return
    seq = []
    tokens_info = []
    for token in tokens:
        info = {"charge":0,"isotope":0,"atom":0,"parentbond":"-","parentbondstereo":"","expliciths":0,"aromatic":False,"chiral":"","noimplicit":False}
        if "%" in token:
            if token[0] in "-#:=":
                info["parentbond"] = token[0]
                if token[1] in "EZ":
                    info["parentbondstereo"] = token[1]
            seq.append(-1*int(token[token.index("%")+1:]))
        else:
            if token[-1] in "1234567890":
                decimals = 1
                while token[-decimals] in "1234567890":
                    decimals+=1
                seq.append(int(token[-decimals+1:]))
                token = token[:-decimals+1]
            else:
                seq.append(2)
            if token[0] in "-#:=":
                info["parentbond"] = token[0]
                token = token[1:]
            if token[0] in "EZ":
                info["parentbondstereo"] = token[0]
                token = token[1:]
            if token[0] == "[":
                info["noimplicit"] = True
                token = token[1:-1]
                if token[0] in "1234567890":
                    decimals = 0
                    while token[decimals] in "1234567890":
                        decimals+=1
                    info["isotope"] = int(token[:decimals])
                    token = token[decimals:]
                if "-" in token or "+" in token:
                    try:
                        idx = token.index("-")
                    except:
                        idx = token.index("+")
                    charge = token[idx:]
                    if len(charge) == 1:
                        charge+="1"
                    info["charge"] = int(charge)
                    token = token[:idx]
                if "H" in token[1:]:
                    idx = token[1:].index("H")+1
                    hcount = token[idx:]
                    if len(hcount) == 1:
                        hcount+="1"
                    info["expliciths"] = int(hcount[1:])
                    token = token[:idx]
                else:
                    info["expliciths"] = 0
                if token[-1]=="@":
                    info["chiral"]="@"
                    token = token[:-1]
                    if token[-1]=="@":
                        info["chiral"]="@@"
                        token = token[:-1] 
            if token[0].islower():
                info['aromatic'] = True
            info['atom'] = pt.GetAtomicNumber(token.capitalize())     
        tokens_info.append(info)
    return(seq,tokens_info)

def create_mol(seq,tokens_info):
    mol = Chem.RWMol()
    n_atoms = [num for num in seq if num>=0]
    atom_info = [info for i,info in enumerate(tokens_info) if seq[i]>=0]
    ring_info = []
    rings_info = [(seq[i],info) for i,info in enumerate(tokens_info) if seq[i]<0]
    open_ring = []
    for ring_seq,info in rings_info:
        if ring_seq not in open_ring:
            open_ring.append(ring_seq)
            ring_info.append(info)
        else:
            open_ring.pop(open_ring.index(ring_seq))
    correct_cip_tags = {}
    chiral_double_bonds = []
    for i,atom in enumerate(n_atoms):
        at = Chem.Atom(0)
        if atom_info[i]["aromatic"]:
            at.SetIsAromatic(True)
        if atom_info[i]["noimplicit"]:
            at.SetNoImplicit(True)
        at.SetFormalCharge(atom_info[i]["charge"])
        at.SetIsotope(atom_info[i]["isotope"])
        at.SetAtomicNum(atom_info[i]["atom"])
        at.SetNumExplicitHs(atom_info[i]["expliciths"])
        if atom_info[i]["chiral"]:
            if atom_info[i]["chiral"] == "@":
                at.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
                correct_cip_tags[i] = "R"
            else:
                at.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
                correct_cip_tags[i] = "S"
        mol.AddAtom(at)
    pairs = seq_to_adj(seq)
    for i,pair in enumerate(pairs[:len(n_atoms)-1]):
        if mol.GetAtomWithIdx(pair[0]).GetIsAromatic() and mol.GetAtomWithIdx(pair[1]).GetIsAromatic():
            mol.AddBond(*pair,bonddict[":"])
        else:
            mol.AddBond(*pair,bonddict[atom_info[pair[1]]["parentbond"]])
            if atom_info[pair[1]]["parentbondstereo"]:
                chiral_double_bonds.append([pair[0],pair[1],atom_info[pair[1]]["parentbondstereo"]])
    for i,pair in enumerate(pairs[len(n_atoms)-1:]):
        if mol.GetAtomWithIdx(pair[0]).GetIsAromatic() and mol.GetAtomWithIdx(pair[1]).GetIsAromatic():
            mol.AddBond(*pair,bonddict[":"])
        else:
            mol.AddBond(*pair,bonddict[ring_info[i]["parentbond"]])
            if atom_info[pair[1]]["parentbondstereo"]:
                chiral_double_bonds.append([pair[0],pair[1],atom_info[pair[1]]["parentbondstereo"]])
    try:
        Chem.SanitizeMol(mol)
    except:
        if VERBOSE:
            print("sanitize fail")
        mol = None
    try:
        Chem.AssignCIPLabels(mol)
        Chem.AssignStereochemistry(mol, flagPossibleStereoCenters=True)
        ciplabs = []
        for a in mol.GetAtoms():
            try:
                ciplabs.append(a.GetProp("_CIPCode"))
            except:
                pass
        for i,a in enumerate(mol.GetAtoms()):
            if i in correct_cip_tags:
                try:
                    if a.GetProp("_CIPCode").lower() != correct_cip_tags[i].lower(): #lowercase due to rR and sE
                        if a.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
                            a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
                        else:
                            a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
                except Exception as e:
                    pass#print("CipCode issue",e)
            if "r" in ciplabs or "s" in ciplabs:
                Chem.AssignCIPLabels(mol) #must be done every step to fix pseudoassymetry
    except Exception as e:
        if VERBOSE:
            print("chirality correction fail",e)
        mol=None
    try:
        for a1,a2,label in chiral_double_bonds:
            a1neighbors = [at for at in mol.GetAtomWithIdx(a1).GetNeighbors() if at.GetIdx()!=a2]
            a2neighbors = [at for at in mol.GetAtomWithIdx(a2).GetNeighbors() if at.GetIdx()!=a1]
            a1highcip = sorted(a1neighbors,key=lambda x:-int(x.GetProp("_CIPRank")))[0].GetIdx()
            a2highcip = sorted(a2neighbors,key=lambda x:-int(x.GetProp("_CIPRank")))[0].GetIdx()
            mol.GetBondBetweenAtoms(a1,a2).SetStereoAtoms(a1highcip,a2highcip) 
            if label=="E":
                mol.GetBondBetweenAtoms(a1,a2).SetStereo(Chem.BondStereo.STEREOE)
                mol.GetBondBetweenAtoms(a1,a1highcip).SetBondDir(Chem.BondDir.ENDUPRIGHT)
                mol.GetBondBetweenAtoms(a2,a2highcip).SetBondDir(Chem.BondDir.ENDUPRIGHT)
            else:
                mol.GetBondBetweenAtoms(a1,a2).SetStereo(Chem.BondStereo.STEREOZ)
                mol.GetBondBetweenAtoms(a1,a1highcip).SetBondDir(Chem.BondDir.ENDUPRIGHT)
                mol.GetBondBetweenAtoms(a2,a2highcip).SetBondDir(Chem.BondDir.ENDDOWNRIGHT)
    except Exception as e:
        if VERBOSE:
            print("bond stereochemistry fail",e)
    return mol

def get_token(a):
    sym = pt.GetElementSymbol(a.GetAtomicNum())
    cha = a.GetFormalCharge()
    iso = a.GetIsotope()
    aro = a.GetIsAromatic()
    hco = a.GetNumExplicitHs()
    chi = str(a.GetChiralTag())
    nre = a.GetNumRadicalElectrons()
    try:
        cla = a.GetProp("_CIPCode")
    except:
        cla = ""
    if cha == 0 and iso == 0 and sym in ["B","C","N","O","F","P","S","Cl","Br","I"] and hco == 0 and cla == "" and nre==0:
        token = sym
        if aro:
            token = token.lower()
    else:
        if aro:
            sym = sym.lower()
        token = sym
        if iso:
            token = str(iso)+token
        if cla == "R":
            token += "@"
        elif cla == "r":
            token += "@"
        elif cla == "S":
            token += "@@"
        elif cla == "s":
            token += "@@"
        if nre > 0:
            hco = a.GetTotalNumHs()
        if hco:
            if hco==1:
                token += "H"
            else:
                token += f"H{hco}"
        else:
            if nre > 0:
                token += "H0"
        if cha:
            if cha>0:
                if cha == 1:
                    token += "+"
                else:
                    token += f"+{cha}"
            else:
                if cha == -1:
                    token += "-"
                else:
                    token += str(cha)
        token = f"[{token}]"
    return token


def adj_to_seq(adj,tokens=None,bond_types={}):
    """
    takes all edges of a graph as the input, encoded as a tuple of edges
    for example:
    `[(0,1),(1,3),(3,0),(3,4)]`
    and outputs a degree sequence with ring operations. for the above:
    `[2,2,3,-1,-1,1]`
    """
    unique_idx = sorted(list(set([p[0] for p in adj] + [p[1] for p in adj])))
    node_visited = {idx:False for idx in unique_idx}
    idx_neighbors = {idx:[[cidx for cidx in p if cidx!=idx][0] for p in adj if idx in p] for idx in unique_idx}
    idx_degree = {idx:str(len(idx_neighbors[idx])) if len(idx_neighbors[idx]) != 2 else "" for idx in unique_idx}
    if not tokens:
        tokens = {idx:"" for idx in unique_idx}
    idx_parent = {idx:None for idx in unique_idx}
    expand_queue = [unique_idx[0]]
    idx_parent[unique_idx[0]] = -1 # root node
    seq = [tokens[unique_idx[0]]+str(idx_degree[unique_idx[0]])]
    ring_to_idx = {}
    ring_counter = 0
    while len(expand_queue)>0:
        curr_idx = expand_queue.pop(0)
        curr_parent = idx_parent[curr_idx]
        curr_neighbors = [idx for idx in idx_neighbors[curr_idx] if idx > 0 and idx != curr_parent]
        if node_visited[curr_idx]:
            if VERBOSE:
                print("this shouldnt happen")
        node_visited[curr_idx] = True
        for curr_neighbor in curr_neighbors:
            bond = tuple(sorted([curr_idx,curr_neighbor]))
            if bond in bond_types:
                bs = bond_types[bond]
            else:
                bs = ""
            if node_visited[curr_neighbor] or curr_neighbor in expand_queue:
                if bond in ring_to_idx:
                    seq.append(ring_to_idx[bond])
                else:
                    ring_counter+=1
                    ring_to_idx[bond]=f"{bs}%{ring_counter}"
                    seq.append(f"{bs}%{ring_counter}")
            else:
                idx_parent[curr_neighbor] = curr_idx 
                expand_queue.append(curr_neighbor)
                seq.append(bs+tokens[curr_neighbor]+idx_degree[curr_neighbor])
    return "".join(seq)

def mol_to_blincs(mol):
    blis = []
    try:
        for frag in Chem.GetMolFrags(mol,asMols=True): #in case of salts etc:
            Chem.AssignCIPLabels(frag)
            frag = Chem.RemoveHs(Chem.AddHs(frag))
            bonds = []
            idx_to_token = {}
            bond_to_symbol = {}
            for a in frag.GetAtoms():
                idx = a.GetIdx()
                idx_to_token[idx] = get_token(a)
                for n in a.GetNeighbors():
                    bond = tuple(sorted([idx,n.GetIdx()]))
                    if bond not in bonds:
                        bonds.append(bond)
            for b in frag.GetBonds():
                bt = str(b.GetBondType())
                bs= ""
                if bt not in ["AROMATIC","SINGLE"]:
                    if bt == "DOUBLE":
                        bs = "="
                        if str(b.GetStereo()) in ["STEREOE","STEREOTRANS"]:
                            bs = "=E"
                        elif str(b.GetStereo()) in ["STEREOZ","STEREOCIS"]:
                            bs = "=Z"
                    elif bt == "TRIPLE":
                        bs = "#"
                    if bs!="":
                        bond = tuple(sorted([b.GetBeginAtomIdx(),b.GetEndAtomIdx()]))
                        bond_to_symbol[bond] = bs
            if len(bonds)>0:
                bli = adj_to_seq(bonds,idx_to_token,bond_to_symbol)
            else:
                bli = [idx_to_token[k] for k in idx_to_token][0]+ "0"
            blis.append(bli)
    except Exception as e:
        if VERBOSE:
            print("mol to blincs failure",e)
        blis = [""]
    return ".".join(sorted(blis))

def blincs_to_mol(bli):
    try:
        blis = bli.split(".")
        frags = []
        for cbli in blis:
            s,m = blincs_parser(cbli)
            frags.append(create_mol(s,m))
        if len(frags) == 1:
            mol = frags[0]
        else:
            m1 = frags.pop(0)
            m2 = frags.pop(0)
            mol = Chem.CombineMols(m1,m2)
            while len(frags)>0:
                mol = Chem.CombineMols(mol,frags.pop(0))
    except Exception as e:
        if VERBOSE:
            print("blincs parse failure",e)
        mol = None
    return mol


