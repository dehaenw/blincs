# BLINCS
Breadth-first Line Notation for Chemical Structures​.

Here is the initial implementation of BLINCS. It is reasonably stable, depends on just RDKit and the python stdlib. A presentation with more info is available here (link will be inserted here in ~1 day)

## Installation

clone repo, go to dir and run `pip install .`

## Usage
Convert a rdkit mol to blincs
```python
from rdkit import Chem
import blincs
mol = Chen.MolFromSmiles("C1CCCC1")
blincs.mol_to_blincs
```
convert blincs to an rdkit mol
```python
blincs.blincs_to_mol("C1CCCC1")
```
## Known issues
Stereochemistry on adamantanes and similar molecules is sometimes inconsistent. Stereochemistry parsing needs to be more efficient and will probably changed to a SMILES-like system.
