import copy
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem, Draw

from rdkit.Chem import rdMolTransforms
from rdkit.Chem.Draw import IPythonConsole

IPythonConsole.drawOptions.addAtomIndices = True
IPythonConsole.molSize = 300, 300


class Environment:
    """
    The Environment class is a wrapper for the RDKit molecule object. It provides methods for generating a 3D
    conformer, enumerating the torsions, and evaluating a set of parameters by generating a new conformation and
    calculating its energy. Environment contains information about its mol: molecule, mol_3D: 3D conformer,
    and torsion_atoms: list of lists, where each list contains the atom indices of the four atoms defining a torsion.

    :param smiles:  SMILES representation of the molecule
    :type smiles:  str
    """
    def __init__(self, smiles):
        """
        Constructor method
        """
        self.mol = Chem.MolFromSmiles(smiles)  # RDKit molecule object
        self.mol_3D = self.generate_conformer()
        self.torsion_atoms = self.enumerate_torsions()

    def draw(self):
        """
        Draw the molecule with the torsion atoms highlighted.

        :return:  a grid of molecular images with the torsion atoms highlighted
        :rtype:  PIL object
        """
        highlight_atoms = [i[1:3] for i in self.torsion_atoms]
        mols_to_draw = [Chem.AddHs(self.mol)]*len(self.torsion_atoms)
        return Draw.MolsToGridImage(mols_to_draw, highlightAtomLists=highlight_atoms)

    def enumerate_torsions(self):
        """
        Enumerate all the torsions in the molecule  (i.e., all the rotatable bonds).

        :return:  list of lists, where each list contains the atom indices of the four atoms defining a torsion
        :rtype:  list
        """
        torsion_smarts = '[!$(*#*)&!D1]-!@[!$(*#*)&!D1]'

        # Implemented from the https://sourceforge.net/p/rdkit/mailman/message/34554615/ with modifications.
        def remove_dihedrals_describing_same_conformation(listoflists):
            """
            Remove dihedrals that describe the same conformation.

            :param listoflists:  list of lists, where each list contains indices of four atoms defining a torsion
            :type listoflists:  list
            :return: list of lists, where each list contains indices of four atoms defining a unique torsion
            :rtype:  list
            """
            final_dihedrals = []
            for item in listoflists:
                if len(final_dihedrals) == 0:
                    final_dihedrals.append(item)
                else:
                    flag = 0
                    for j in final_dihedrals:
                        if len(set(item[1:3]) & set(j[1:3])) == 2:
                            flag = 1
                    if not flag:
                        final_dihedrals.append(item)
            return final_dihedrals

        matches = self.mol_3D.GetSubstructMatches(Chem.MolFromSmarts(torsion_smarts))
        torsion_list = []
        for match in matches:
            idx2, idx3 = match[0], match[1]
            bond = self.mol_3D.GetBondBetweenAtoms(idx2, idx3)
            j_atom, k_atom = self.mol_3D.GetAtomWithIdx(idx2), self.mol_3D.GetAtomWithIdx(idx3)
            if (((j_atom.GetHybridization() != Chem.HybridizationType.SP2)
                 and (j_atom.GetHybridization() != Chem.HybridizationType.SP3))
                    or ((k_atom.GetHybridization() != Chem.HybridizationType.SP2)
                        and (k_atom.GetHybridization() != Chem.HybridizationType.SP3))):
                continue
            for b1 in j_atom.GetBonds():
                if b1.GetIdx() == bond.GetIdx():
                    continue
                idx1 = b1.GetOtherAtomIdx(idx2)
                for b2 in k_atom.GetBonds():
                    if ((b2.GetIdx() == bond.GetIdx())
                            or (b2.GetIdx() == b1.GetIdx())):
                        continue
                    idx4 = b2.GetOtherAtomIdx(idx3)
                    # skip 3-membered rings
                    if idx4 == idx1:
                        continue
                    atoms = [idx1, idx2, idx3, idx4]
                    torsion_list.append(atoms)
        return remove_dihedrals_describing_same_conformation(torsion_list)

    def generate_conformer(self):
        """
        Generate a 3D conformer of the molecule.

        :return:  3D RDKit molecule object with a single conformer
        :rtype:  RDKit molecule object
        """
        m = copy.deepcopy(self.mol)
        m_h = Chem.AddHs(m)
        AllChem.EmbedMolecule(m_h)
        new_mol = Chem.Mol(m_h)
        new_mol.RemoveAllConformers()
        new_mol.AddConformer(m_h.GetConformer(0))
        return new_mol

    @staticmethod
    def get_conformer_energy(conformer):
        """
        Get the energy of a conformer using the MMFF force field.

        :param conformer:  conformation of the molecule
        :type conformer:  RDKit molecule object
        :return:  energy of the conformer
        :rtype:  float
        """
        mp = AllChem.MMFFGetMoleculeProperties(conformer)
        ff = AllChem.MMFFGetMoleculeForceField(conformer, mp)
        return ff.CalcEnergy()

    def create_root_parameters(self):
        """
        Create a set of random parameters for the torsions in the molecule.

        :return:  list of random parameters for the torsions
        :rtype:  list
        """
        parameters = np.random.uniform(0, 360, len(self.torsion_atoms))
        return parameters

    def mol_from_parameters(self, mol, parameters):
        """
        Generate a new molecule with a given set of torsion angles.

        :param mol:  RDKit molecule object to be modified
        :type mol:  RDKit molecule object
        :param parameters:  list of torsion angles to be applied to the molecule
        :type parameters:  list
        :return:  new conformation created by applying the torsion angles
        :rtype:  RDKit molecule object
        """
        for n, indices in enumerate(self.torsion_atoms):
            at1, at2, at3, at4 = indices
            rdMolTransforms.SetDihedralDeg(mol.GetConformer(0), at1, at2, at3, at4, parameters[n])
        return mol

    def evaluate(self, parameters):
        """
        Evaluate a set of parameters by generating a new conformation and calculating its energy.

        :param parameters:  list of torsion angles to be applied to the molecule
        :type parameters:  list
        :return:  list of parameters and the energy of the new conformation
        :rtype:  list
        """
        new_conformation = self.mol_from_parameters(self.mol_3D, parameters)
        return [parameters, self.get_conformer_energy(new_conformation)]
